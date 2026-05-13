#include "../inc/ISD.h"
#include "../inc/quantum_utils.h"
#include "../inc/grover.h"
#include <iostream>
#include <map>
#include <algorithm>

int ISD::get_combinations(int n, int k) {
    if (k > n) return 1;
    if (k == 0 || k == n) return 1;
    if (k > n / 2) k = n - k;

    int res = 1;
    for (int i = 1; i <= k; ++i) {
        res = res * (n - i + 1) / i;
    }
    return res;
}

Operator ISD::controlled_solve_from_back(int n, int k, int start_index) {
    Operator solve;
    int x_start = start_index + (n + 1) * k;
    
    // Shift needed to target the bottom-right submatrix
    int row_shift = k - n; 

    // NEW: Calculate the safe stopping point.
    // If n=4, k=3, start_col = 1. If n=2, k=3, start_col = 0.
    int start_col = std::max(0, n - k);

    // FIX 1: Iterate down to start_col instead of 0 to prevent p_row < 0
    for (int i = n - 1; i >= start_col; i--) { 
        int p_row = i + row_shift; // Calculate correct pivot row
        Operator pivot_search_and_row_reduce;
        
        // pivot search
        for (int j = p_row - 1; j >= 0; j--) {
            Operator add_row_op;
            // Note: 'it' loops should also safely stay within bounds of valid columns if needed,
            // but since they rely on 'i' they will now stay safe.
            for (int it = i - 1; it >= 0; it--)
                add_row_op.add_gate<CNOTMatrix>({get_qubit_index(j, it, k, start_index), get_qubit_index(p_row, it, k, start_index)});
            add_row_op.add_gate<CNOTMatrix>({get_qubit_index(j, n, k, start_index), get_qubit_index(p_row, n, k, start_index)});

            for (int it = p_row; it > j; it--)
                pivot_search_and_row_reduce.add_gate<NotMatrix>({get_qubit_index(it, i, k, start_index)});

            vector<int> target_cnot_qubits;
            for (int it = p_row; it > j; it--)
                target_cnot_qubits.push_back(get_qubit_index(it, i, k, start_index));

            for (auto const &gate: add_row_op.gates) {
                Matrix gate_with_CNOTs = create_multiple_controlled_matrix(gate.matrix, p_row - j);
                vector<int> target_qubits = target_cnot_qubits;
                for (auto const &index: gate.target_indices)
                    target_qubits.push_back(index);
                
                pivot_search_and_row_reduce.add_gate({gate_with_CNOTs, target_qubits});
            }

            for (int it = p_row; it > j; it--)
                pivot_search_and_row_reduce.add_gate<NotMatrix>({get_qubit_index(it, i, k, start_index)});
        }

        // row reduce
        for (int j = p_row - 1; j >= 0; j--) {
            Operator add_row_op;
            for (int it = i - 1; it >= 0; it--)
                add_row_op.add_gate<CNOTMatrix>({get_qubit_index(p_row, it, k, start_index), get_qubit_index(j, it, k, start_index)});
            add_row_op.add_gate<CNOTMatrix>({get_qubit_index(p_row, n, k, start_index), get_qubit_index(j, n, k, start_index)});

            for (auto const &gate: add_row_op.gates)             
                pivot_search_and_row_reduce.add_gate(create_controlled_gate(gate, get_qubit_index(j, i, k, start_index)));
        }

        for (const auto &gate: pivot_search_and_row_reduce.gates) 
            solve.add_gate(create_controlled_gate(gate, x_start + i));
    }

    // back substitution 
    // FIX 2: Start from start_col so back substitution only hits valid, pivoted rows
    for (int i = start_col; i < n; i++) { 
        int p_row = i + row_shift;
        Operator back_substitution;
        // Target rows strictly below the pivot
        for (int j = p_row + 1; j < k; j++)
            back_substitution.add_gate<CCNOTMatrix>({get_qubit_index(j, i, k, start_index), get_qubit_index(p_row, n, k, start_index), get_qubit_index(j, n, k, start_index)});

        for (const auto &gate: back_substitution.gates) 
            solve.add_gate(create_controlled_gate(gate, x_start + i));
    }

    return solve;
}


Operator ISD::get_uniform_superposition_operator(int n, int k, int ancilla_start_index) {
    Operator create_uniform_superposition_op;
    
    int c_bits = bit_width((unsigned int)(k));
    create_uniform_superposition_op.add_operator(QuantumUtils::get_NOTs(c_bits, k, ancilla_start_index, 1));

    int ancilla_after_c_index = ancilla_start_index + c_bits;
    Operator check_c_ones = QuantumUtils::check_if_all_ones(c_bits, ancilla_start_index, ancilla_after_c_index);
    Operator reverse_check_c_ones = reverse_operator(check_c_ones);

    Operator c_minus_one = QuantumUtils::get_minus_one_operator(c_bits, ancilla_start_index);

    for (int i = 0; i < n; i++) {
        for (int j = k; j >= 1; j--) {
            Operator nots_to_j = QuantumUtils::get_NOTs(c_bits, j, ancilla_start_index, 0);            
            create_uniform_superposition_op.add_operator(nots_to_j);

            create_uniform_superposition_op.add_operator(check_c_ones);

            double one_probability = double(get_combinations(n - i - 1, j - 1)) / get_combinations(n - i, j);
            double theta = 2.0 * asin(sqrt(one_probability));
            create_uniform_superposition_op.add_gate(create_controlled_gate({RotationMatrix(theta), {i}}, ancilla_after_c_index + c_bits - 2));

            create_uniform_superposition_op.add_operator(reverse_check_c_ones);
            create_uniform_superposition_op.add_operator(reverse_operator(nots_to_j));
        }

        for (const auto &gate: c_minus_one.gates)
            create_uniform_superposition_op.add_gate(create_controlled_gate(gate, i));
    }

    return create_uniform_superposition_op;
}


int ISD::get_qubit_index(int row, int column, int k, int start_index){
    return start_index + column * k + row;
}


Operator ISD::get_optimisation_operator(int n, int k, int ancilla_start_index) {
    Operator optimisation;

    int c_bits = bit_width((unsigned int)(n - k));
    optimisation.add_operator(QuantumUtils::get_NOTs(c_bits, n - k, ancilla_start_index, 1));

    int ancilla_after_c_index = ancilla_start_index + c_bits;

    Operator check_c_ones = QuantumUtils::check_if_all_ones(c_bits, ancilla_start_index, ancilla_after_c_index);
    Operator reverse_check_c_ones = reverse_operator(check_c_ones);
    Operator c_minus_one = QuantumUtils::get_minus_one_operator(c_bits, ancilla_start_index);

    int x_start = ancilla_start_index - k;

    for (int i = 0; i < n - k; i++) {
        for (int j = 0; j < i; j++) {
            Operator nots_to_n_minus_k_minus_j = QuantumUtils::get_NOTs(c_bits, n - k - j, ancilla_start_index, 0);
            Operator reverse_nots = reverse_operator(nots_to_n_minus_k_minus_j);

            optimisation.add_operator(nots_to_n_minus_k_minus_j);
            optimisation.add_operator(check_c_ones);

            Operator swap_rows;
            for (int it = 0; it < k + 1; it++)
                swap_rows.add_gate<SwapMatrix>({get_qubit_index(i, it, n - k, n + 1), get_qubit_index(j, it, n - k, n + 1)});

            for (const auto &gate: swap_rows.gates) 
                optimisation.add_gate(create_controlled_gate(create_controlled_gate(gate, ancilla_after_c_index + c_bits - 2), i));

            optimisation.add_operator(reverse_check_c_ones);
            optimisation.add_operator(reverse_nots);
        }

        for (const auto &gate: c_minus_one.gates)
            optimisation.add_gate(create_controlled_gate(gate, i));
    }

    for (int i = k - 1; i >= 0; i--) {
        Operator op;
        op.add_operator(c_minus_one);
        op.add_gate<NotMatrix>({x_start + i});

        for (int j = i; j < k - 1; j++) {
            for (int it = 0; it < n - k; it++)
                op.add_gate<SwapMatrix>({get_qubit_index(it, j, n - k, n + 1), get_qubit_index(it, j + 1, n - k, n + 1)});

            op.add_gate<SwapMatrix>({x_start + j, x_start + j + 1});
        }

        for (const auto &gate: op.gates)
            optimisation.add_gate(create_controlled_gate(gate, n - k + i));
    }

    return optimisation;
}

void ISD::run_superposition_solver(ISDInfo isd_info, int it = 1000) {
    cout << "Running superposition solver..." << endl;;
    cout << "Expected probability of weight one: " << double(isd_info.good_columns_choices.size()) / get_combinations(isd_info.n, isd_info.n - isd_info.k) * 100 << "%" << endl << endl;

    int matrix_size = (isd_info.k + 1) * (isd_info.n - isd_info.k);
    int c_bits = bit_width((unsigned int)(isd_info.n - isd_info.k)); 
    int circuit_size = isd_info.n + 1 + matrix_size + isd_info.k + c_bits + c_bits - 1; 

    int ancilla_start = isd_info.n + 1 + (isd_info.k + 1) * (isd_info.n - isd_info.k) + isd_info.k;

    Circuit superposition_solver(circuit_size);
    superposition_solver.add_operator(get_uniform_superposition_operator(isd_info.n, isd_info.n - isd_info.k, ancilla_start));
    superposition_solver.add_operator(isd_info.matrix_init);
    superposition_solver.add_operator(get_optimisation_operator(isd_info.n, isd_info.k, ancilla_start));
    superposition_solver.add_operator(controlled_solve_from_back(isd_info.k, isd_info.n - isd_info.k, isd_info.n + 1));
    superposition_solver.run();

    int syndrome_start = isd_info.n + 1 + isd_info.k * (isd_info.n - isd_info.k);
    map<int, int> M;
    for (int i = 0; i < it; i++) {
        vector<bool> read = superposition_solver.magic_read();

        int res = 0;
        for (int j = 0; j < isd_info.n - isd_info.k; j++)
            res += (read[syndrome_start + j]);

        M[res]++;
    }

    for (auto &[res, cnt]: M) 
        cout << res << ": " << double(cnt) / it * 100 << "%\n";
}

ISDInfo ISD::get_small_isd_info() {
    /*
    1 0 0 1 1 | 1
    0 1 0 1 1 | 1
    0 0 1 0 1 | 0
    */

    int n = 5, k = 2;
    int matrix_start = n + 1;

    Operator init;
    init.add_gate<NotMatrix>({matrix_start + 0});
    init.add_gate<NotMatrix>({matrix_start + 1});
    init.add_gate<NotMatrix>({matrix_start + 3});
    init.add_gate<NotMatrix>({matrix_start + 4});
    init.add_gate<NotMatrix>({matrix_start + 5});
    init.add_gate<NotMatrix>({matrix_start + 6});
    init.add_gate<NotMatrix>({matrix_start + 7});


    vector<int> good_columns_subsets {13, 14, 25, 26};
    return {n, k, good_columns_subsets, init};
}

ISDInfo ISD::get_isd_info() {
    /*
    1 0 0 0 1 1 | 0
    0 1 0 1 0 1 | 0
    0 0 1 1 1 0 | 1
    */

    cout << "Performing ISD for H|s:" << endl;
    cout << "1 0 0 1 0 1 | 1" << endl;
    cout << "0 1 0 1 1 0 | 1" << endl;
    cout << "0 0 1 0 1 1 | 0" << endl;

    int n = 6, k = 3;
    int matrix_start = n + 1;
    Operator init;
    init.add_gate<NotMatrix>({matrix_start + 0});
    init.add_gate<NotMatrix>({matrix_start + 1});
    init.add_gate<NotMatrix>({matrix_start + 4});
    init.add_gate<NotMatrix>({matrix_start + 5});
    init.add_gate<NotMatrix>({matrix_start + 6});
    init.add_gate<NotMatrix>({matrix_start + 8});
    init.add_gate<NotMatrix>({matrix_start + 9});
    init.add_gate<NotMatrix>({matrix_start + 10});
    
    vector<int> good_columns_subsets {13, 14, 25, 26, 28, 41, 42, 44};
    return {n, k, good_columns_subsets, init};
}

GroverInfo ISD::get_grover_info(ISDInfo isd_info) {
    cout << "Running amplitude amplification for weight one..." << endl;
    cout << "Expected high probabilities of: ";
    for (auto &column: isd_info.good_columns_choices)
        cout << column << " ";
    cout << endl << endl;

    auto &n = isd_info.n;
    auto &k = isd_info.k;

    int c_bits = bit_width((unsigned int)(n - k));
    int syndrome_start = n + 1 + k * (n - k);
    int ancilla_start = n + 1 + (k + 1) * (n - k) + k;

    Operator isd_function;
    isd_function.add_operator(get_optimisation_operator(n, k, ancilla_start));
    isd_function.add_operator(controlled_solve_from_back(k, n - k, n + 1));

    Operator grover_function = grover.get_Fn_for_ISD(n, isd_function, syndrome_start, ancilla_start + c_bits);
    Operator superposition_op = get_uniform_superposition_operator(n, n - k, ancilla_start);

    return {n, (int)(isd_info.good_columns_choices.size()), ancilla_start + c_bits + c_bits - 1 - n - 1, get_combinations(n, n - k), grover_function, superposition_op, isd_info.matrix_init};
}

ISDInfo ISD::get_big_isd_info() {
    /*
    1 0 0 1 1 1 0 | 0
    0 1 0 0 1 1 1 | 0
    0 0 1 1 0 1 1 | 1
    */

    int n = 7, k = 4;
    int matrix_start = n + 1;
    
    Operator init;
    init.add_gate<NotMatrix>({n + 0});
    init.add_gate<NotMatrix>({n + 2});
    init.add_gate<NotMatrix>({n + 3});
    init.add_gate<NotMatrix>({n + 4});
    init.add_gate<NotMatrix>({n + 6});
    init.add_gate<NotMatrix>({n + 7});
    init.add_gate<NotMatrix>({n + 8});
    init.add_gate<NotMatrix>({n + 10});
    init.add_gate<NotMatrix>({n + 11});
    init.add_gate<NotMatrix>({n + 14});

    return {n, k, {13, 14, 25, 26, 28, 41, 42, 44, 73, 74, 88, 104}, init};
}


Operator ISD::sus(int dim, int start_index, int all_columns ) {
    int x_start = start_index + (dim * (all_columns + 1));

    int offset = max(all_columns - dim, 0);

    int start_column = max(0, dim - all_columns);

    Operator solve;
    for (int i = 0 + start_column; i < dim - 1; i++) {
        Operator pivor_search_and_row_reduce_op;
        //pivot search
        for (int j = i + 1; j < dim; j++) {            
            Operator add_row_op;
            for (int it = i + 1; it < dim + 1; it++)
                add_row_op.add_gate<CNOTMatrix>({get_qubit_index(j, it + offset - start_column, dim, start_index), get_qubit_index(i, it + offset - start_column, dim, start_index)});

            for (int it = i; it < j; it++)
                pivor_search_and_row_reduce_op.add_gate<NotMatrix>({get_qubit_index(it, i + offset - start_column, dim, start_index)});
                //solve.add_gate<NotMatrix>({i * n + it});

            vector<int> target_cnot_qubits;
            for (int it = i; it < j; it++)
                target_cnot_qubits.push_back(get_qubit_index(it, i + offset - start_column, dim, start_index));

            for (auto const &gate: add_row_op.gates) {
                Matrix gate_with_CNOTs = create_multiple_controlled_matrix(gate.matrix, j - i);
                vector<int> target_qubits = target_cnot_qubits;
                for (auto const &index: gate.target_indices)
                    target_qubits.push_back(index);
                
                pivor_search_and_row_reduce_op.add_gate({gate_with_CNOTs, target_qubits});
            }

            for (int it = i; it < j; it++)
                pivor_search_and_row_reduce_op.add_gate<NotMatrix>({get_qubit_index(it, i + offset - start_column, dim, start_index)});
        }

        //row reduce
        for (int j = i + 1; j < dim; j++) {
            Operator add_row_op;
            for (int it = i + 1; it < dim + 1; it++)
                add_row_op.add_gate<CNOTMatrix>({get_qubit_index(i, it + offset - start_column, dim, start_index), get_qubit_index(j, it + offset - start_column, dim, start_index)});

            for (auto const &gate: add_row_op.gates)             
                pivor_search_and_row_reduce_op.add_gate(create_controlled_gate(gate, get_qubit_index(j, i + offset - start_column, dim, start_index)));
        }

         for (const auto &gate: pivor_search_and_row_reduce_op.gates) 
            solve.add_gate(create_controlled_gate(gate, x_start + i + offset - start_column));
  
    }

    //back substitution
    for (int i = dim - 1; i >= 1 + start_column; i--) {
        Operator back_substitution_op;
        for (int j = i - 1; j >= 0; j--)
            back_substitution_op.add_gate<CCNOTMatrix>({get_qubit_index(j, i + offset - start_column, dim, start_index), get_qubit_index(i, dim + offset - start_column, dim, start_index), get_qubit_index(j, dim + offset - start_column, dim, start_index)});
        
        for (const auto &gate: back_substitution_op.gates) 
            solve.add_gate(create_controlled_gate(gate, x_start + i + offset - start_column));
    }

    return solve;
}


void ISD::run_example_system_of_equations() {
    cout << "Solving the following system of binary equations:\n";
    // cout << "x2 + x3 = 1\n";
    // cout << "x1 + x3 = 0\n";
    // cout << "x1 + x2 + x3 = 1\n\n";
    // cout << "Expected values: x1 = 0, x2 = 1, x3 = 0\n\n";

    // cout << "Input matrix:\n";
    // cout << "0 1 1 | 1\n";
    // cout << "1 0 1 | 0\n";
    // cout << "1 1 1 | 1\n\n";

    // Circuit ISD(12);
    // ISD.add_gate<NotMatrix>({1});
    // ISD.add_gate<NotMatrix>({2});
    // ISD.add_gate<NotMatrix>({3});
    // ISD.add_gate<NotMatrix>({5});
    // ISD.add_gate<NotMatrix>({6});
    // ISD.add_gate<NotMatrix>({7});
    // ISD.add_gate<NotMatrix>({8});
    // ISD.add_gate<NotMatrix>({9});
    // ISD.add_gate<NotMatrix>({11});
    // int dim = 3, all_columns = 3;
    // int offset = 0;

    // cout << "Input matrix:\n";
    // cout << "0 1 1 | 1\n";
    // cout << "1 0 1 | 0\n";
    // cout << "1 1 1 | 1\n\n";

    // Circuit ISD(12 + 3);
    // ISD.add_gate<NotMatrix>({1});
    // ISD.add_gate<NotMatrix>({2});
    // ISD.add_gate<NotMatrix>({3});
    // ISD.add_gate<NotMatrix>({5});
    // ISD.add_gate<NotMatrix>({6});
    // ISD.add_gate<NotMatrix>({7});
    // ISD.add_gate<NotMatrix>({8});
    // ISD.add_gate<NotMatrix>({9});
    // ISD.add_gate<NotMatrix>({11});

    // ISD.add_gate<NotMatrix>({12});
    // ISD.add_gate<NotMatrix>({13});
    // ISD.add_gate<NotMatrix>({14});
    
    // int dim = 3, all_columns = 3;
    // int offset = 0;

    // cout << "Input matrix:\n";
    // cout << "1 0 1 1 | 1\n";
    // cout << "0 1 0 1 | 0\n";
    // cout << "0 1 1 1 | 1\n\n";

    // Circuit ISD(15 + 4);
    // ISD.add_gate<NotMatrix>({0});
    // ISD.add_gate<NotMatrix>({4});
    // ISD.add_gate<NotMatrix>({5});
    // ISD.add_gate<NotMatrix>({6});
    // ISD.add_gate<NotMatrix>({8});
    // ISD.add_gate<NotMatrix>({9});
    // ISD.add_gate<NotMatrix>({10});
    // ISD.add_gate<NotMatrix>({11});
    // ISD.add_gate<NotMatrix>({12});
    // ISD.add_gate<NotMatrix>({14});

    // ISD.add_gate<NotMatrix>({16});
    // ISD.add_gate<NotMatrix>({17});
    // ISD.add_gate<NotMatrix>({18});
    // int dim = 3, all_columns = 4;
    // int offset = 3;


    cout << "Input matrix:\n";
    cout << "1 1 | 1\n";
    cout << "0 1 | 0\n";
    cout << "1 0 | 1\n\n";

    Circuit ISD(9 + 2);
    ISD.add_gate<NotMatrix>({0});
    ISD.add_gate<NotMatrix>({2});
    ISD.add_gate<NotMatrix>({3});
    ISD.add_gate<NotMatrix>({4});
    ISD.add_gate<NotMatrix>({6});
    ISD.add_gate<NotMatrix>({8});

    ISD.add_gate<NotMatrix>({9});
    ISD.add_gate<NotMatrix>({10});
    int dim = 3, all_columns = 2;
    int offset = -3;


    // cout << "Input matrix:\n";
    // cout << "0 0 0 | 0\n";
    // cout << "1 0 1 | 1\n";
    // cout << "0 1 1 | 1\n\n";

    // Circuit ISD(12);
    // ISD.add_gate<NotMatrix>({1});
    // ISD.add_gate<NotMatrix>({5});
    // ISD.add_gate<NotMatrix>({7});
    // ISD.add_gate<NotMatrix>({8});
    // ISD.add_gate<NotMatrix>({10});
    // ISD.add_gate<NotMatrix>({11});

  
    ISD.add_operator(sus(dim, 0, all_columns));
    ISD.run();
    vector<bool> read = ISD.magic_read();

    cout << "Output matrix:\n";
    cout << /*to_string(read[0 + offset]) +*/ " " + to_string(read[3 + offset]) + " " + to_string(read[6 + offset]) + " | " + to_string(read[9 + offset]) + "\n";
    cout << /*to_string(read[1 + offset]) +*/ " " + to_string(read[4 + offset]) + " " + to_string(read[7 + offset]) + " | " + to_string(read[10 + offset]) + "\n";
    cout << /*to_string(read[2 + offset]) + */" " + to_string(read[5 + offset]) + " " + to_string(read[8 + offset]) + " | " + to_string(read[11 + offset]) + "\n\n";

    cout << "Calculated values: x1 = " + to_string(read[9 + offset]) + ", x2 = " + to_string(read[10 + offset]) + ", x3 = " + to_string(read[11 + offset]) + "\n\n";
}

int ISD::stupid(vector<int> B, vector<int> M, int n, int k) {
    //cout << "Performing superposition equation solver for H|s:\n";
    // cout << "1 0 0 1 0 0 | 1\n";
    // cout << "0 1 0 1 0 0 | 1\n";
    // cout << "0 0 1 0 1 1 | 0\n\n";


    // cout << "1 0 0 0 1 | 0\n";
    // cout << "0 1 0 1 0 | 0\n";
    // cout << "0 0 1 1 1 | 1\n\n";

    int k_bits = bit_width((unsigned int)(n - k));

    Circuit superposition_solver(n + 1 + (k + 1) * (n - k) + k + k_bits + k_bits - 1);

    superposition_solver.add_gate<NotMatrix>({B[0]});
    superposition_solver.add_gate<NotMatrix>({B[1]});
    superposition_solver.add_gate<NotMatrix>({B[2]});


    // superposition_solver.add_operator(get_optimized_example_init());
    int matrix_start = n + 1;
    for (auto &m: M)
        superposition_solver.add_gate<NotMatrix>({matrix_start + m});
    // superposition_solver.add_gate<NotMatrix>({matrix_start + 1});
    // superposition_solver.add_gate<NotMatrix>({matrix_start + 2});
    // superposition_solver.add_gate<NotMatrix>({matrix_start + 3});
    // superposition_solver.add_gate<NotMatrix>({matrix_start + 5});
    // superposition_solver.add_gate<NotMatrix>({matrix_start + 8});


    superposition_solver.add_operator(get_optimisation_operator(n, k, n + 1 + (k + 1) * (n - k) + k));
    // superposition_solver.add_operator(controlled_solve_from_back(k, n - k, n + 1));
    superposition_solver.add_operator(sus(n - k, n + 1, k));

    superposition_solver.run();
    vector<bool> read = superposition_solver.magic_read();
    int res = 0;
    for (int i = 0; i < n - k; i++)
        res += (read[n + 1 + k * (n - k) + i]);

    //cout << res << endl << endl;
    return res;

    // for (int i = 0; i <n + 1 + (k + 1) * (n - k) + k + k_bits + k_bits - 1; i++ )
    // cout << read[i] << " ";
    // cout << endl;

}


void ISD::brut() {
    vector<vector<int>> K {{1, 1, 0}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

    vector<vector<vector<int>>> RES(21);
    int counter = 0;
    int n = 3, k = 3;

    for (int first_column = 0; first_column < 4; first_column++) {
        for (int second_column = 0; second_column < 4; second_column++) {
            if (second_column == first_column) continue;
            for (int third_column = 0; third_column < 4; third_column++) {
                if (third_column == first_column || third_column == second_column) continue;

                for (int syndrom =0; syndrom < 6; syndrom++) {

                    vector<int> M;
                    if (K[first_column][0] == 1) M.push_back(0);
                    if (K[first_column][1] == 1) M.push_back(1);
                    if (K[first_column][2] == 1) M.push_back(2);

                    if (K[second_column][0] == 1) M.push_back(3);
                    if (K[second_column][1] == 1) M.push_back(4);
                    if (K[second_column][2] == 1) M.push_back(5);

                    if (K[third_column][0] == 1) M.push_back(6);
                    if (K[third_column][1] == 1) M.push_back(7);
                    if (K[third_column][2] == 1) M.push_back(8);

                    switch (syndrom) {
                    case 0:
                        M.push_back(9);
                        break;
                    case 1:
                        M.push_back(10);
                        break;
                    case 2:
                        M.push_back(11);
                        break;
                    case 3:
                        if (K[first_column][0] == 1) M.push_back(9);
                        if (K[first_column][1] == 1) M.push_back(10);
                        if (K[first_column][2] == 1) M.push_back(11);
                        break;
                    case 4:
                        if (K[second_column][0] == 1) M.push_back(9);
                        if (K[second_column][1] == 1) M.push_back(10);
                        if (K[second_column][2] == 1) M.push_back(11);
                        break;
                    case 5:
                        if (K[third_column][0] == 1) M.push_back(9);
                        if (K[third_column][1] == 1) M.push_back(10);
                        if (K[third_column][2] == 1) M.push_back(11);
                        break;
                    default:
                        break;
                    }


                    if (first_column == 0 && second_column == 2 && third_column == 1 && syndrom == 0) {
                        cout << "";
                    } else continue;

                    int one_cnt = 0;
                    for (int mask = 0; mask < (1 << 6); mask++) {
                        if (popcount((unsigned int)(mask)) != 3) continue;

                        vector<int> B;
                        for (int i = 0; i < 6; i++) {
                            if ((mask >> i) & 1)
                                B.push_back(i);
                        }

                        cout << ++counter << "/" << 2880 << endl;
                        if (stupid(B, M, n, k) == 1)
                            one_cnt++;
                    }

                    RES[one_cnt].push_back({first_column, second_column, third_column, syndrom});               
                }
            }
        }
    }

    for (int i = 20; i >= 0; i--) {
        cout << "================ " << i << " ================" << endl;
        for (auto &v: RES[i]) {
            cout << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << endl;
        }
        cout << endl;
    }
}

void ISD::apply_classic_optimisation_operator(vector<bool>B, vector<bool> &M, vector<bool> &X, int n, int k) {
    int c = n - k;

    for (int i = 0; i < n - k; i++) {
        for (int j = 0; j < i; j++)
            if (B[i] == 1 && c == n - k - j)
                for (int it = 0; it < k + 1; it++)
                    swap(M[get_qubit_index(i, it, n - k, 0)], M[get_qubit_index(j, it, n - k, 0)]);

        if (B[i] == 1) c--;
    }

    for (int i = k - 1; i >= 0; i--) {
        if (B[n - k + i] == 1) {
            c--;
            X[i] = 1;
            for (int j = i; j < k - 1; j++) {
                for (int it = 0; it < n - k; it++)
                    swap(M[get_qubit_index(it, j, n - k, 0)], M[get_qubit_index(it, j + 1, n - k, 0)]);

                swap(X[j], X[j + 1]);
            }
        }
    }
}

void ISD::apply_classic_sus(int dim, vector<bool> &M, vector<bool> &X, int all_columns) {
    int offset = max(all_columns - dim, 0);
    int start_column = max(0, dim - all_columns);

    for (int i = 0 + start_column; i < dim - 1; i++) {
        if (X[i + offset - start_column]) {
            //pivot search
            for (int j = i + 1; j < dim; j++) {
                bool all_zeros = 1;
                for (int it = i; it < j; it++) 
                    if (M[get_qubit_index(it, i + offset - start_column, dim, 0)] == 1)
                        all_zeros = 0;

                if (all_zeros) {
                    for (int it = i + 1; it < dim + 1; it++) {
                        //M[get_qubit_index(i, it + offset - start_column, dim, 0)] ^= M[get_qubit_index(j, it + offset - start_column, dim , 0)];
                        auto idx_i_695 = get_qubit_index(i, it + offset - start_column, dim, 0);
                        auto idx_j_695 = get_qubit_index(j, it + offset - start_column, dim, 0);
                        M[idx_i_695] = M[idx_i_695] ^ M[idx_j_695];
                    }
                }
            }

            //row reduce
            for (int j = i + 1; j < dim; j++) {
                if (M[get_qubit_index(j, i + offset - start_column, dim, 0)] == 1) {
                    for (int it = i + 1; it < dim + 1; it++) {
                        //M[get_qubit_index(j, it + offset - start_column, dim, 0)] ^= M[get_qubit_index(i, it + offset - start_column, dim , 0)];
                        auto idx_j_702 = get_qubit_index(j, it + offset - start_column, dim, 0);
                        auto idx_i_702 = get_qubit_index(i, it + offset - start_column, dim, 0);
                        M[idx_j_702] = M[idx_j_702] ^ M[idx_i_702];
                    }
                }
            }

        }  
    }

    //back substitution
    for (int i = dim - 1; i >= max(1, start_column); i--) {
        if (X[i + offset - start_column]) {
            for (int j = i - 1; j >= 0; j--) {
                if (M[get_qubit_index(j, i + offset - start_column, dim, 0)] == 1) {
                    //M[get_qubit_index(j, dim + offset - start_column, dim, 0)] ^= M[get_qubit_index(i, dim + offset - start_column, dim, 0)];
                    auto idx_j_712 = get_qubit_index(j, dim + offset - start_column, dim, 0);
                    auto idx_i_712 = get_qubit_index(i, dim + offset - start_column, dim, 0);
                    M[idx_j_712] = M[idx_j_712] ^ M[idx_i_712];
                }
            }
        }
    }

}

int ISD::classic_stupid(vector<bool>B, vector<bool> M, int n, int k) {
    vector<bool> X(k);
    apply_classic_optimisation_operator(B, M, X, n, k);
    apply_classic_sus(n - k, M, X, k);


    int res = 0;
    for (int i = 0; i < n - k; i++)
        res += (int)(M[k * (n - k) + i]);

    return res;
}

void ISD::classic_smol_brut() {
    vector<vector<int>> K {{1, 1, 0}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

    vector<vector<vector<int>>> RES(11);
    int counter = 0;
    int n = 5, k = 2;

    for (int first_column = 0; first_column < 4; first_column++) {
        for (int second_column = 0; second_column < 4; second_column++) {
            if (second_column == first_column) continue;

            for (int syndrom =0; syndrom < 5; syndrom++) {

                vector<bool> M;
                M.push_back(K[first_column][0]);
                M.push_back(K[first_column][1]);
                M.push_back(K[first_column][2]);

                M.push_back(K[second_column][0]);
                M.push_back(K[second_column][1]);
                M.push_back(K[second_column][2]);

                switch (syndrom) {
                case 0:
                    M.push_back(1);
                    M.push_back(0);
                    M.push_back(0);
                    break;
                case 1:
                    M.push_back(0);
                    M.push_back(1);
                    M.push_back(0);
                    break;
                case 2:
                    M.push_back(0);
                    M.push_back(0);
                    M.push_back(1);
                    break;
                case 3:
                    M.push_back(K[first_column][0]);
                    M.push_back(K[first_column][1]);
                    M.push_back(K[first_column][2]);
                    break;
                case 4:
                    M.push_back(K[second_column][0]);
                    M.push_back(K[second_column][1]);
                    M.push_back(K[second_column][2]);
                    break;
                default:
                    break;
                }

                int one_cnt = 0;
                for (int mask = 0; mask < (1 << 5); mask++) {
                    if (popcount((unsigned int)(mask)) != 3) continue;

                    vector<bool> B;
                    for (int i = 0; i < 5; i++)
                        B.push_back((mask >> i) & 1);


                    cout << ++counter << "/" << 600 << endl;

                    if (classic_stupid(B, M, n, k) == 1)
                        one_cnt++;
                }

                RES[one_cnt].push_back({first_column, second_column, syndrom});               
            }
        }
    }

    for (int i = 10; i >=0; i--) {
        cout << "================ " << i << " ================" << endl;
        for (auto &v: RES[i]) {
            cout << v[0] << " " << v[1] << " " << v[2] << endl;
        }
        cout << endl;
    }
}

void ISD::classic_brut() {
    vector<vector<int>> K {{1, 1, 0}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

    vector<vector<vector<int>>> RES(21);
    int counter = 0;
    int n = 6, k = 3;

    for (int first_column = 0; first_column < 4; first_column++) {
        for (int second_column = 0; second_column < 4; second_column++) {
            if (second_column == first_column) continue;
            for (int third_column = 0; third_column < 4; third_column++) {
                if (third_column == first_column || third_column == second_column) continue;

                for (int syndrom =0; syndrom < 6; syndrom++) {

                    vector<bool> M;
                    M.push_back(K[first_column][0]);
                    M.push_back(K[first_column][1]);
                    M.push_back(K[first_column][2]);

                    M.push_back(K[second_column][0]);
                    M.push_back(K[second_column][1]);
                    M.push_back(K[second_column][2]);

                    M.push_back(K[third_column][0]);
                    M.push_back(K[third_column][1]);
                    M.push_back(K[third_column][2]);

                    switch (syndrom) {
                    case 0:
                        M.push_back(1);
                        M.push_back(0);
                        M.push_back(0);
                        break;
                    case 1:
                        M.push_back(0);
                        M.push_back(1);
                        M.push_back(0);
                        break;
                    case 2:
                        M.push_back(0);
                        M.push_back(0);
                        M.push_back(1);
                        break;
                    case 3:
                        M.push_back(K[first_column][0]);
                        M.push_back(K[first_column][1]);
                        M.push_back(K[first_column][2]);
                        break;
                    case 4:
                        M.push_back(K[second_column][0]);
                        M.push_back(K[second_column][1]);
                        M.push_back(K[second_column][2]);
                        break;
                    case 5:
                        M.push_back(K[third_column][0]);
                        M.push_back(K[third_column][1]);
                        M.push_back(K[third_column][2]);
                        break;
                    default:
                        break;
                    }

                    int one_cnt = 0;
                    for (int mask = 0; mask < (1 << 6); mask++) {
                        if (popcount((unsigned int)(mask)) != 3) continue;

                        vector<bool> B;
                        for (int i = 0; i < 6; i++)
                            B.push_back((mask >> i) & 1);


                        cout << ++counter << "/" << 2880 << endl;
                        if (classic_stupid(B, M, n, k) == 1)
                            one_cnt++;
                    }

                    RES[one_cnt].push_back({first_column, second_column, third_column, syndrom});               
                }
            }
        }
    }

    for (int i = 20; i >= 0; i--) {
        cout << "================ " << i << " ================" << endl;
        for (auto &v: RES[i]) {
            cout << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << endl;
        }
        cout << endl;
    }
}


void ISD::classic_big_brut() {
    vector<vector<int>> K {{1, 1, 0}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

    vector<vector<vector<int>>> RES(36);
    int counter = 0;
    int n = 7, k = 4;

    for (int first_column = 0; first_column < 4; first_column++) {
        for (int second_column = 0; second_column < 4; second_column++) {
            if (second_column == first_column) continue;
            for (int third_column = 0; third_column < 4; third_column++) {
                if (third_column == first_column || third_column == second_column) continue;
                for (int forth_column = 0; forth_column < 4; forth_column++) {
                    if (forth_column == first_column || forth_column == second_column || forth_column == third_column) continue;
                    for (int syndrom = 0; syndrom < 7; syndrom++) {

                        vector<bool> M;
                        M.push_back(K[first_column][0]);
                        M.push_back(K[first_column][1]);
                        M.push_back(K[first_column][2]);

                        M.push_back(K[second_column][0]);
                        M.push_back(K[second_column][1]);
                        M.push_back(K[second_column][2]);

                        M.push_back(K[third_column][0]);
                        M.push_back(K[third_column][1]);
                        M.push_back(K[third_column][2]);

                        M.push_back(K[forth_column][0]);
                        M.push_back(K[forth_column][1]);
                        M.push_back(K[forth_column][2]);

                        switch (syndrom) {
                        case 0:
                            M.push_back(1);
                            M.push_back(0);
                            M.push_back(0);
                            break;
                        case 1:
                            M.push_back(0);
                            M.push_back(1);
                            M.push_back(0);
                            break;
                        case 2:
                            M.push_back(0);
                            M.push_back(0);
                            M.push_back(1);
                            break;
                        case 3:
                            M.push_back(K[first_column][0]);
                            M.push_back(K[first_column][1]);
                            M.push_back(K[first_column][2]);
                            break;
                        case 4:
                            M.push_back(K[second_column][0]);
                            M.push_back(K[second_column][1]);
                            M.push_back(K[second_column][2]);
                            break;
                        case 5:
                            M.push_back(K[third_column][0]);
                            M.push_back(K[third_column][1]);
                            M.push_back(K[third_column][2]);
                            break;
                        case 6:
                            M.push_back(K[forth_column][0]);
                            M.push_back(K[forth_column][1]);
                            M.push_back(K[forth_column][2]);
                            break;
                        default:
                            break;
                        }
                        

                        if (first_column == 0 && second_column == 3 && third_column == 1 && forth_column == 2 && syndrom == 3) {
                            cout << "";
                        }

                        int one_cnt = 0;
                        for (int mask = 0; mask < (1 << 7); mask++) {
                            if (popcount((unsigned int)(mask)) != 3) continue;

                            vector<bool> B;
                            for (int i = 0; i < 7; i++)
                                B.push_back((mask >> i) & 1);

                            cout << ++counter << "/" << 5880 << endl;
                            if (classic_stupid(B, M, n, k) == 1)
                                one_cnt++;
                        }

                        RES[one_cnt].push_back({first_column, second_column, third_column, forth_column, syndrom});               
                    }
                }
            }
        }
    }

    for (int i = 35; i >= 0; i--) {
        cout << "================ " << i << " ================" << endl;
        for (auto &v: RES[i]) {
            cout << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << " " << v[4] << endl;
        }
        cout << endl;
    }
}

void ISD::run_defaults() {
    cout << "======================= ISD ======================\n";

    classic_smol_brut();
    //classic_brut();
    //classic_big_brut();
    // ISDInfo small_isd_info = get_small_isd_info();
    // run_superposition_solver(small_isd_info);
    // grover.run_grover(get_grover_info(small_isd_info));

    // ISDInfo isd_info = get_isd_info();
    // run_superposition_solver(isd_info);
    // grover.run_grover(get_grover_info(isd_info));

    // ISDInfo big_isd_info = get_big_isd_info();
    // run_superposition_solver(big_isd_info);
    // grover.run_grover(get_grover_info(big_isd_info));
}