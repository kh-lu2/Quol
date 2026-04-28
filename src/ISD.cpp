#include "../inc/ISD.h"
#include "../inc/quantum_utils.h"
#include <iostream>
#include <map>

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

Operator ISD::solve(int n, int k, int start_index) {
    Operator solve;

    for (int i = 0; i < n - 1; i++) {
        //pivot search
        for (int j = i + 1; j < k; j++) {            
            Operator add_row_op;
            for (int it = i + 1; it < n + 1; it++)
                add_row_op.add_gate<CNOTMatrix>({get_qubit_index(j, it, k, start_index), get_qubit_index(i, it, k, start_index)});

            for (int it = i; it < j; it++)
                solve.add_gate<NotMatrix>({get_qubit_index(it, i, k, start_index)});
                //solve.add_gate<NotMatrix>({i * n + it});

            vector<int> target_cnot_qubits;
            for (int it = i; it < j; it++)
                target_cnot_qubits.push_back(get_qubit_index(it, i, k, start_index));

            for (auto const &gate: add_row_op.gates) {
                Matrix gate_with_CNOTs = createMultipleControlledMatrix(gate.matrix, j - i);
                vector<int> target_qubits = target_cnot_qubits;
                for (auto const &index: gate.target_indices)
                    target_qubits.push_back(index);
                
                solve.add_gate({gate_with_CNOTs, target_qubits});
            }

            for (int it = i; it < j; it++)
                solve.add_gate<NotMatrix>({get_qubit_index(it, i, k, start_index)});
        }

        //row reduce
        for (int j = i + 1; j < n; j++) {
            Operator add_row_op;
            for (int it = i + 1; it < n + 1; it++)
                add_row_op.add_gate<CNOTMatrix>({get_qubit_index(i, it, k, start_index), get_qubit_index(j, it, k, start_index)});

            for (auto const &gate: add_row_op.gates)             
                solve.add_gate(createControlledGate(gate, get_qubit_index(j, i, k, start_index)));
        }
    }

    //back substitution
    for (int i = n - 1; i >= 1; i--)
        for (int j = i - 1; j >= 0; j--)
            solve.add_gate<CCNOTMatrix>({get_qubit_index(j, i, k, start_index), get_qubit_index(i, n, k, start_index), get_qubit_index(j, n, k, start_index)});

    return solve;
}

Operator ISD::solve_from_back(int n, int k, int start_index) {
    Operator solve;

    for (int i = n - 1; i >= 1; i--) {
        //pivot search
        for (int j = i - 1; j >= 0; j--) {            
            Operator add_row_op;
            for (int it = i - 1; it >= 0; it--)
                add_row_op.add_gate<CNOTMatrix>({get_qubit_index(j, it, k, start_index), get_qubit_index(i, it, k, start_index)});
            add_row_op.add_gate<CNOTMatrix>({get_qubit_index(j, n, k, start_index), get_qubit_index(i, n, k, start_index)});

            for (int it = i; it > j; it--)
                solve.add_gate<NotMatrix>({get_qubit_index(it, i, k, start_index)});

            vector<int> target_cnot_qubits;
            for (int it = i; it > j; it--)
                target_cnot_qubits.push_back(get_qubit_index(it, i, k, start_index));

            for (auto const &gate: add_row_op.gates) {
                Matrix gate_with_CNOTs = createMultipleControlledMatrix(gate.matrix, i - j);
                vector<int> target_qubits = target_cnot_qubits;
                for (auto const &index: gate.target_indices)
                    target_qubits.push_back(index);
                
                solve.add_gate({gate_with_CNOTs, target_qubits});
            }

            for (int it = i; it > j; it--)
                solve.add_gate<NotMatrix>({get_qubit_index(it, i, k, start_index)});
        }

        //row reduce
        for (int j = i - 1; j >= 0; j--) {
            Operator add_row_op;
            for (int it = i - 1; it >= 0; it--)
                add_row_op.add_gate<CNOTMatrix>({get_qubit_index(i, it, k, start_index), get_qubit_index(j, it, k, start_index)});
            add_row_op.add_gate<CNOTMatrix>({get_qubit_index(i, n, k, start_index), get_qubit_index(j, n, k, start_index)});


            for (auto const &gate: add_row_op.gates)             
                solve.add_gate(createControlledGate(gate, get_qubit_index(j, i, k, start_index)));
        }
    }

    //back substitution 
    for (int i = 0; i < n - 1; i++)
        for (int j = i + 1; j < n; j++)
            solve.add_gate<CCNOTMatrix>({get_qubit_index(j, i, k, start_index), get_qubit_index(i, n, k, start_index), get_qubit_index(j, n, k, start_index)});

    return solve;
}

Operator ISD::controlled_solve_from_back(int n, int k, int start_index) {
    Operator solve;
    int x_start = start_index + (k + 1) * (n - k);

    for (int i = n - 1; i >= 1; i--) {
        Operator pivot_search_and_row_reduce;
        //pivot search
        for (int j = i - 1; j >= 0; j--) {
            Operator add_row_op;
            for (int it = i - 1; it >= 0; it--)
                add_row_op.add_gate<CNOTMatrix>({get_qubit_index(j, it, k, start_index), get_qubit_index(i, it, k, start_index)});
            add_row_op.add_gate<CNOTMatrix>({get_qubit_index(j, n, k, start_index), get_qubit_index(i, n, k, start_index)});

            for (int it = i; it > j; it--)
                pivot_search_and_row_reduce.add_gate<NotMatrix>({get_qubit_index(it, i, k, start_index)});

            vector<int> target_cnot_qubits;
            for (int it = i; it > j; it--)
                target_cnot_qubits.push_back(get_qubit_index(it, i, k, start_index));

            for (auto const &gate: add_row_op.gates) {
                Matrix gate_with_CNOTs = createMultipleControlledMatrix(gate.matrix, i - j);
                vector<int> target_qubits = target_cnot_qubits;
                for (auto const &index: gate.target_indices)
                    target_qubits.push_back(index);
                
                pivot_search_and_row_reduce.add_gate({gate_with_CNOTs, target_qubits});
            }

            for (int it = i; it > j; it--)
                pivot_search_and_row_reduce.add_gate<NotMatrix>({get_qubit_index(it, i, k, start_index)});
        }

        //row reduce
        for (int j = i - 1; j >= 0; j--) {
            Operator add_row_op;
            for (int it = i - 1; it >= 0; it--)
                add_row_op.add_gate<CNOTMatrix>({get_qubit_index(i, it, k, start_index), get_qubit_index(j, it, k, start_index)});
            add_row_op.add_gate<CNOTMatrix>({get_qubit_index(i, n, k, start_index), get_qubit_index(j, n, k, start_index)});


            for (auto const &gate: add_row_op.gates)             
                pivot_search_and_row_reduce.add_gate(createControlledGate(gate, get_qubit_index(j, i, k, start_index)));
        }

        for (const auto &gate: pivot_search_and_row_reduce.gates) 
            solve.add_gate(createControlledGate(gate, x_start + i));
    }

    //back substitution 
    for (int i = 0; i < n - 1; i++) {
        Operator back_substitution;
        for (int j = i + 1; j < n; j++)
            back_substitution.add_gate<CCNOTMatrix>({get_qubit_index(j, i, k, start_index), get_qubit_index(i, n, k, start_index), get_qubit_index(j, n, k, start_index)});

        for (const auto &gate: back_substitution.gates) 
            solve.add_gate(createControlledGate(gate, x_start + i));
    }

    return solve;
}


void ISD::run_example_system_of_equations() {
    cout << "Solving the following system of binary equations:\n";
    cout << "x2 + x3 = 1\n";
    cout << "x1 + x3 = 0\n";
    cout << "x1 + x2 + x3 = 1\n\n";
    cout << "Expected values: x1 = 0, x2 = 1, x3 = 0\n\n";

    cout << "Input matrix:\n";
    cout << "0 1 1 | 1\n";
    cout << "1 0 1 | 0\n";
    cout << "1 1 1 | 1\n\n";

    Circuit ISD(12);
    ISD.add_gate<NotMatrix>({1});
    ISD.add_gate<NotMatrix>({2});
    ISD.add_gate<NotMatrix>({3});
    ISD.add_gate<NotMatrix>({5});
    ISD.add_gate<NotMatrix>({6});
    ISD.add_gate<NotMatrix>({7});
    ISD.add_gate<NotMatrix>({8});
    ISD.add_gate<NotMatrix>({9});
    ISD.add_gate<NotMatrix>({11});

    int n = 3, k = 3;
    ISD.add_operator(solve(n, k, 0));
    ISD.run();
    vector<bool> read = ISD.magic_read();

    cout << "Output matrix:\n";
    cout << to_string(read[0]) + " " + to_string(read[3]) + " " + to_string(read[6]) + " | " + to_string(read[9]) + "\n";
    cout << to_string(read[1]) + " " + to_string(read[4]) + " " + to_string(read[7]) + " | " + to_string(read[10]) + "\n";
    cout << to_string(read[2]) + " " + to_string(read[5]) + " " + to_string(read[8]) + " | " + to_string(read[11]) + "\n\n";

    cout << "Calculated values: x1 = " + to_string(read[9]) + ", x2 = " + to_string(read[10]) + ", x3 = " + to_string(read[11]) + "\n\n";
}

void ISD::run_example_system_of_equations_from_back() {
    cout << "Solving the following system of binary equations:\n";
    cout << "x2 + x3 = 1\n";
    cout << "x1 + x3 = 0\n";
    cout << "x1 + x2 = 1\n\n";
    cout << "Expected values: x1 = 0, x2 = 1, x3 = 0\n\n";

    cout << "Input matrix:\n";
    cout << "0 1 1 | 1\n";
    cout << "1 0 1 | 0\n";
    cout << "1 1 0 | 1\n\n";

    Circuit ISD(12);
    ISD.add_gate<NotMatrix>({1});
    ISD.add_gate<NotMatrix>({2});
    ISD.add_gate<NotMatrix>({3});
    ISD.add_gate<NotMatrix>({5});
    ISD.add_gate<NotMatrix>({6});
    ISD.add_gate<NotMatrix>({7});
    ISD.add_gate<NotMatrix>({9});
    ISD.add_gate<NotMatrix>({11});

    int n = 3, k = 3;
    ISD.add_operator(solve_from_back(n, k, 0));
    ISD.run();
    vector<bool> read = ISD.magic_read();

    cout << "Output matrix:\n";
    cout << to_string(read[0]) + " " + to_string(read[3]) + " " + to_string(read[6]) + " | " + to_string(read[9]) + "\n";
    cout << to_string(read[1]) + " " + to_string(read[4]) + " " + to_string(read[7]) + " | " + to_string(read[10]) + "\n";
    cout << to_string(read[2]) + " " + to_string(read[5]) + " " + to_string(read[8]) + " | " + to_string(read[11]) + "\n\n";

    cout << "Calculated values: x1 = " + to_string(read[9]) + ", x2 = " + to_string(read[10]) + ", x3 = " + to_string(read[11]) + "\n\n";
}


void ISD::run_toy_ISD() {
    cout << "Performing ISD for...\n";
    cout << "Matrix H:\n";
    cout << "1 1 0 1 0 0\n";
    cout << "1 0 1 0 1 0\n";
    cout << "0 1 1 0 0 1\n\n";
    cout << "Syndrome: 1 0 0\n\n";
    cout << "Exprected error vector: 0 0 0 1 0 0\n\n";
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
            create_uniform_superposition_op.add_gate(createControlledGate({RotationMatrix(theta), {i}}, ancilla_after_c_index + c_bits - 2));

            create_uniform_superposition_op.add_operator(reverse_check_c_ones);
            create_uniform_superposition_op.add_operator(reverse_operator(nots_to_j));
        }

        for (const auto &gate: c_minus_one.gates)
            create_uniform_superposition_op.add_gate(createControlledGate(gate, i));
    }

    return create_uniform_superposition_op;
}

int ISD::get_qubit_index(int row, int column, int k, int start_index){
    return start_index + column * k + row;
}

Operator ISD::swap_chosen_columns_to_front_operator(int n, int k) {
    Operator column_swap_op;
    for (int i = 0; i < n; i++) {
        for (int j = i - 1; j >= 0; j--) {
            Operator swap_col;
            for (int it = 0; it < k; it++)
                //swap_col.add_gate<SwapMatrix>({get_qubit_index(it, j, k, 2 * n), get_qubit_index(it, j + 1, k, 2 * n)});
                swap_col.add_gate<SwapMatrix>({get_qubit_index(it, j, k, n), get_qubit_index(it, j + 1, k, n)});

            for (auto &gate: swap_col.gates)
                column_swap_op.add_gate(createControlledGate(gate, i));
        }
    }

    return column_swap_op;
}

Operator ISD::get_example_init() {
    // cout << "1 1 0 1 0 0 | 1\n";
    // cout << "1 0 1 0 1 0 | 0\n";
    // cout << "0 1 1 0 0 1 | 0\n\n";
    int n = 3;
    Operator init;


    init.add_gate<NotMatrix>({n + 0});
    init.add_gate<NotMatrix>({n + 2});
    init.add_gate<NotMatrix>({n + 4});
    init.add_gate<NotMatrix>({n + 5});
    init.add_gate<NotMatrix>({n + 6});
    init.add_gate<NotMatrix>({n + 10});
    init.add_gate<NotMatrix>({n + 14});
    init.add_gate<NotMatrix>({n + 16});
    init.add_gate<NotMatrix>({n + 17});

    return init;
}

Operator ISD::get_optimized_example_init() {
    // cout << "1 0 0 0 1 1 | 1\n";
    // cout << "0 1 0 1 0 1 | 0\n";
    // cout << "0 0 1 1 1 0 | 1\n\n";
    int n = 6;
    Operator init;
    init.add_gate<NotMatrix>({n + 1});
    init.add_gate<NotMatrix>({n + 2});
    init.add_gate<NotMatrix>({n + 3});
    init.add_gate<NotMatrix>({n + 5});
    init.add_gate<NotMatrix>({n + 6});
    init.add_gate<NotMatrix>({n + 7});
    init.add_gate<NotMatrix>({n + 9});
    init.add_gate<NotMatrix>({n + 11});

    return init;
}

int ISD::uniform_solver() {
    int n = 5, k = 3;
    int k_bits = bit_width((unsigned int)(k));

    //Circuit superposition_solver(n + n + k_bits + k_bits - 1 + (n + 1) * k);
    Circuit superposition_solver(n + k_bits + k_bits - 1 + (n + 1) * k);

    //superposition_solver.add_operator(get_uniform_superposition_operator(n, k, 2 * n + (n + 1) * k));
    superposition_solver.add_operator(get_uniform_superposition_operator(n, k, n + (n + 1) * k));
    superposition_solver.add_operator(get_example_init());
    superposition_solver.add_operator(swap_chosen_columns_to_front_operator(n, k));
    //superposition_solver.add_operator(solve(n, k, 2 * n));
    superposition_solver.add_operator(solve(n, k, n));

    superposition_solver.run();
    vector<bool> read = superposition_solver.magic_read();

    int res = 0;
    for (int i = 0; i < k; i++)
        res += (read[get_qubit_index(0, n, k, n) + i]);

    return res;
}

Operator ISD::get_optimisation_operator(int n, int k, int ancilla_start_index) {
    Operator optimisation;

    int c_bits = bit_width((unsigned int)(n - k));
    optimisation.add_operator(QuantumUtils::get_NOTs(c_bits, n - k, ancilla_start_index, 1));

    int ancilla_after_c_index = ancilla_start_index + c_bits;

    Operator check_c_ones = QuantumUtils::check_if_all_ones(c_bits, ancilla_start_index, ancilla_after_c_index);
    Operator reverse_check_c_ones = reverse_operator(check_c_ones);
    Operator c_minus_one = QuantumUtils::get_minus_one_operator(c_bits, ancilla_start_index);

    int x_start = ancilla_start_index - (n - k);

    for (int i = 0; i < n - k; i++) {
        for (int j = 0; j < i - 1; j++) {
            Operator nots_to_n_minus_k_minus_j_plus_one = QuantumUtils::get_NOTs(c_bits, n - k - j + 1, ancilla_start_index, 0);
            Operator reverse_nots = reverse_operator(nots_to_n_minus_k_minus_j_plus_one);

            optimisation.add_operator(nots_to_n_minus_k_minus_j_plus_one);
            optimisation.add_operator(check_c_ones);

            Operator swap_rows;
            for (int it = 0; it < k + 1; it++)
                swap_rows.add_gate<SwapMatrix>({get_qubit_index(i, it, k, n), get_qubit_index(j, it, k, n)});

            for (const auto &gate: swap_rows.gates) 
                optimisation.add_gate(createControlledGate(createControlledGate(gate, ancilla_after_c_index + c_bits - 2), i));
        }

        for (const auto &gate: c_minus_one.gates)
            optimisation.add_gate(createControlledGate(gate, i));
    }

    for (int i = n - 1; i >= n - k; i--) {
        Operator op;
        op.add_operator(c_minus_one);
        op.add_gate<NotMatrix>({x_start + i});

        for (int j = i; j < n - 1; j++) {
            for (int it = 0; it < n - k; it++)
                op.add_gate<SwapMatrix>({get_qubit_index(it, j, k, n), get_qubit_index(it, j + 1, k, n)});

            op.add_gate<SwapMatrix>({x_start + j, x_start + j + 1});
        }

        for (const auto &gate: op.gates)
            optimisation.add_gate(createControlledGate(gate, i));
    }

    return optimisation;
}

int ISD::optimized_uniform_solver() {
    int n = 6, k = 3;
    int k_bits = bit_width((unsigned int)(k));

    Circuit superposition_solver(n + (k + 1) * (n - k) + n - k + k_bits + k_bits - 1);

    superposition_solver.add_operator(get_uniform_superposition_operator(n, k, n + (k + 1) * (n - k) + n - k));
    superposition_solver.add_operator(get_optimized_example_init());
    superposition_solver.add_operator(get_optimisation_operator(n, k, n + (k + 1) * (n - k) + n - k));
    superposition_solver.add_operator(solve_from_back(n, k, n));

    superposition_solver.run();
    vector<bool> read = superposition_solver.magic_read();

    int res = 0;
    for (int i = 0; i < k; i++)
        res += (read[get_qubit_index(0, n, k, n) + i]);

    return res;
}

void ISD::run_uniform_solver() {
    cout << "Performing superposition equation solver for H|s:\n";
    cout << "1 1 0 1 0 0 | 1\n";
    cout << "1 0 1 0 1 0 | 0\n";
    cout << "0 1 1 0 0 1 | 0\n\n";

    int it = 1;

    map<int, int> M;
    for (int i = 0; i < it; i++)
        M[uniform_solver()]++;

    for (auto &[res, cnt]: M) 
        cout << res << ": " << double(cnt) / it * 100 << "%\n";
}

void ISD::run_optimized_uniform_solver() {
    cout << "Performing superposition equation solver for H|s:\n";
    cout << "1 0 0 0 1 1 | 1\n";
    cout << "0 1 0 1 0 1 | 0\n";
    cout << "0 0 1 1 1 0 | 1\n\n";

    int it = 1;

    map<int, int> M;
    for (int i = 0; i < it; i++)
        M[optimized_uniform_solver()]++;

    for (auto &[res, cnt]: M) 
        cout << res << ": " << double(cnt) / it * 100 << "%\n";
}


void ISD::run_defaults() {
    cout << "======================= ISD ======================\n";
    //run_example_system_of_equations();
    //run_example_system_of_equations_from_back();
    run_optimized_uniform_solver();
}