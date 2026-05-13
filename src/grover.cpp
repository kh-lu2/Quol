#define _USE_MATH_DEFINES
#include <map>
#include <iostream>
#include "../inc/grover.h"
#include "../inc/quantum_utils.h"

Operator Grover::get_Fn(int n, const Operator& fx) {
    Operator Fn;

    Fn.add_operator(fx);
    Operator if_one_op = check_if_all_ones(n);
    Fn.add_operator(if_one_op);

    Fn.add_gate<CNOTMatrix>({2 * n - 1, n});

    Fn.add_operator(reverse_operator(if_one_op));
    Fn.add_operator(reverse_operator(fx));

    return Fn;
}

Operator Grover::get_Fn_for_ISD(int n, const Operator& fx, int syndrome_start, int ancilla_start) {
    Operator Fn;

    Fn.add_operator(fx);

    Operator is_true = check_if_weight_one(syndrome_start, ancilla_start);
    Fn.add_operator(is_true);

    Fn.add_gate<CNOTMatrix>({ancilla_start, n});

    Fn.add_operator(reverse_operator(is_true));
    Fn.add_operator(reverse_operator(fx));

    return Fn;
}

Operator Grover::check_if_weight_one(int syndrome_start, int ancilla_start) {
    Operator check;
    check.add_gate<CNOTMatrix>({syndrome_start, ancilla_start});
    check.add_gate<CNOTMatrix>({syndrome_start + 1, ancilla_start});
    check.add_gate<CNOTMatrix>({syndrome_start + 2, ancilla_start});

    check.add_gate<CCCNOTMatrix>({syndrome_start, syndrome_start + 1, syndrome_start + 2, ancilla_start});

    return check;
}

Operator Grover::check_if_all_ones(int n) {
    return QuantumUtils::check_if_all_ones(n, 0, n + 1);
    //return QuantumUtils::check_if_all_ones(n, 0 , 22);
}

Operator Grover::minus_one_to_power_of_f(int n, const Operator& Fn) {
    Operator res;
    res.add_gate<HadamardMatrix>({n});
    res.add_operator(Fn);
    res.add_gate<HadamardMatrix>({n});
    return res;
}

Operator Grover::get_f0_operator(int n) {
    Operator f0x;
     for (int i = 0; i < n; i++) 
        f0x.add_gate<NotMatrix>({i});
    
    return f0x;    
}

Operator Grover::get_F0n(int n) {
    Operator f0x = get_f0_operator(n);

    NotMatrix not_matrix;
    Matrix multiCNOT = create_multiple_controlled_matrix(not_matrix, n);
    vector<int> target_qubits;
    target_qubits.reserve(n + 1);
    for (int i = 0; i < n; i++)
        target_qubits.push_back(i);
    target_qubits.push_back(n);

    Operator res;
    res.add_operator(f0x);
    res.add_gate({multiCNOT, target_qubits});
    res.add_operator(reverse_operator(f0x));
    return res;
}

Circuit Grover::grover(const GroverInfo& grover_info) {
    double theta = asin(sqrt(double(grover_info.k) / grover_info.superposition_size));
    int r = M_PI / (4 * theta);

    Circuit grover_circuit(grover_info.n + 1 + grover_info.f_anc);
    // q1 q2 .. qn y f1 f2... fanc
    grover_circuit.add_operator(grover_info.init);
    grover_circuit.add_gate<NotMatrix>({grover_info.n});

    grover_circuit.add_operator(grover_info.superposition);

    for (int i = 0; i < r; i++) {
        //Vf
        grover_circuit.add_operator(minus_one_to_power_of_f(grover_info.n, grover_info.Fn));

        //Hn
        grover_circuit.add_operator(reverse_operator(grover_info.superposition));

        //Rn
        Operator F0n = get_F0n(grover_info.n);
        grover_circuit.add_operator(minus_one_to_power_of_f(grover_info.n, F0n));

        //Hn
        grover_circuit.add_operator(grover_info.superposition);
    }

    grover_circuit.run();
    return grover_circuit;
}

Operator Grover::get_f_for_indexes(int n, vector<int> indexes) {
    cout << "Expected one of the following: ";
    Operator fx;
    for (auto &id: indexes) {
        cout << id << " ";
        Operator f_for_id;
        for (int i = 0; i < n; i++)
            if (!((id >> i) & 1))
                f_for_id.add_gate<NotMatrix>({i});
        fx.add_operator(get_Fn(n, f_for_id));
    }

    cout << "\n";
    return fx;
}

Operator Grover::get_example_function_f() {
    cout << "Expected result: 2\n";
    Operator fx;
    fx.add_gate<SwapMatrix>({1, 3});
    fx.add_gate<NotMatrix>({0});
    fx.add_gate<NotMatrix>({1});

    //checks if output is 11
    fx.add_gate<NotMatrix>({2});

    return get_Fn(4, fx);
}

void Grover::run_grover(GroverInfo grover_info, int it)  {
    Circuit grover_circuit = grover(grover_info);

    map<int, int> M;
    for (int i = 0; i < it; i++) {
        vector<bool> read = grover_circuit.magic_read();
        int res = 0;
        for (int i = 0; i < grover_info.n; i++)
            res += (read[i] << i);
        
        M[res]++;
    }

    for (auto &[res, cnt]: M) 
        cout << res << ": " << double(cnt) / it * 100 << "%\n";
}

void Grover::run_examples() {
    cout << "===================== GROVER =====================\n";

    int n = 4;
    Operator basic_superposition_op;
    for (int i = 0; i < n; i++)
        basic_superposition_op.add_gate<HadamardMatrix>({i});

    GroverInfo grover_info {4, 1, 3, (1 << n), get_f_for_indexes(4, {11}), basic_superposition_op};
    run_grover(grover_info);
    cout << "\n";

    GroverInfo grover_info2 {4, 3, 3, (1 << n), get_f_for_indexes(4, {5, 8, 1}), basic_superposition_op};
    run_grover(grover_info2);
    cout << "\n";

    GroverInfo grover_info3 {4, 1, 3, (1 << n), get_example_function_f(), basic_superposition_op};
    run_grover(grover_info3);
    cout << "\n";
}