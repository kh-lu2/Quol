#define _USE_MATH_DEFINES
#include <map>
#include <iostream>
#include "../inc/grover.h"

Operator Grover::get_f0_operator(int n) {
    Operator f0x;
     for (int i = 0; i < n; i++) 
        f0x.add_gate<NotGate>({i});
    
    return f0x;    
}

Operator Grover::check_if_one(int n) {
    Operator op;
    op.add_gate<Gate>({0, 1, n}, CCNOT_matrix);
    for (int i = 0; i < n - 2; i++)
        op.add_gate<Gate>({n + i, i + 2, n + i + 1 }, CCNOT_matrix);

    return op;
}


Operator Grover::get_Fn(int n, const Operator& fx) {
    Operator Fn;

    Fn.add_operator(fx);
    Operator if_one_op = check_if_one(n);
    Fn.add_operator(if_one_op);

    Fn.add_gate<Gate>({2 * n - 2, 2 * n - 1}, CNOT_matrix);

    Fn.add_operator(reverse_operator(if_one_op));
    Fn.add_operator(reverse_operator(fx));

    return Fn;
}


Operator Grover::wrap_Fn(int n, const Operator& Fn) {
    Operator Vf;
    Vf.add_gate<HadamardGate>({2 * n - 1});

    Vf.add_operator(Fn);
    Vf.add_gate<HadamardGate>({2 * n - 1});

    return Vf;
}

Operator Grover::get_F0n(int n) {
    Operator f0x = get_f0_operator(n);
    return get_Fn(n, f0x);
}

int Grover::grover(const GroverInfo& grover_info) {
    double theta = asin(sqrt(double(grover_info.k) / (1 << grover_info.n)));
    int r = M_PI / (4 * theta);

    Circuit grover_circuit(2 * grover_info.n + grover_info.f_anc);
    // q1 q2 .. qn a1 a2 ... an-1 y f1 f2... fanc
    grover_circuit.add_gate<NotGate>({2 * grover_info.n - 1});

    for (int i = 0; i < grover_info.n; i++)
        grover_circuit.add_gate<HadamardGate>({i});

    for (int i = 0; i < r; i++) {
        //Vf
        grover_circuit.add_operator(wrap_Fn(grover_info.n, grover_info.Fn));

        //Hn
        for (int i = 0; i < grover_info.n; i++)
            grover_circuit.add_gate<HadamardGate>({i});

        //Rn
        Operator F0n = get_F0n(grover_info.n);
        grover_circuit.add_operator(wrap_Fn(grover_info.n, F0n));

        //Hn
        for (int i = 0; i < grover_info.n; i++)
            grover_circuit.add_gate<HadamardGate>({i});
    }

    grover_circuit.run();
    vector<bool> read = grover_circuit.magic_read();

    int res = 0;
    for (int i = 0; i < grover_info.n; i++)
        res += (read[i] << i);

    return res;
}

Operator Grover::get_Fn_for_indexes(int n, vector<int> indexes) {
    cout << "Expected one of the following: ";
    Operator fx;
    for (auto &id: indexes) {
        cout << id << " ";
        Operator f_for_id;
        for (int i = 0; i < n; i++)
            if (!((id >> i) & 1))
                f_for_id.add_gate<NotGate>({i});
        fx.add_operator(get_Fn(n, f_for_id));
    }

    cout << "\n";
    return fx;
}

Operator Grover::get_example_function_Fn() {
    cout << "Expected result: 2\n";
    Operator fx;
    fx.add_gate<SwapGate>({1, 3});
    fx.add_gate<NotGate>({0});
    fx.add_gate<NotGate>({1});

    //checks if output is 11
    fx.add_gate<NotGate>({2});

    return get_Fn(4, fx);
}

void Grover::run_grover(GroverInfo grover_info, int it)  {
    map<int, int> M;
    for (int i = 0; i < it; i++)
        M[grover(grover_info)]++;

    for (auto &[res, cnt]: M) 
        cout << res << ": " << double(cnt) / it * 100 << "%\n";
}

void Grover::run_examples() {
    GroverInfo grover_info {4, 1, 0, get_Fn_for_indexes(4, {11})};
    run_grover(grover_info);
    cout << "\n";

    GroverInfo grover_info2 {4, 3, 0, get_Fn_for_indexes(4, {5, 8, 1})};
    run_grover(grover_info2);
    cout << "\n";

    GroverInfo grover_info3 {4, 1, 0, get_example_function_Fn()};
    run_grover(grover_info3);
    cout << "\n";
}