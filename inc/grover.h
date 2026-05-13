#pragma once
#include "quantum_utils.h"
struct GroverInfo {
    int n;
    int k;
    int f_anc;
    int superposition_size;

    Operator Fn;
    Operator superposition;
    Operator init;
};

struct Grover {
private:
    static Operator check_if_all_ones(int n);
    static Operator check_if_weight_one(int syndrome_start, int ancilla_start);
    static Operator minus_one_to_power_of_f(int n, const Operator&);

    static Operator get_f0_operator(int n);
    static Operator get_F0n(int n);

    static Circuit grover(const GroverInfo& grover_info);
    static Operator get_f_for_indexes(int n, vector<int> indexes);
    static Operator get_example_function_f();
    static Operator get_Fn(int n, const Operator& fx);

public:
    static Operator get_Fn_for_ISD(int n, const Operator& fx, int syndrome_start, int ancilla_start);
    static void run_grover(GroverInfo grover_info, int it = 1000);
    void run_examples();
};