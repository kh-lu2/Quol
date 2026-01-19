#include "quantum_utils.h"
struct GroverInfo {
    int n;
    int k;
    int f_anc;

    Operator Fn;
};

struct Grover {
private:
    Operator get_f0_operator(int n);
    Operator check_if_one(int n);
    Operator get_Fn(int n, const Operator& fx);
    Operator wrap_Fn(int n, const Operator& Fn);
    Operator get_F0n(int n);
    int grover(const GroverInfo& grover_info);
    Operator get_Fn_for_indexes(int n, vector<int> indexes);
    Operator get_example_function_Fn();

public:
    void run_grover(GroverInfo grover_info, int it = 50);
    void run_examples();
};