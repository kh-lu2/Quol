#include "quantum_utils.h"
struct GroverInfo {
    int n;
    int k;
    int f_anc;

    Operator Fn;
};

struct Grover {
private:
    Operator get_Fn(int n, const Operator& fx);
    Operator check_if_all_ones(int n);
    Operator minus_one_to_power_of_f(int n, const Operator&);

    Operator get_f0_operator(int n);
    Operator get_F0n(int n);

    int grover(const GroverInfo& grover_info);
    Operator get_f_for_indexes(int n, vector<int> indexes);
    Operator get_example_function_f();

public:
    void run_grover(GroverInfo grover_info, int it = 50);
    void run_examples();
};