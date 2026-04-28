#include "quantum_utils.h"

using namespace std;

struct ISD {
private:
    int get_qubit_index(int row, int column, int k, int start_index);
    int get_combinations(int n, int k);
    Operator get_uniform_superposition_operator(int n, int k, int ancilla_start_index);
    int uniform_solver();
    int optimized_uniform_solver();
    void run_example_system_of_equations();
    void run_example_system_of_equations_from_back();
    void run_uniform_solver();
    void run_optimized_uniform_solver();
    void run_toy_ISD();
    Operator swap_chosen_columns_to_front_operator(int n, int k);
    Operator solve(int n, int k, int start_index);
    Operator solve_from_back(int n, int k, int start_index);
    Operator controlled_solve_from_back(int n, int k, int start_index);
    Operator get_example_init();
    Operator get_optimized_example_init();
    Operator get_optimisation_operator(int n, int k, int ancilla_start_index);
public:
    void run_defaults();
};