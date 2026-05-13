#pragma once
#include "quantum_utils.h"
#include "grover.h"

using namespace std;

struct ISDInfo {
    int n;
    int k;

    vector<int> good_columns_choices;
    Operator matrix_init;
};


struct ISD {
private:
    Operator sus(int dim, int start_index, int all_columns);
    void classic_smol_brut();
    int classic_stupid(vector<bool>B, vector<bool> M, int n, int k);
    void apply_classic_sus(int dim, vector<bool> &M, vector<bool> &X, int all_columns);
    void apply_classic_optimisation_operator(vector<bool>B, vector<bool> &M, vector<bool> &X, int n, int k);
    void brut();
    void classic_brut();
    void classic_big_brut();
    int stupid(vector<int> B, vector<int> M, int n, int k);
    void run_example_system_of_equations();

    Grover grover;
    GroverInfo get_grover_info(ISDInfo isd_info);
    ISDInfo get_small_isd_info();
    ISDInfo get_isd_info();
    ISDInfo get_big_isd_info();
    void run_superposition_solver(ISDInfo isd_info, int it);

    int get_qubit_index(int row, int column, int k, int start_index);
    int get_combinations(int n, int k);
    Operator get_uniform_superposition_operator(int n, int k, int ancilla_start_index);
    Operator controlled_solve_from_back(int n, int k, int start_index);
    Operator get_optimisation_operator(int n, int k, int ancilla_start_index);
public:
    void run_defaults();
};