#pragma once
#include <complex>
#include <vector>
#include <random>
using namespace std;

typedef complex<double> CD;
typedef vector<CD> StateVector;
typedef vector<vector<CD>> MCD;

extern mt19937 rng;
extern uniform_real_distribution<double> uniform_dist;

struct Matrix {
    MCD matrix;
    Matrix(const MCD& matrix);
    Matrix get_reverse() const;
};

struct NotMatrix : Matrix {
    NotMatrix();
};

struct HadamardMatrix : Matrix {
    HadamardMatrix();
};

struct PhaseFlipMatrix : Matrix
{
    PhaseFlipMatrix(double theta);
};

struct RotationMatrix : Matrix {
    RotationMatrix(double theta);
};

struct SwapMatrix : Matrix {
    SwapMatrix();
};

Matrix create_controlled_matrix(const Matrix& matrix);

struct CNOTMatrix : Matrix {
    CNOTMatrix() : Matrix(create_controlled_matrix(Matrix({{0, 1}, {1, 0}})).matrix) {}
};

struct CCNOTMatrix : Matrix {
    CCNOTMatrix() : Matrix(create_controlled_matrix(create_controlled_matrix(Matrix({{0, 1}, {1, 0}}))).matrix) {}
};

struct CCCNOTMatrix : Matrix {
    CCCNOTMatrix() : Matrix(create_controlled_matrix(create_controlled_matrix(create_controlled_matrix(Matrix({{0, 1}, {1, 0}})))).matrix) {}
};

Matrix create_multiple_controlled_matrix(const Matrix& matrix, int n);

struct Gate {
    Matrix matrix;
    vector<int> target_indices;

    Gate(const Matrix& matrix, const vector<int>& t);
};

Gate create_controlled_gate(const Gate& gate, int control_index);

struct Operator {
    vector<Gate> gates;

    template<typename T, typename... Args>
    void add_gate(vector<int> indexes, Args&&... args) {
        T matrix{forward<Args>(args)...};
        gates.push_back({matrix, indexes});
    }

    void add_gate(const Gate& gate);
    void add_operator(const Operator& op);
};

Operator reverse_operator(const Operator& op);

struct Circuit {
    int num_qubits;
    StateVector state;
    vector<Gate> gates;

    Circuit(int n_qubits);
    ~Circuit();

    template<typename T, typename... Args>
    void add_gate(vector<int> indexes, Args&&... args) {
        T matrix{forward<Args>(args)...};
        gates.push_back({matrix, indexes});
    }

    void add_gate(const Gate& gate);
    void add_operator(const Operator& op);
    void apply_gate(const Gate& gate);
    void run();
    vector<bool> magic_read();
};

StateVector get_result_state(const StateVector& state, const Matrix& matrix);

struct QuantumUtils {
    static Operator check_if_all_ones(int n, int start_index, int ancilla_start_index);
    static Operator get_NOTs(int bits, int number, int start_index, bool compare);
    static Operator get_minus_one_operator(int bits_count, int start_index);
};
