#ifndef QUANTUM_UTILS_H
#define QUANTUM_UTILS_H
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

struct SwapMatrix : Matrix {
    SwapMatrix();
};

Matrix createControlledMatrix(const Matrix& matrix);

struct CNOTMatrix : Matrix {
    CNOTMatrix() : Matrix(createControlledMatrix(Matrix({{0, 1}, {1, 0}})).matrix) {}
};

struct CCNOTMatrix : Matrix {
    CCNOTMatrix() : Matrix(createControlledMatrix(createControlledMatrix(Matrix({{0, 1}, {1, 0}}))).matrix) {}
};

Matrix addMultiCNOTs(const Matrix& matrix, int n);

struct Gate {
    Matrix matrix;
    vector<int> target_indices;

    Gate(const Matrix& matrix, const vector<int>& t);
};

Gate createControlledGate(const Gate& gate, int control_index);

struct Operator {
    vector<Gate> gates;

    template<typename T, typename... Args>
    void add_gate(vector<int> indexes, Args&&... args) {
        T matrix{forward<Args>(args)...};
        gates.push_back({matrix, indexes});
    }

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

#endif // QUANTUM_UTILS_H