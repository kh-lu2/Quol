#ifndef QUANTUM_UTILS_H
#define QUANTUM_UTILS_H
#include <complex>
#include <vector>
#include <random>
using namespace std;

typedef complex<double> CD;
typedef vector<CD> VCD;
typedef vector<vector<CD>> MCD;

extern std::mt19937 rng;
extern std::uniform_real_distribution<double> uniform_dist;

struct Gate {
    MCD matrix;
    Gate(MCD matrix);
};

struct NotGate : Gate {
    NotGate();
};

struct HadamardGate : Gate {
    HadamardGate();
};

struct PhaseFlipGate : Gate
{
    PhaseFlipGate(double theta);
};

struct SwapGate : Gate {
    SwapGate();
};

struct Operator {
    vector<MCD> gates;
    vector<vector<int>> gate_indexes;

    template<typename T, typename... Args>
    void add_gate(vector<int> indexes, Args&&... args) {
        T gate{forward<Args>(args)...};
        gates.push_back(gate.matrix);
        gate_indexes.push_back(indexes);
    }

    void add_operator(const Operator& op);
};

Operator reverse_operator(const Operator& op);

MCD createControlledMatrix(const MCD& U);

extern MCD CNOT_matrix;
extern MCD CCNOT_matrix;

VCD get_result_state (VCD state, MCD matrix);

struct GateInfo {
    MCD matrix;
    vector<int> target_indices;
    
    GateInfo(MCD m, vector<int> t);
};
    
struct Circuit {
    int num_qubits;
    VCD state;
    vector<GateInfo> gates;

    Circuit(int n_qubits);

    ~Circuit();

    template<typename T, typename... Args>
    void add_gate(vector<int> indexes, Args&&... args) {
        auto gate = new T{forward<Args>(args)...};
        gates.push_back(GateInfo(gate->matrix, indexes));
        delete gate;
    }

    void add_operator(const Operator& op);

    void apply_gate(const MCD& U, const vector<int>& targets);

    void run();

    vector<bool> magic_read();
};

#endif // QUANTUM_UTILS_H