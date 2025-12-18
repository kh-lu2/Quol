#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include <queue>
#include <cmath>
#include <map>
using namespace std;

typedef complex<double> CD;
typedef vector<CD> VCD;
typedef vector<vector<CD>> MCD;

struct Gate {
    MCD matrix;
    Gate(MCD matrix) : matrix(matrix) {};
};

struct NotGate : Gate {
    NotGate() : Gate (
        {
            {0, 1},
            {1, 0}
        }
    ) {}
};

struct HadamardGate : Gate {
    HadamardGate() : Gate (
        {
            {1/sqrt(2), 1/sqrt(2)},
            {1/sqrt(2), -1/sqrt(2)}
        }
    ) {}

};

struct PhaseFlipGate : Gate
{
    PhaseFlipGate(double theta) : Gate (
        {
            {1, 0},
            {0, CD(polar(1.0, theta))}
        }
    ) {}
};

struct SwapGate : Gate {
    SwapGate() : Gate (
        {
            {1, 0, 0, 0},
            {0, 0, 1, 0},
            {0, 1, 0, 0},
            {0, 0, 0, 1}
        }
    ) {}
};

struct Operator {
    vector<Gate*> gates;
    vector<vector<int>> gate_indexes;

    ~Operator() {
        for (auto gate : gates)
            delete gate;
    }

    template<typename T, typename... Args>
    void add_gate(vector<int> indexes, Args&&... args) {
        auto gate = new T{forward<Args>(args)...};
        gates.push_back(gate);
        gate_indexes.push_back(indexes);
    }
};

MCD createControlledMatrix(const MCD& U) {
    size_t target_dim = U.size();
    size_t dim = target_dim * 2;
    MCD controlled_matrix(dim, VCD(dim, 0));
    
    for (size_t r = 0; r < dim; r++) {
        for (size_t c = 0; c < dim; c++) {
            if ((r & 1) == 0) {
                if (r == c)
                    controlled_matrix[r][c] = 1;
            } else {
                size_t r_sub = r >> 1;
                size_t c_sub = c >> 1;
                if ((r & 1) == (c & 1))
                    controlled_matrix[r][c] = U[r_sub][c_sub];
            }
        }
    }
    
    return controlled_matrix;
}

struct GateInfo {
    MCD matrix;
    vector<int> target_indices;
    
    GateInfo(MCD m, vector<int> t) : matrix(m), target_indices(t) {}
};
    
struct Circuit {
    int num_qubits;
    VCD state;
    vector<GateInfo> gates;

    Circuit(int n_qubits) : num_qubits(n_qubits) {
        size_t dim = 1 << num_qubits;
        state = VCD(dim, 0);
        state[0] = 1;
    }

    ~Circuit() {}

    template<typename T, typename... Args>
    void add_gate(vector<int> indexes, Args&&... args) {
        auto gate = new T{forward<Args>(args)...};
        gates.push_back(GateInfo(gate->matrix, indexes));
        delete gate;
    }

    void add_operator(Operator& op) {
        for (size_t g = 0; g < op.gates.size(); g++) {
            gates.push_back(GateInfo(op.gates[g]->matrix, op.gate_indexes[g]));
        }
        op.gates.clear();
        op.gate_indexes.clear();
    }

    void apply_gate(const MCD& gate_matrix, const vector<int>& target_indices) {
        size_t full_dim = 1 << num_qubits;
        VCD new_state(full_dim, 0);
        
        for (size_t i = 0; i < full_dim; i++) {
            for (size_t j = 0; j < full_dim; j++) {
                bool non_target_match = true;
                for (int q = 0; q < num_qubits; q++) {
                    bool is_target = false;
                    for (int t : target_indices) {
                        if (t == q) {
                            is_target = true;
                            break;
                        }
                    }
                    if (!is_target) {
                        bool i_bit = (i >> q) & 1;
                        bool j_bit = (j >> q) & 1;
                        if (i_bit != j_bit) {
                            non_target_match = false;
                            break;
                        }
                    }
                }
                
                if (!non_target_match) continue;
                
                size_t i_target = 0, j_target = 0;
                for (size_t t = 0; t < target_indices.size(); t++) {
                    int qubit = target_indices[t];
                    if ((i >> qubit) & 1) i_target |= (1 << t);
                    if ((j >> qubit) & 1) j_target |= (1 << t);
                }
                
                new_state[i] += gate_matrix[i_target][j_target] * state[j];
            }
        }
        
        state = new_state;
    }

    void run() {
        for (auto& gate : gates) {
            apply_gate(gate.matrix, gate.target_indices);
        }
    }

    vector<bool> read() {
        vector<bool> res(num_qubits);
        
        for (int q = 0; q < num_qubits; q++) {
            double p1 = 0;
            size_t full_dim = 1 << num_qubits;
            
            for (size_t i = 0; i < full_dim; i++) {
                if ((i >> q) & 1) {
                    p1 += norm(state[i]);
                }
            }
            
            double r = (double) rand() / RAND_MAX;
            res[q] = (r < p1);
            
            VCD collapsed_state(full_dim, 0);
            double norm_factor = 0;
            
            for (size_t i = 0; i < full_dim; i++) {
                bool bit = (i >> q) & 1;
                if (bit == res[q]) {
                    collapsed_state[i] = state[i];
                    norm_factor += norm(state[i]);
                }
            }
            
            if (norm_factor > 0) {
                for (size_t i = 0; i < full_dim; i++) {
                    collapsed_state[i] /= sqrt(norm_factor);
                }
            }
            state = collapsed_state;
        }
        
        return res;
    }
};


int main() {
    int z = 30;
    map<int, int> counts;
    while (z--)
    {    
        Circuit sus(8);
        sus.add_gate<NotGate>({4});

        for (int qubit = 0; qubit < 4; qubit++)
        {
            sus.add_gate<HadamardGate>({qubit});
        }

        Operator U;
        U.add_gate<SwapGate>({4, 5});
        U.add_gate<SwapGate>({5, 6});
        U.add_gate<SwapGate>({6, 7});
        U.add_gate<NotGate>({4});
        U.add_gate<NotGate>({5});
        U.add_gate<NotGate>({6});
        U.add_gate<NotGate>({7});

        for (int i = 0; i < 4; i++) {
            int reps = 1 << i;
            for (int r = 0; r < reps; r++) {
                for (size_t g = 0; g < U.gates.size(); g++) {
                    MCD controlled_matrix = createControlledMatrix(U.gates[g]->matrix);
                    vector<int> controlled_indexes = {i};
                    for (auto &idx : U.gate_indexes[g])
                        controlled_indexes.push_back(idx);
                    sus.add_gate<Gate>(controlled_indexes, controlled_matrix);
                }
            }
        }

        for (int qubit = 3; qubit >= 0; qubit--) {
            sus.add_gate<HadamardGate>({qubit});
            
            for (int qubit2 = qubit - 1; qubit2 >= 0; qubit2--) {
                double theta = M_PI / (1 << (qubit - qubit2));
                MCD phaseMatrix = {{1, 0}, {0, CD(polar(1.0, theta))}};
                MCD controlledPhase = createControlledMatrix(phaseMatrix);
                sus.add_gate<Gate>({qubit2, qubit}, controlledPhase);
            }
        }
        
        sus.add_gate<SwapGate>({0, 3});
        sus.add_gate<SwapGate>({1, 2});

        sus.run();
        vector<bool> read = sus.read();
        int res = 0;
        for (int i = 0; i < 4; i++)
            res += (read[i] << i);

        counts[res]++;        
    }

    for (auto& [count, value] : counts) {
        cout << count << ": " << value << endl;
    }

}