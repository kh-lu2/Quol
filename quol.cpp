#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include <queue>
#include <cmath>
#include <map>
#include <numeric>
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
    int k = U.size();

    MCD controlled_matrix(2 * k);

    for (int i = 0; i < k; i++)
        for (int j = 0; j < 2 * k; j++)
            controlled_matrix[2 * i].push_back((2 * i == j) ? 1 : 0);

    for (int i = 0; i < k; i++)
        for (int j = 0; j < 2 * k; j++)
            controlled_matrix[2 * i + 1].push_back((j % 2 == 0) ? 0 : U[i][j / 2]);

    return controlled_matrix;
}

VCD get_result_state (VCD state, MCD matrix) {
    VCD result(state.size());
    for (int i = 0; i < state.size(); i++) {
        for (int j = 0; j < state.size(); j++) {
            result[i] += state[j] * matrix[j][i];
        }
    }
    return result;
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
        int dim = 1 << num_qubits;
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
        for (int g = 0; g < op.gates.size(); g++)
            gates.push_back(GateInfo(op.gates[g]->matrix, op.gate_indexes[g]));
    }

    void apply_gate(const MCD& U, const vector<int>& targets) {
        int num_spectators = num_qubits - targets.size(); 
        int target_mask = 0;
        for (auto &target: targets)
            target_mask |= (1 << target);
        
        int spectator_mask = 0;
        for (int i = 0; i < num_qubits; i++)
            if (!((target_mask >> i) & 1))
                spectator_mask |= (1 << i);

        for (int spectator_state = 0; spectator_state < (1 << num_spectators); spectator_state++) {
            VCD working_state;
            vector<int> working_state_ids;
            int long_spectator_state = 0;
            int spect_id = 0;
            for (int i = 0; i < num_qubits; i++) {
                if ((spectator_mask >> i) & 1) {
                    long_spectator_state |= ((spectator_state >> spect_id) & 1) << i;
                    spect_id++;
                }
            }

            for (int target_state = 0; target_state < (1 << targets.size()); target_state++) {
                int state_id = long_spectator_state;

                for (int target_id = 0; target_id < targets.size(); target_id++)
                    if ((target_state >> target_id) & 1)
                        state_id |= (1 << targets[target_id]);

                working_state_ids.push_back(state_id);
                working_state.push_back(state[state_id]);
            }

            VCD result = get_result_state(working_state, U);

            for (int i = 0; i < working_state_ids.size(); i++)
                state[working_state_ids[i]] = result[i];
        }
    }

    void run() {
        for (auto& gate : gates)
            apply_gate(gate.matrix, gate.target_indices);
    }

    vector<bool> magic_read() {
        vector<bool> res(num_qubits);
        
        for (int q = 0; q < num_qubits; q++) {
            double p1 = 0;
            int full_dim = 1 << num_qubits;
            
            for (int i = 0; i < full_dim; i++)
                if ((i >> q) & 1)
                    p1 += norm(state[i]);
            
            double r = (double) rand() / RAND_MAX;
            res[q] = (r < p1);
        }
        
        return res;
    }
};

int phase_estimation()  {
    Circuit sus(8);
    sus.add_gate<NotGate>({4});

    for (int qubit = 0; qubit < 4; qubit++)
        sus.add_gate<HadamardGate>({qubit});

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
            for (int g = 0; g < U.gates.size(); g++) {
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
            MCD phase_matrix = {{1, 0}, {0, CD(polar(1.0, theta))}};
            MCD controlled_phase = createControlledMatrix(phase_matrix);
            sus.add_gate<Gate>({qubit2, qubit}, controlled_phase);
        }
    }
    
    sus.add_gate<SwapGate>({0, 3});
    sus.add_gate<SwapGate>({1, 2});

    sus.run();
    vector<bool> read = sus.magic_read();
    int res = 0;
    for (int i = 0; i < 4; i++)
        res += (read[i] << i);

    return res;       
}

int get_r_from_fraction_expansion(int p, int n, int a, int m) {
    int temp_m = m, temp_p = p;
    vector<int> wholes {0}; 

    while (temp_p > 0) {
        wholes.push_back(temp_m / temp_p);
        int mod = temp_m % temp_p;
        temp_m = temp_p;
        temp_p = mod;
    }

    int q_prev = 0, q_prevprev = 1;

    for (int i = 0; i < wholes.size(); i++) {
        int r = wholes[i] * q_prev + q_prevprev;
        q_prevprev = q_prev;
        q_prev = r;
        if (! (r % 2) && ((int)pow(a, r / 2) % n != n - 1 ) )
            return r;
    }


    return -1;
}

void shor(int n = 15, int a = 7, int m = 16) {
    bool finished = 0;
    do {
        int p = phase_estimation();
        cout << "got " << p << " from phase estimation\n";
        int r = get_r_from_fraction_expansion(p, n, a, m);

        int q1 = gcd((int)pow(a, r / 2) - 1, n);
        int q2 = gcd((int)pow(a, r / 2) + 1, n);
        if (q1 != 1 && q1 != n) {
            cout << "found divisor! " << q1 << "\n";
            finished = 1;
        }
        if (q2 != 1 && q2 != n) {
            cout << "found divisor! " << q2 << "\n";
            finished = 1;
        }
    } while (!finished);
}


int main() {
    shor();
}