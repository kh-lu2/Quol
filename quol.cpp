#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include <queue>
#include <cmath>
#include <map>
#include <numeric>
#include <algorithm>
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
    vector<MCD> gates;
    vector<vector<int>> gate_indexes;

    template<typename T, typename... Args>
    void add_gate(vector<int> indexes, Args&&... args) {
        T gate{forward<Args>(args)...};
        gates.push_back(gate.matrix);
        gate_indexes.push_back(indexes);
    }

    void add_operator(const Operator& op) {
        for (int i = 0; i < op.gates.size(); i++) {
            gates.push_back(op.gates[i]);
            gate_indexes.push_back(op.gate_indexes[i]);
        }
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

MCD CNOT_matrix = createControlledMatrix({{0, 1}, {1, 0}});
MCD CCNOT_matrix = createControlledMatrix(CNOT_matrix);

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

    void add_operator(const Operator& op) {
        for (int g = 0; g < op.gates.size(); g++)
            gates.push_back(GateInfo(op.gates[g], op.gate_indexes[g]));
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
        int full_dim = 1 << num_qubits;
        
        vector<double> cumulative_prob(full_dim);
        cumulative_prob[0] = norm(state[0]);
        for (int i = 1; i < full_dim; i++)
            cumulative_prob[i] = cumulative_prob[i-1] + norm(state[i]);
        
        double r = (double) rand() / RAND_MAX;
        int measured_state = 0;
        for (int i = 0; i < full_dim; i++) {
            if (r <= cumulative_prob[i]) {
                measured_state = i;
                break;
            }
        }
        
        for (int q = 0; q < num_qubits; q++)
            res[q] = (measured_state >> q) & 1;
        
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
                MCD controlled_matrix = createControlledMatrix(U.gates[g]);
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
    cout << "Factorizing 15...\n";
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

struct GroverInfo {
    int n;
    int k;
    int f_anc;


    Operator Fn;
};

Operator get_f0_operator(int n) {
    Operator f0x;
     for (int i = 0; i < n; i++) 
        f0x.add_gate<NotGate>({i});
    
    return f0x;    
}

Operator check_if_one(int n) {
    Operator op;
    op.add_gate<Gate>({0, 1, n}, CCNOT_matrix);
    for (int i = 0; i < n - 2; i++)
        op.add_gate<Gate>({n + i, i + 2, n + i + 1 }, CCNOT_matrix);

    return op;
}

Operator reverse_operator(const Operator& op) {
    Operator result;
    for (int i = op.gates.size() - 1; i >= 0; i--) {
        result.gates.push_back(op.gates[i]);
        result.gate_indexes.push_back(op.gate_indexes[i]);
    }
    return result;
}

Operator get_Fn(int n, const Operator& fx) {
    Operator Fn;

    Fn.add_operator(fx);
    Operator if_one_op = check_if_one(n);
    Fn.add_operator(if_one_op);

    Fn.add_gate<Gate>({2 * n - 2, 2 * n - 1}, CNOT_matrix);

    Fn.add_operator(reverse_operator(if_one_op));
    Fn.add_operator(reverse_operator(fx));

    return Fn;
}


Operator wrap_Fn(int n, const Operator& Fn) {
    Operator Vf;
    Vf.add_gate<HadamardGate>({2 * n - 1});

    Vf.add_operator(Fn);
    Vf.add_gate<HadamardGate>({2 * n - 1});

    return Vf;
}

Operator get_F0n(int n) {
    Operator f0x = get_f0_operator(n);
    return get_Fn(n, f0x);
}

int grover(const GroverInfo& grover_info) {
    double theta = asin(sqrt(double(grover_info.k) / (1 << grover_info.n)));
    int r = M_PI / (4 * theta);

    Circuit grover_circuit(2 * grover_info.n + grover_info.f_anc);
    // q1 q2 .. qn a1 a2 ... an-1 y f1 f2... fanc
    grover_circuit.add_gate<NotGate>({2 * grover_info.n - 1});

    for (int i = 0; i < grover_info.n; i++)
        grover_circuit.add_gate<HadamardGate>({i});

    for (int i = 0; i < r; i++) {
        //Vf
        grover_circuit.add_operator(wrap_Fn(grover_info.n, grover_info.Fn));

        //Hn
        for (int i = 0; i < grover_info.n; i++)
            grover_circuit.add_gate<HadamardGate>({i});

        //Rn
        Operator F0n = get_F0n(grover_info.n);
        grover_circuit.add_operator(wrap_Fn(grover_info.n, F0n));

        //Hn
        for (int i = 0; i < grover_info.n; i++)
            grover_circuit.add_gate<HadamardGate>({i});
    }

    grover_circuit.run();
    vector<bool> read = grover_circuit.magic_read();

    int res = 0;
    for (int i = 0; i < grover_info.n; i++)
        res += (read[i] << i);

    return res;
}

Operator get_Fn_for_indexes(int n, vector<int> indexes) {
    cout << "Expected one of the following: ";
    Operator fx;
    for (auto &id: indexes) {
        cout << id << " ";
        Operator f_for_id;
        for (int i = 0; i < n; i++)
            if (!((id >> i) & 1))
                f_for_id.add_gate<NotGate>({i});
        fx.add_operator(get_Fn(n, f_for_id));
    }

    cout << "\n";
    return fx;
}

Operator get_example_function_Fn() {
    cout << "Expected result: 2\n";
    Operator fx;
    fx.add_gate<SwapGate>({1, 3});
    fx.add_gate<NotGate>({0});
    fx.add_gate<NotGate>({1});

    //checks if output is 11
    fx.add_gate<NotGate>({2});

    return get_Fn(4, fx);
}


void run_grover(GroverInfo grover_info, int it = 20)  {
    map<int, int> M;
    for (int i = 0; i < it; i++)
        M[grover(grover_info)]++;

    for (auto &[res, cnt]: M) 
        cout << res << ": " << double(cnt) / it * 100 << "%\n";
}

int main() {
    srand(time(0));
    shor();
    cout << "\n";

    GroverInfo grover_info {4, 1, 0, get_Fn_for_indexes(4, {11})};
    run_grover(grover_info);
    cout << "\n";

    GroverInfo grover_info2 {4, 3, 0, get_Fn_for_indexes(4, {5, 8, 1})};
    run_grover(grover_info2);
    cout << "\n";

    GroverInfo grover_info3 {4, 1, 0, get_example_function_Fn()};
    run_grover(grover_info3);
    cout << "\n";
}
