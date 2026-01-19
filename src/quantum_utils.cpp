#include <random>
#include "../inc/quantum_utils.h"

std::mt19937 rng(std::random_device{}());
std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

Gate::Gate(MCD matrix) : matrix(matrix) {};

NotGate::NotGate() : Gate (
    {
        {0, 1},
        {1, 0}
    }
) {}

HadamardGate::HadamardGate() : Gate (
    {
        {1/sqrt(2), 1/sqrt(2)},
        {1/sqrt(2), -1/sqrt(2)}
    }
) {}

PhaseFlipGate::PhaseFlipGate(double theta) : Gate (
    {
        {1, 0},
        {0, CD(polar(1.0, theta))}
    }
) {}

SwapGate::SwapGate() : Gate (
    {
        {1, 0, 0, 0},
        {0, 0, 1, 0},
        {0, 1, 0, 0},
        {0, 0, 0, 1}
    }
) {}


void Operator::add_operator(const Operator& op) {
    for (int i = 0; i < op.gates.size(); i++) {
        gates.push_back(op.gates[i]);
        gate_indexes.push_back(op.gate_indexes[i]);
    }
}

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

GateInfo::GateInfo(MCD m, vector<int> t) : matrix(m), target_indices(t) {}
    
Circuit::Circuit(int n_qubits) : num_qubits(n_qubits) {
    int dim = 1 << num_qubits;
    state = VCD(dim, 0);
    state[0] = 1;
}

Circuit::~Circuit() {}


void Circuit::add_operator(const Operator& op) {
    for (int g = 0; g < op.gates.size(); g++)
        gates.push_back(GateInfo(op.gates[g], op.gate_indexes[g]));
}

void Circuit::apply_gate(const MCD& U, const vector<int>& targets) {
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

void Circuit::run() {
    for (auto& gate : gates)
        apply_gate(gate.matrix, gate.target_indices);
}

vector<bool> Circuit::magic_read() {
    vector<bool> res(num_qubits);
    int full_dim = 1 << num_qubits;
    
    vector<double> cumulative_prob(full_dim);
    cumulative_prob[0] = norm(state[0]);
    for (int i = 1; i < full_dim; i++)
        cumulative_prob[i] = cumulative_prob[i-1] + norm(state[i]);
    
    double r = uniform_dist(rng);
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

Operator reverse_operator(const Operator& op) {
    Operator result;
    for (int i = op.gates.size() - 1; i >= 0; i--) {
        result.gates.push_back(op.gates[i]);
        result.gate_indexes.push_back(op.gate_indexes[i]);
    }
    return result;
}
