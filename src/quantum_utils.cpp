#include <random>
#include "../inc/quantum_utils.h"

mt19937 rng(random_device{}());
uniform_real_distribution<double> uniform_dist(0.0, 1.0);

Matrix::Matrix(const MCD& matrix) : matrix(matrix) {}

NotMatrix::NotMatrix() : Matrix (
    {
        {0, 1},
        {1, 0}
    }
) {}

HadamardMatrix::HadamardMatrix() : Matrix (
    {
        {1/sqrt(2), 1/sqrt(2)},
        {1/sqrt(2), -1/sqrt(2)}
    }
) {}

PhaseFlipMatrix::PhaseFlipMatrix(double theta) : Matrix (
    {
        {1, 0},
        {0, CD(polar(1.0, theta))}
    }
) {}

SwapMatrix::SwapMatrix() : Matrix (
    {
        {1, 0, 0, 0},
        {0, 0, 1, 0},
        {0, 1, 0, 0},
        {0, 0, 0, 1}
    }
) {}

Matrix createControlledMatrix(const Matrix& matrix) {
    const MCD& U = matrix.matrix;
    int k = U.size();
    MCD controlled_matrix(2 * k);
    for (int i = 0; i < k; i++)
        for (int j = 0; j < 2 * k; j++)
            controlled_matrix[2 * i].push_back((2 * i == j) ? 1 : 0);
    for (int i = 0; i < k; i++)
        for (int j = 0; j < 2 * k; j++)
            controlled_matrix[2 * i + 1].push_back((j % 2 == 0) ? 0 : U[i][j / 2]);
    return {controlled_matrix};
}

Matrix createMultipleControlledMatrix(const Matrix& matrix, int n) {
    Matrix oldMatrix = matrix, newMatrix = matrix;
    for (int i = 0; i < n; i++) {
        newMatrix = createControlledMatrix(oldMatrix);
        oldMatrix = newMatrix;
    }

    return newMatrix;
}

Gate::Gate(const Matrix& matrix, const vector<int>& t) : matrix(matrix), target_indices(t) {}

Gate createControlledGate(const Gate& gate, int control_index) {
    Matrix controlled_matrix = createControlledMatrix(gate.matrix);
    vector<int> controlled_indexes = {control_index};
    for (auto &idx : gate.target_indices)
        controlled_indexes.push_back(idx);
    return {controlled_matrix, controlled_indexes};
}

void Operator::add_operator(const Operator& op) {
    for (int i = 0; i < op.gates.size(); i++)
        gates.push_back(op.gates[i]);
}

Operator reverse_operator(const Operator& op) {
    Operator result;
    for (int i = op.gates.size() - 1; i >= 0; i--)
        result.gates.push_back(op.gates[i]);
    return result;
}

Circuit::Circuit(int n_qubits) : num_qubits(n_qubits) {
    int dim = 1 << num_qubits;
    state = StateVector(dim, 0);
    state[0] = 1;
}

Circuit::~Circuit() {}

void Circuit::add_gate(const Gate& gate) {
    gates.push_back(gate);
}

void Circuit::add_operator(const Operator& op) {
    for (int g = 0; g < op.gates.size(); g++)
        gates.push_back(op.gates[g]);
}

void Circuit::apply_gate(const Gate& gate) {
    auto &targets = gate.target_indices;

    int num_spectators = num_qubits - targets.size(); 
    int target_mask = 0;
    for (auto &target: targets)
        target_mask |= (1 << target);

    int spectator_mask = 0;
    for (int i = 0; i < num_qubits; i++)
        if (!((target_mask >> i) & 1))
            spectator_mask |= (1 << i);

    for (int spectator_state = 0; spectator_state < (1 << num_spectators); spectator_state++) {
        StateVector working_state;
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

        StateVector result = get_result_state(working_state, gate.matrix);
        for (int i = 0; i < working_state_ids.size(); i++)
            state[working_state_ids[i]] = result[i];
    }
}

void Circuit::run() {
    for (auto& gate : gates)
        apply_gate(gate);
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

StateVector get_result_state(const StateVector& state, const Matrix& matrix) {
    StateVector result(state.size());
    for (int i = 0; i < state.size(); i++) {
        for (int j = 0; j < state.size(); j++) {
            result[i] += state[j] * matrix.matrix[j][i];
        }
    }
    return result;
}
