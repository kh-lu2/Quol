#include <random>
#include <bit>
#include <iostream>
#include "../inc/quantum_utils.h"

mt19937 rng(random_device{}());
uniform_real_distribution<double> uniform_dist(0.0, 1.0);

Matrix::Matrix(const MCD& matrix) : matrix(matrix) {}
Matrix Matrix::get_reverse() const {
    if (matrix.empty())
        return {MCD{}};

    size_t rows = matrix.size();
    size_t cols = matrix[0].size();
    MCD reversed(cols, vector<CD>(rows, 0));
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            reversed[j][i] = std::conj(matrix[i][j]);
    return {reversed};
}

NotMatrix::NotMatrix() : Matrix (
    {
        {0, 1},
        {1, 0}
    }
) {}

HadamardMatrix::HadamardMatrix() : Matrix (
    {
        {1 / sqrt(2), 1 / sqrt(2)},
        {1 / sqrt(2), -1 / sqrt(2)}
    }
) {}

PhaseFlipMatrix::PhaseFlipMatrix(double theta) : Matrix (
    {
        {1, 0},
        {0, CD(polar(1.0, theta))}
    }
) {}

RotationMatrix::RotationMatrix(double theta) : Matrix (
    {
        {cos(theta / 2), -sin(theta / 2)},
        {sin(theta / 2), cos(theta / 2)}
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

Matrix create_controlled_matrix(const Matrix& matrix) {
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

Matrix create_multiple_controlled_matrix(const Matrix& matrix, int n) {
    Matrix oldMatrix = matrix, newMatrix = matrix;
    for (int i = 0; i < n; i++) {
        newMatrix = create_controlled_matrix(oldMatrix);
        oldMatrix = newMatrix;
    }

    return newMatrix;
}

Gate::Gate(const Matrix& matrix, const vector<int>& t) : matrix(matrix), target_indices(t) {}

Gate create_controlled_gate(const Gate& gate, int control_index) {
    Matrix controlled_matrix = create_controlled_matrix(gate.matrix);
    vector<int> controlled_indexes = {control_index};
    for (auto &idx : gate.target_indices)
        controlled_indexes.push_back(idx);
    return {controlled_matrix, controlled_indexes};
}

void Operator::add_gate(const Gate& gate) {
    gates.push_back(gate);
}

void Operator::add_operator(const Operator& op) {
    for (int i = 0; i < op.gates.size(); i++)
        gates.push_back(op.gates[i]);
}

Operator reverse_operator(const Operator& op) {
    Operator result;
    for (int i = op.gates.size() - 1; i >= 0; i--)
        result.gates.push_back({op.gates[i].matrix.get_reverse(), op.gates[i].target_indices});
    return result;
}

Circuit::Circuit(int n_qubits) : num_qubits(n_qubits) {
    size_t dim = 1ULL << num_qubits;
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

    size_t num_spectators = num_qubits - targets.size(); 
    size_t target_mask = 0;
    for (auto &target: targets)
        target_mask |= (1ULL << target);

    size_t spectator_mask = 0;
    for (size_t i = 0; i < num_qubits; i++)
        if (!((target_mask >> i) & 1))
            spectator_mask |= (1ULL << i);

    for (size_t spectator_state = 0; spectator_state < (1ULL << num_spectators); spectator_state++) {
        StateVector working_state;
        vector<size_t> working_state_ids;
        size_t long_spectator_state = 0;
        size_t spect_id = 0;

        for (size_t i = 0; i < num_qubits; i++) {
            if ((spectator_mask >> i) & 1) {
                long_spectator_state |= ((spectator_state >> spect_id) & 1ULL) << i;
                spect_id++;
            }
        }

        for (size_t target_state = 0; target_state < (1 << targets.size()); target_state++) {
            size_t state_id = long_spectator_state;
            for (size_t target_id = 0; target_id < targets.size(); target_id++)
                if ((target_state >> target_id) & 1)
                    state_id |= (1 << targets[target_id]);
            working_state_ids.push_back(state_id);
            working_state.push_back(state[state_id]);
        }

        StateVector result = get_result_state(working_state, gate.matrix);
        for (size_t i = 0; i < working_state_ids.size(); i++)
            state[working_state_ids[i]] = result[i];
    }
}

void Circuit::run() {
    for (auto& gate : gates)
        apply_gate(gate);
}

vector<bool> Circuit::magic_read() {
    vector<bool> res(num_qubits);
    size_t full_dim = 1 << num_qubits;
    
    vector<double> cumulative_prob(full_dim);
    cumulative_prob[0] = norm(state[0]);
    for (size_t i = 1; i < full_dim; i++)
        cumulative_prob[i] = cumulative_prob[i-1] + norm(state[i]);
    
    double r = uniform_dist(rng);
    size_t measured_state = 0;
    for (size_t i = 0; i < full_dim; i++) {
        if (r <= cumulative_prob[i]) {
            measured_state = i;
            break;
        }
    }
    
    for (size_t q = 0; q < num_qubits; q++)
        res[q] = (measured_state >> q) & 1;
    
    return res;
}

StateVector get_result_state(const StateVector& state, const Matrix& matrix) {
    StateVector result(state.size());
    for (size_t i = 0; i < state.size(); i++)
        for (size_t j = 0; j < state.size(); j++) 
            result[i] += state[j] * matrix.matrix[j][i];

    return result;
}

Operator QuantumUtils::check_if_all_ones(int n, int start_index, int ancilla_start_index) {
    Operator op;
    if (n == 1) return op;

    op.add_gate<CCNOTMatrix>({start_index, start_index + 1, ancilla_start_index});
    for (int i = 0; i < n - 2; i++)
        op.add_gate<CCNOTMatrix>({ancilla_start_index + i, i + 2, ancilla_start_index + i + 1 }); 
    return op;
}

Operator QuantumUtils::get_NOTs(int bits, int number, int start_index, bool compare) {
    Operator op;

    for (int bit = 0; bit < bits; bit++)
        if (((number >> bit) & 1) == compare)
            op.add_gate<NotMatrix>({start_index + bit});

    return op;
}

Operator QuantumUtils::get_minus_one_operator(int bits_count, int start_index) {
    Operator op;
    op.add_gate<NotMatrix>({start_index});

    for (int i = 1; i < bits_count; i++) {
        NotMatrix not_matrix = NotMatrix();
        Matrix multiCNOT = create_multiple_controlled_matrix(not_matrix, i);

        vector<int> target_qubits;
        for (int j = start_index; j <= start_index + i; j++)
            target_qubits.push_back(j);

        op.add_gate({multiCNOT, target_qubits});
    }

    return op;
}