#include "../inc/ISD.h"
#include "../inc/quantum_utils.h"
#include <iostream>

void ISD::run_example_system_of_equations() {
    cout << "======================= ISD ======================\n";
    cout << "Solving the following system of binary equations:\n";
    cout << "x2 + x3 = 1\n";
    cout << "x1 + x3 = 0\n";
    cout << "x1 + x2 + x3 = 1\n\n";
    cout << "Expected values: x1 = 0, x2 = 1, x3 = 0\n\n";

    cout << "Input matrix:\n";
    cout << "0 1 1 | 1\n";
    cout << "1 0 1 | 0\n";
    cout << "1 1 1 | 1\n\n";

    Circuit ISD(12);
    ISD.add_gate<NotMatrix>({1});
    ISD.add_gate<NotMatrix>({2});
    ISD.add_gate<NotMatrix>({3});
    ISD.add_gate<NotMatrix>({5});
    ISD.add_gate<NotMatrix>({6});
    ISD.add_gate<NotMatrix>({7});
    ISD.add_gate<NotMatrix>({8});
    ISD.add_gate<NotMatrix>({9});
    ISD.add_gate<NotMatrix>({11});

    int n = 3;

    for (int i = 0; i < n - 1; i++) {
        //pivot search
        for (int j = i + 1; j < n; j++) {            
            Operator add_row_op;
            for (int it = i + 1; it < n + 1; it++)
                add_row_op.add_gate<CNOTMatrix>({it * n + j, it * n + i});

            for (int it = i; it < j; it++)
                ISD.add_gate<NotMatrix>({i * n + it});

            vector<int> target_cnot_qubits;
            for (int it = i; it < j; it++)
                target_cnot_qubits.push_back(i * n + it);

            for (auto const &gate: add_row_op.gates) {
                Matrix gate_with_CNOTs = createMultipleControlledMatrix(gate.matrix, j - i);
                vector<int> target_qubits = target_cnot_qubits;
                for (auto const &index: gate.target_indices)
                    target_qubits.push_back(index);
                
                ISD.add_gate({gate_with_CNOTs, target_qubits});
            }

            for (int it = i; it < j; it++)
                ISD.add_gate<NotMatrix>({i * n + it});
        }

        //row reduce
        for (int j = i + 1; j < n; j++) {
            Operator add_row_op;
            for (int it = i + 1; it < n + 1; it++)
                add_row_op.add_gate<CNOTMatrix>({it * n + i, it * n + j});

            for (auto const &gate: add_row_op.gates)             
                ISD.add_gate(createControlledGate(gate, (i * n + j)));
        }
    }

    //back substitution
    for (int i = n - 1; i >= 1; i--)
        for (int j = i - 1; j >= 0; j--)
            ISD.add_gate<CCNOTMatrix>({i * n + j, n * n + i, n * n + j});


    ISD.run();
    vector<bool> read = ISD.magic_read();

    cout << "Output matrix:\n";
    cout << to_string(read[0]) + " " + to_string(read[3]) + " " + to_string(read[6]) + " | " + to_string(read[9]) + "\n";
    cout << to_string(read[1]) + " " + to_string(read[4]) + " " + to_string(read[7]) + " | " + to_string(read[10]) + "\n";
    cout << to_string(read[2]) + " " + to_string(read[5]) + " " + to_string(read[8]) + " | " + to_string(read[11]) + "\n\n";

    cout << "Calculated values: x1 = " + to_string(read[9]) + ", x2 = " + to_string(read[10]) + ", x3 = " + to_string(read[11]) + "\n\n";
}
