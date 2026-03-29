#include "../inc/ISD.h"
#include "../inc/quantum_utils.h"
#include <iostream>

void ISD::run_example_system_of_equations() {
    cout << "======================= ISD ======================\n";

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
                Matrix gate_with_CNOTs = addMultiCNOTs(gate.matrix, j - i);
                vector<int> target_qubits = target_cnot_qubits;
                for (auto const &index: gate.target_indices)
                    target_qubits.push_back(index);
                
                ISD.add_gate({gate_with_CNOTs, target_qubits});
            }

            for (int it = i; it < j; it++)
                ISD.add_gate<NotMatrix>({i * n + it});
        }
    }

    ISD.run();
    vector<bool> read = ISD.magic_read();

    for (int i = 0; i < 12; i++) 
        cout << i << " " << read[i] << "\n"; 

}