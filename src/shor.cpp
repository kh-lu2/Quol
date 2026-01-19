#define _USE_MATH_DEFINES
#include <iostream>
#include "../inc/shor.h"

int Shor::phase_estimation()  {
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
                Gate controlled_gate = createControlledGate(U.gates[g]);
                vector<int> controlled_indexes = {i};
                for (auto &idx : U.gate_indexes[g])
                    controlled_indexes.push_back(idx);
                sus.add_gate<Gate>(controlled_indexes, controlled_gate.matrix);
            }
        }
    }

    for (int qubit = 3; qubit >= 0; qubit--) {
        sus.add_gate<HadamardGate>({qubit});
        
        for (int qubit2 = qubit - 1; qubit2 >= 0; qubit2--) {
            double theta = M_PI / (1 << (qubit - qubit2));
            Gate phase_gate({{1, 0}, {0, CD(polar(1.0, theta))}});
            Gate controlled_phase_gate = createControlledGate(phase_gate);
            sus.add_gate<Gate>({qubit2, qubit}, controlled_phase_gate.matrix);
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

int Shor::get_r_from_fraction_expansion(int p, int n, int a, int m) {
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

void Shor::run_default(int n, int a, int m) {
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
    cout << "\n";
}
