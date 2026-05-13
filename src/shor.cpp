#define _USE_MATH_DEFINES
#include <iostream>
#include <numeric>
#include "../inc/shor.h"

int Shor::phase_estimation()  {
    Circuit shor(8);
    shor.add_gate<NotMatrix>({4});

    for (int qubit = 0; qubit < 4; qubit++)
        shor.add_gate<HadamardMatrix>({qubit}); 
    Operator U;
    U.add_gate<SwapMatrix>({4, 5});
    U.add_gate<SwapMatrix>({5, 6});
    U.add_gate<SwapMatrix>({6, 7});
    U.add_gate<NotMatrix>({4});
    U.add_gate<NotMatrix>({5});
    U.add_gate<NotMatrix>({6});
    U.add_gate<NotMatrix>({7});

    for (int i = 0; i < 4; i++) {
        int reps = 1 << i;
        for (int r = 0; r < reps; r++)
            for (int g = 0; g < U.gates.size(); g++)
                shor.add_gate(create_controlled_gate(U.gates[g], i));
    }

    for (int qubit = 3; qubit >= 0; qubit--) {
        shor.add_gate<HadamardMatrix>({qubit});
        
        for (int qubit2 = qubit - 1; qubit2 >= 0; qubit2--) {
            double theta = M_PI / (1 << (qubit - qubit2));
            shor.add_gate(create_controlled_gate({PhaseFlipMatrix(theta), {qubit}}, qubit2));
        }
    }
    
    shor.add_gate<SwapMatrix>({0, 3});
    shor.add_gate<SwapMatrix>({1, 2});

    shor.run();
    vector<bool> read = shor.magic_read();
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
    cout << "====================== SHOR ======================\n";
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
