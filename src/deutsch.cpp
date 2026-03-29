#include "../inc/deutsch.h"
#include <iostream>

Deutsch::Deutsch() {
    Uf[1].add_gate<CNOTMatrix>({0, 1});

    Uf[2].add_gate<NotMatrix>({0});
    Uf[2].add_gate<CNOTMatrix>({0, 1});

    Uf[3].add_gate<NotMatrix>({0});
}

bool Deutsch::is_balanced(Operator U) {
    Circuit deutsch(2);
    deutsch.add_gate<NotMatrix>({1});
    deutsch.add_gate<HadamardMatrix>({0});
    deutsch.add_gate<HadamardMatrix>({1});
    deutsch.add_operator(U);
    deutsch.add_gate<HadamardMatrix>({0});
    deutsch.run();

    return deutsch.magic_read().front();
}

void Deutsch::run_deutsch() {
    cout << "===================== DEUTSCH ====================\n";
    for (int i = 0; i < 4; i++) {
        cout << "Checking if 0 -> " << ((i >> 1) & 1) << ", 1 -> " << (i & 1) << " is balanced\n";
        cout << "Result: " << boolalpha << is_balanced(Uf[i]) << "\n\n";
    }
}