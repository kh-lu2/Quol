#include "quantum_utils.h"
#include <array>

struct Deutsch {
    array<Operator, 4> Uf;

    Deutsch();

    bool is_balanced(Operator U);
    void run_deutsch();
};