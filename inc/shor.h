#include "quantum_utils.h"

struct Shor {
private:
    int phase_estimation();
    int get_r_from_fraction_expansion(int p, int n, int a, int m);

public:
    void run_default(int n = 15, int a = 7, int m = 16);
};