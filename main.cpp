#include "inc/shor.h"
#include "inc/grover.h"

int main() {
    Shor shor;
    shor.run_default();

    Grover grover;
    grover.run_examples();
}
