#include "inc/shor.h"
#include "inc/grover.h"
#include "inc/deutsch.h"

int main() {
    Deutsch deutsch;
    deutsch.run_deutsch();

    Shor shor;
    shor.run_default();

    Grover grover;
    grover.run_examples();
}
