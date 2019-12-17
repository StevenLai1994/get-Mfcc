#include <cstdio>
#include "Mfcc.h"
#include "Mfcc.cpp"

int main() {
    Mfcc m = Mfcc();
    m.init();

    clock_t s_clock = clock();
    m.mfcc_from_file(m.test_path);
    clock_t e_clock = clock();
    printf("time cose %.4f ms\n", get_time(e_clock, s_clock));;
    return 0;
}