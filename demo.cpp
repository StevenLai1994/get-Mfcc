#include <cstdio>
#include "Mfcc.h"
#include "Mfcc.cpp"

int main() {
    Mfcc m = Mfcc();
    m.init();

    timeval s_tv, e_tv;
    set_time(s_tv);
    for(int i=0; i<1; i++) {
        float** features = m.mfcc_from_file(m.test_path);
        delete_2ddata((void**)features, 99);
    }
    set_time(e_tv);
    printf("all time cose %.4f ms\n", get_time(e_tv, s_tv));;
    return 0;
}