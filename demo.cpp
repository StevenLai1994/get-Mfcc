#include <cstdio>
#include "Mfcc.h"

int main() {
    Mfcc m = Mfcc();
    m.init();

    m.mfcc_file();
    return 0;
}