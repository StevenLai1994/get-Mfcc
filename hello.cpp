#include <iostream>
#include <ctime>
#include "Mfcc.cpp"

bool not_equel(float* d1, float* d2, int len) {
    for(int i=0; i<len; i++) {
        if(abs(d1[i] - d2[i]) > 0.01)
            return true;
    }
    return false;
}

// bool not_equel2(float** d1, float** d2, int h, int w) {
//     for(int i=0; i<h; i++) {
//         for(int j=0; j<w; j++)
//             if(d1[i][j] != d2[i][j])
//                 return true;
//     }
//     return false;
// }

// extern char* get_wav_string(const char* path);
int main() {
    float* test_data = new float[8];
    for(int i=0; i<8; i++) {
        test_data[i] = i*1.0 + 1;
    }

    Mfcc m = Mfcc();
    m.init();
    clock_t start = clock();
    for(int i=0; i<1; i++) {
        float** dct_f = m.mfcc_file(m.test_path);
    }

    // float* ft1 = m.single_dft(test_data, 1024);
    // float* ft2 = m.single_fft(test_data, 1024);
    // not_equel(ft1, ft2, 513);

    clock_t end = clock();
    double timecost = (double)(end-start) / CLOCKS_PER_SEC;
    printf("耗时%.4lf\n", timecost);
    return 0;
}