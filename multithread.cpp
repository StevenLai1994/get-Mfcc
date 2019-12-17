#include <thread>
#include <ctime>
#include <cstdio>
using namespace std;
#define get_time(e, s) ((e - s) * 1000.0 / CLOCKS_PER_SEC)
const int num_cpu = 4;
typedef unsigned long int64;
typedef short int int16;

void compute(float* arr, int64 begin, int64 end, int64& ans) {
    ans = 0;
    for(int64 i=begin; i<end; i++) {
        ans = (ans + i) % 10000;
    }
}

void multi_compute(float* arr, long size) {
    //=============1=============
    clock_t s_clock = clock();
    clock_t e_clock;
    int64 avg = size / num_cpu;
    long remain = size % num_cpu;
    int64 shift = 0;
    int64 anses[4] = {0, 0, 0, 0};
    std::thread* threads = new std::thread[num_cpu-1];
    for(int16 i=0; i<num_cpu-1; i++) {
        int64 batch_size = (remain > 0)?(avg+1):avg;
        e_clock = clock();
        printf("线程%d开始时间%.4f ms\n", i, get_time(e_clock, s_clock));
        threads[i] = std::thread(compute, arr, shift, shift+batch_size, ref(anses[i]));
        shift += batch_size;
        --remain;
    }

    compute(arr, shift, size-1, anses[4]);
    int64 ans = 0;
    for(int i=0; i<num_cpu-1; i++) {
        threads[i].join();
        ans = (ans + anses[i]) % 10000;
    }
}

int main() {
    clock_t s_clock, e_clock;
    const int64 size = 1024 * 1024 * 1024;

    float* arr = new float[size];
    s_clock = clock();
    int64 ans = 0;
    for(int i=0; i<10; i++)
        compute(arr, 0, size, ans);
    e_clock = clock();
    printf("ans = %d, 单线程耗时%.4f ms\n", ans, get_time(e_clock, s_clock));
    delete[] arr;

    arr = new float[size];
    s_clock = clock();
    for(int i=0; i<10; i++)
        multi_compute(arr, size);
    e_clock = clock();
    printf("ans = %d, 多线程耗时%.4f ms\n", ans, get_time(e_clock, s_clock));
    delete[] arr;
    return 0;
}