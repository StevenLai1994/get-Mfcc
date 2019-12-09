#pragma once
#include<cstdio>

typedef unsigned char uchar;
typedef short int16;
typedef char byte;

//2s 16k/s16bit 单声道音频
const float PI = 3.1415926535898;

struct complex {
    float x,y;
    complex(float xx=0, float yy=0):x(xx), y(yy){}
    complex operator + (complex a) {
        return complex(x + a.x, y + a.y);
    }
    complex operator - (complex a) {
        return complex(x - a.x, y - a.y);
    }
    complex operator * (complex a) {
        return complex(x * a.x - y * a.y, x * a.y + y * a.x);
    }
};

class Mfcc {
    public:
        int sample_rate;    //16000
        int voice_samples;   //32000
        int16 window_size;    //640
        int16 step_size;  //320
        int16 time_step;  //99
        float lf, hf;   //20, 4000
        int16 fbank_count;    //40
        int16 dct_count;  //10
        int16 nfft;   //1024
        int16 ceplifter; //22
        const char* test_path;
        float * Hamming;
        
        void show_info();
        void init();

        //functions
        float* get_strings_from_bytes(byte* bytes);    //单声道，双精度， 16K, 一次处理32000个采样点, 64000bytes
        float* get_strings_from_wav(const char* path);    //单声道，双精度， 16K, 一次处理2秒， wav文件44bytes头
        void preemp(float* datas, float preemph=0.97);
        float** frame_seg(float* datas, int16 window_size, int16 step_size);
        float* single_dft(float* frames, int16 nfft); //单帧的dft
        float* single_fft(float* frames, int16 nfft); //单帧的fft
        float mel2hz(float hz);
        float hz2mel(float mel);
        float** get_melfilter(float lf, float hf, int16 fbank_count, int16 nfft);
        float** filte_and_log(float** frames, float** filters, int fbank_count, int nfft);
        float** dct2(float** mfcc_logf, int16 dct_count, int16 fbank_count, int16 ceplifter);
        float** mfcc(float* datas);
        float** mfcc_string(byte* bytes);
        float** mfcc_file(const char* file_path);
};