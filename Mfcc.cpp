#pragma once
#include "Mfcc.h"

//utils
void** create_2ddata(uint h, uint w, size_t psize) {  
    void **array2D = new void*[h];  
    for(uint i=0; i<h; ++i)  
    {  
        array2D[i] = (void*)new char[w * psize]();  
    }
    return array2D;
}

void delete_2ddata(void** buff, uint h) {
    for(uint i=0; i<h; ++i)  
    {  
        delete[] ((uchar*)buff[i]);
    }  
    // 然后回收高一级的动态数组.  
    delete[] buff;  
}

void* copy1d(void* datas, size_t size, size_t psize) {
    //把datas复制到指定长度;
    void* copy_datas =  (void*)new char[size * psize];
    memcpy(copy_datas, datas, size * psize);
    return copy_datas;
}

void hadman_mul(float* datas1, float* datas2, int len) {
    for(int i=0; i<len; i++) {
        datas1[i] *= datas2[i];
    }
}

float dot_mul(float* datas1, float* datas2, int len) {
    float ret = 0;
    for(int i=0; i<len; i++) {
        ret += (datas1[i] * datas2[i]);
    }
    return ret;
}

//functions
void Mfcc::show_info() {
    printf("这是一个测试%d", this->voice_samples);
}

void Mfcc::init() {
    this->test_path = "./test_datas/key_word.wav";
    this->sample_rate = 16000;
    this->voice_samples = 32000;
    this->window_size = 640;
    this->step_size = 320;
    this->time_step = (this->voice_samples - this->window_size) / this->step_size + 1;
    this->lf = 20;
    this->hf = 4000;
    this->fbank_count = 40;
    this->dct_count = 10;
    this->nfft = 1024;
    this->ceplifter = 22;
    this->Hamming = new float[this->window_size];
    for(int i=0; i<window_size; i++) {
        this->Hamming[i] = 0.54 - 0.46 * cos(2 * i * PI / (window_size-1));
    }
}

// functions
float* Mfcc::get_datas_from_string(byte* voice_string) {
    float* ret = new float[voice_samples];
    int16* data = (int16*)voice_string;
    for(int i=0; i<voice_samples; i++) {
        ret[i] = (float)data[i];
    }
    return ret;
}

float* Mfcc::get_datas_from_file(const char* file_path) {
    FILE* fp = fopen(file_path, "rb");
    char* buffer = NULL;
    int bbytes = voice_samples*2 + 44;
    if(fp != NULL) {
        buffer = new char[bbytes];
        fread(buffer, bbytes, 1, fp);
    }
    fclose(fp);
    float* ret = get_datas_from_string(buffer + 44);
    delete[] buffer;
    return ret;
}

void Mfcc::preemp(float* datas, float preemph) {
    //预加重
    float* tmp = (float*)copy1d(datas, voice_samples, sizeof(float));
    for(int i=1; i<voice_samples; i++) {
        datas[i] = tmp[i] - preemph * tmp[i-1];
    }
    delete[] tmp;
}

float** Mfcc::frame_seg(float* datas, int16 window_size, int16 step_size) {
    //分帧、加窗。
    //return shape (99， 640)
    float** ret = (float**)create_2ddata(time_step, window_size, sizeof(float));
    for(int i=0; i<time_step; i++) {
        memcpy(ret[i], datas + (i * step_size), window_size * sizeof(float));
        hadman_mul(ret[i], Hamming, window_size);
    }
    return ret;
}

float* Mfcc::single_dft(float* frames, int16 window_size, int16 nfft) {
    //太慢了没用
    //对分帧加窗后的各帧信号进行DFT变换得到各帧的频谱
	//并对语音信号的频谱取模平方得到语音信号的功率谱
    //return shape (99, 513)
    float* tmp_datas = new float[nfft]();
    memcpy(tmp_datas, frames, window_size*sizeof(float));
    float* ret = new float[nfft / 2 + 1];
    float* data_r = new float[nfft];
    float* data_i = new float[nfft];
    for(int k=0; k<nfft; k++) {
        float f = (2 * PI * k) / nfft;
        for(int n=0; n<nfft; n++) {
            data_r[k] += tmp_datas[n] * cos(n * f);
            data_i[k] -= tmp_datas[n] * sin(n * f);
        }
    }
    for(int i=0; i<nfft / 2 + 1; i++) {
        ret[i] = ((data_r[i]*data_r[i] + data_i[i]*data_i[i])) / nfft;
    }
    return ret;
}

float* Mfcc::single_fft(float* frames, int16 window_size, int16 nfft) {   
    //对单帧的快速傅里叶变换
    //return shape (99, 513)
    float* tmp_datas = new float[nfft];
    memset(tmp_datas, 0, nfft * sizeof(float));
    memcpy(tmp_datas, frames, window_size * sizeof(float));
    float* ret = new float[nfft / 2 + 1];
    uint bit = 1;
    while(nfft >> bit) ++bit;
    --bit;
    uint *rev = new uint[nfft];
    rev[0] = 0;
    for(int i=0; i<nfft; i++)
        rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
    complex* d_fft = new complex[nfft];
    memset(d_fft, 0, nfft * sizeof(complex));
    for(uint i=0; i<nfft; i++) {
        d_fft[i].x = tmp_datas[rev[i]];
    }
    delete[] tmp_datas;
    delete rev;

    for (int mid=1; mid<nfft; mid*=2)//mid是准备合并序列的长度的二分之一
    {
        complex temp(cos(PI/mid),sin(PI/mid));//单位根，pi的系数2已经约掉了
        for (int i=0; i<nfft; i+=mid*2)//mid*2是准备合并序列的长度，i是合并到了哪一位
		{
            complex omega(1,0);
            for (int j=0; j<mid; j++,omega=omega*temp)//扫描一半
            {
                complex x = d_fft[i + j], y = omega*d_fft[i + j + mid];
                d_fft[i + j] = x + y;
                d_fft[i + j + mid] = x - y;//蝴蝶变换
            }
        }
    }
    for(int i=0; i<nfft / 2 + 1; i++) {
        float x = d_fft[i].x;
        float y = d_fft[i].y;
        ret[i] = (x*x + y*y) / nfft;
    }
    delete[] d_fft;
    return ret;
}

float Mfcc::mel2hz(float mel) {
    return 700*(powf(10, mel/2595.0)-1);
}

float Mfcc::hz2mel(float hz) {
    return 2595 * log10(1 + hz/700);
}

float** Mfcc::get_melfilter(float lf, float hf, int sample_rate, int16 fbank_count, int16 nfft) {
    //获取梅尔滤波器
    float mell = hz2mel(lf);
    float melh = hz2mel(hf);
    float seg = (melh - mell) / (fbank_count + 1);
    int16* bin = new int16[fbank_count+2];
    float frac = (nfft + 1.0) / sample_rate;
    float nowf = mell;
    for(int i=0; i<fbank_count+2; i++) {
        bin[i] = floor(frac * mel2hz(nowf));
        nowf += seg;
    }
    float** ret = (float**)create_2ddata(fbank_count, nfft/2+1, sizeof(float));
    for(int j=0; j<fbank_count; j++) {
        for(int i=bin[j]; i<bin[j+1]; i++) {
            ret[j][i] = (i - bin[j])*1.0 / (bin[j+1]-bin[j]);
        }
        for(int i=bin[j+1]; i<bin[j+2]; i++) {
            ret[j][i] = (bin[j+2]-i)*1.0 / (bin[j+2]-bin[j+1]);
        }
    }
    delete[] bin;
    return ret;
}

void Mfcc::filte_and_log(float** frames, float** filters, float** mfcc_logf,
        int16 window_size, int16 time_step, int16 fbank_count, int16 nfft) {
    //return shape (99, 40)
    int16 fft_size = nfft /2 + 1;
    for(int step=0; step<time_step; step++) {
        float* step_fft = single_fft(frames[step], window_size, nfft);   //1 * 513;
        for(int f=0; f < fbank_count; f++) {
            float tmp_val = dot_mul(step_fft, filters[f], fft_size);
            tmp_val = log(tmp_val);
            if(tmp_val == 0)
                tmp_val = 0.00001;
            mfcc_logf[step][f] = tmp_val;
        }
        delete[] step_fft;
    }
    clock_t e_clock = clock();
    // printf("%d个fft, filte, log 多线程耗时%.4f ms\n", time_step, get_time(e_clock, s_clock));
}

void Mfcc::filte_and_log_and_dct2(float** frames, float** filters, float** dct_mfcc_logf, float* lifter,
        int16 window_size, int16 time_step, int16 fbank_count, int16 nfft, int16 dct_count, int16 ceplifter, float f0, float f1) {
    //对timestep个帧进行 fft --> 滤波 --> 取对数 --> dct2
    //return shape (time_step, dct_count)
    int16 fft_size = nfft /2 + 1;
    float* mfcc_logf = new float[fbank_count];      //1 * 40 存放单个时间步的临时mfcc_logf
    for(int step=0; step<time_step; step++) {
        float* step_fft = single_fft(frames[step], window_size, nfft);   //1 * 513;
        for(int f=0; f < fbank_count; f++) {
            float tmp_val = dot_mul(step_fft, filters[f], fft_size);
            tmp_val = log(tmp_val);
            if(tmp_val == 0)
                tmp_val = 0.00001;
            mfcc_logf[f] = tmp_val;
        }
        //dct2
        for(int i=0; i<dct_count; i++) {
            float sum = 0;
            for(int j=0; j<fbank_count; j++) {
                sum += mfcc_logf[j] * cos(PI * i * (2 * j + 1)/(2 * fbank_count));
            }
            if(i == 0) {
                dct_mfcc_logf[step][i] = f0 * 2 * sum;
            }
            else {
                dct_mfcc_logf[step][i] = f1 * 2 * sum;
            }
            dct_mfcc_logf[step][i] *= lifter[i];

        }
        delete[] step_fft;
    }
    delete[] mfcc_logf;
}

void Mfcc::multi_filte_and_log_and_dct2(float** frames, float** filters, float** dct_mfcc_logf,
        int16 window_size, int16 time_step, int16 fbank_count, int16 nfft, int16 dct_count, int16 ceplifter) {
    //多线程 fft + 梅尔滤波 + log + dct2
    //result shape (99, 10)
    int16 avg = time_step / num_cpu;
    int16 remain = time_step % num_cpu;
    int16 shift = 0;

    //make lifter
    float* lifter = new float[dct_count];
    for(int i=0; i<dct_count; i++) {
        if(ceplifter > 0) {
            lifter[i] = 1 + ceplifter / 2.0 * sin(PI * i / ceplifter);
        }
        else {
            lifter[i] = 1;
        }
    }
    //dct2系数
    float f0 = sqrt(1.0 / (4 * fbank_count));
    float f1 = sqrt(1.0 / (2 * fbank_count));

    // =============1=============
    std::thread* threads = new std::thread[num_cpu-1];
    for(int16 i=0; i<num_cpu-1; i++) {
        int16 batch_size = (remain > 0)?(avg+1):avg;
        float** frames_batch = frames + shift;
        float** mfcc_logf_batch = dct_mfcc_logf + shift;
        threads[i] = std::thread(filte_and_log_and_dct2, frames_batch, filters, mfcc_logf_batch, lifter, window_size, batch_size,
                fbank_count, nfft, dct_count, ceplifter, f0, f1);
        shift += batch_size;
        --remain;
    }
    
    int16 batch_size = (remain > 0)?(avg+1):avg;
    float** frames_batch = frames + shift;
    float** mfcc_logf_batch = dct_mfcc_logf + shift;
    filte_and_log_and_dct2(frames_batch, filters, mfcc_logf_batch, lifter, window_size, batch_size,
            fbank_count, nfft, dct_count, ceplifter, f0, f1);
    for(int16 i=0; i<num_cpu-1; i++) {
        threads[i].join();
    }
    delete[] lifter;
    delete[] threads;
}


float** Mfcc::dct2(float** mfcc_logf, int16 dct_count, int16 fbank_count, int16 ceplifter) {
    //离散余弦变换II, ortho
    //return shape (99, 10)
    float* lifter = new float[dct_count];
    float f0 = sqrt(1.0 / (4 * fbank_count));
    float f1 = sqrt(1.0 / (2 * fbank_count));
    float** ret = (float**)create_2ddata(time_step, dct_count, sizeof(float));

    for(int i=0; i<dct_count; i++) {
        if(ceplifter > 0) {
            lifter[i] = 1 + ceplifter / 2.0 * sin(PI * i / ceplifter);
        }
        else {
            lifter[i] = 1;
        }
    }

    for(int step=0; step<time_step; step++) {
        for(int i=0; i<dct_count; i++) {
            float sum = 0;
            for(int j=0; j<fbank_count; j++) {
                sum += mfcc_logf[step][j] * cos(PI * i * (2 * j + 1)/(2 * fbank_count));
            }
            if(i == 0) {
                ret[step][i] = f0 * 2 * sum;
            }
            else {
                ret[step][i] = f1 * 2 * sum;
            }
            ret[step][i] *= lifter[i];

        }
    }
    delete[] lifter;
    return ret;
}

float** Mfcc::mfcc(float* datas) {
    preemp(datas);//预加重
    float** frames = frame_seg(datas, window_size, step_size);
    float** melfilters = get_melfilter(lf, hf, sample_rate, fbank_count, nfft);
    float** dct_mfcc_logf = (float**)create_2ddata(time_step, dct_count, sizeof(float));
    // filte_and_log(frames, melfilters, mfcc_logf, window_size, time_step, fbank_count, nfft);
    multi_filte_and_log_and_dct2(frames, melfilters, dct_mfcc_logf, window_size,
            time_step, fbank_count, nfft, dct_count, ceplifter);
    delete_2ddata((void**)melfilters, fbank_count);
    delete_2ddata((void**)frames, time_step);
    // float** ret = dct2(mfcc_logf, dct_count, fbank_count, ceplifter);
    // delete_2ddata((void**)mfcc_logf, time_step);
    return dct_mfcc_logf;
}

float** Mfcc::mfcc_from_string(byte* bytes) {
    float* datas = get_datas_from_string(bytes);
    float** ret = mfcc(datas);
    delete[] datas;
    return ret;
}

float** Mfcc::mfcc_from_file(const char* file_path) {
    float* datas = get_datas_from_file(file_path);
    float** ret = mfcc(datas);
    delete[] datas;
    return ret;
}