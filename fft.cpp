#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstring>

#define get_time(a) ((a)*1.0 / CLOCKS_PER_SEC)

const float PI=3.141592653589;

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

uint get_smaller_2_pow(uint val) {
    uint tmp = 1 << 31;
    int count = 32;
    while(!(tmp & val)) {
        tmp >>= 1;
        count -= 1;
    }
    tmp >>= 1;
    while(tmp) {
        if(tmp & val) {
            count += 1;
            break;
        }
        tmp >>= 1;
    }
    return (1 << (count-1));
}

float* dft(float* datas, uint data_len) {
    uint nfft = get_smaller_2_pow(data_len);
    float* tmp_datas = new float[nfft]();
    memcpy(tmp_datas, datas, data_len*sizeof(float));
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
        ret[i] = sqrt((data_r[i]*data_r[i] + data_i[i]*data_i[i]));
    }
    return ret;
}

complex* fft_unit(float* x, uint N) {
    if(N < 2) {
        return new complex(*(x), 0);
    }
    uint K = N / 2;
    float* xl = new float[K];
    float* xr = new float[K];
    for(int i=0; i<K; i++) {
        xl[i] = x[2*i];
        xr[i] = x[2*i + 1];
    }
    complex* cl = fft_unit(xl, K);
    complex* cr = fft_unit(xr, K);
    delete[] xl;
    delete[] xr;

    complex* c = new complex[N];
    for(int k=0; k<K; k++) {
        float dr = cos(2 * PI * (k*1.0 / N));
        float di = -sin(2 * PI * (k*1.0 / N));
        complex W = complex(dr, di);
        c[k] = cl[k] + W * cr[k];
        c[k + K] = cl[k] - W * cr[k];
    }
    delete[] cl;
    delete[] cr;
    return c;
}

float* fft(float* datas, uint data_len) {
    uint nfft = get_smaller_2_pow(data_len);
    float* tmp_datas = new float[nfft];
    memset(tmp_datas, 0, nfft * sizeof(float));
    memcpy(tmp_datas, datas, data_len*sizeof(float));
    float* ret = new float[nfft / 2 + 1];
    complex* d_fft = fft_unit(tmp_datas, nfft);

    for(int i=0; i<nfft / 2 + 1; i++) {
        float x = d_fft[i].x;
        float y = d_fft[i].y;
        ret[i] = sqrt(x*x + y*y);
    }
    delete[] d_fft;
    delete[] tmp_datas;
    return ret;
}

float* fft2(float* datas, uint data_len) {
    uint nfft = get_smaller_2_pow(data_len);
    float* tmp_datas = new float[nfft];
    memset(tmp_datas, 0, nfft * sizeof(float));
    memcpy(tmp_datas, datas, data_len*sizeof(float));
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
    delete[] rev;
    float* ret = new float[nfft / 2 + 1];
    for(uint N=2; N<=nfft; N<<=1) {
        uint K = N>>1;
        complex* W = new complex[K];
        for(int k=0; k<K; k++) {
            float dr = cos(2 * PI * (k*1.0 / N));
            float di = -sin(2 * PI * (k*1.0 / N));
            W[k] = complex(dr, di);
        }
        uint shift =0;
        for(int i=0; i<nfft/N; i++) {
            for(int k=0; k<K; k++) {
                uint idx1 = shift + k;
                uint idx2 = idx1 + K;
                complex v1 = d_fft[idx1] + W[k] * d_fft[idx2];
                complex v2 = d_fft[idx1] - W[k] * d_fft[idx2];
                d_fft[idx1] = v1;
                d_fft[idx2] = v2;
            }
            shift += N;
        }
        delete[] W;
    }
    for(int i=0; i<nfft / 2 + 1; i++) {
        float x = d_fft[i].x;
        float y = d_fft[i].y;
        ret[i] = sqrt(x*x + y*y);
    }
    delete[] d_fft;
    return ret;
}

float* fft3(float* datas, uint data_len)
{
    uint nfft = get_smaller_2_pow(data_len);
    float* tmp_datas = new float[nfft];
    memset(tmp_datas, 0, nfft * sizeof(float));
    memcpy(tmp_datas, datas, data_len*sizeof(float));
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

    for (int mid=1;mid<nfft;mid*=2)//mid是准备合并序列的长度的二分之一
    {
        complex temp(cos(PI/mid),sin(PI/mid));//单位根，pi的系数2已经约掉了
        for (int i=0;i<nfft;i+=mid*2)//mid*2是准备合并序列的长度，i是合并到了哪一位
		{
            complex omega(1,0);
            for (int j=0;j<mid;j++,omega=omega*temp)//只扫左半部分，得到右半部分的答案
            {
                complex x=d_fft[i+j],y=omega*d_fft[i+j+mid];
                d_fft[i+j]=x+y;
                d_fft[i+j+mid]=x-y;//这个就是蝴蝶变换什么的
            }
        }
    }
    for(int i=0; i<nfft / 2 + 1; i++) {
        float x = d_fft[i].x;
        float y = d_fft[i].y;
        ret[i] = sqrt(x*x + y*y);
    }
    delete[] d_fft;
    return ret;
}


int main() {
    clock_t c_start, c_end;
    ulong times = 99;
    uint data_len = 512;
    float* x = new float[data_len];
    for(uint i=0; i<data_len; i++) {
        x[i] = i;
    }
    c_start = clock();
    for(int time=0; time<times; time++) {
        float* dft_datas = fft3(x, data_len);
        delete[] dft_datas;
    }
    c_end = clock();
    printf("耗时%.5f\n", get_time(c_end - c_start));
    return 0;
}