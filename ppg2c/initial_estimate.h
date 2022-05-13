#include<stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.14159265358979



typedef struct Complex
{
    double real;
    double image;
}CPLX;
/*************************************** initial_estimate  *****************************************************/
void _fft(CPLX buf[], CPLX out[], int n, int step);//傅里叶变换第二步
int fft(CPLX buf[], int n);//傅里叶变换第一步
void shi2fu(CPLX buf[], double PPG_final[], int COUNT);//实数转复数
void c_abs(CPLX f[], double out[], int n);//复数数组取模
double findmax(double array[], int n);//求数组最大值
void initial_estimate(CPLX buf[], double PPG_final[], int COUNT, double srate, double Y_N_prev[]);//初始心率估计

/*************************************** data_preprocess  *****************************************************/
void data_preprocess(CPLX buf[], double PPG1[], double srate);
void conjugate_complex(int n, CPLX in[], CPLX out[]);
void ifft(CPLX f[], int N);
void clean_up(CPLX buf[], double PPG1[], double srate, double PPG[]);

/*************************************** bpmFft3  *****************************************************/
void bpmFft3(double E[], double Y0, double srate, double N_prev, double BPM_N_prev[]);