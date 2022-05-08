#include<stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.14159265358979



typedef struct Complex
{
    double real;
    double image;
}CPLX;

void _fft(CPLX buf[], CPLX out[], int n, int step);//傅里叶变换第二步
int fft(CPLX buf[], int n);//傅里叶变换第一步
void shi2fu(CPLX buf[], double ppg[], int COUNT);//实数转复数
void c_abs(CPLX f[], double out[], int n);//复数数组取模
int findmax(double array[], int n);