#include<stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.14159265358979



typedef struct Complex
{
    double real;
    double image;
}CPLX;

void _fft(CPLX buf[], CPLX out[], int n, int step);//����Ҷ�任�ڶ���
int fft(CPLX buf[], int n);//����Ҷ�任��һ��
void shi2fu(CPLX buf[], double ppg[], int COUNT);//ʵ��ת����
void c_abs(CPLX f[], double out[], int n);//��������ȡģ
int findmax(double array[], int n);