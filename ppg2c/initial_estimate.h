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
void _fft(CPLX buf[], CPLX out[], int n, int step);//����Ҷ�任�ڶ���
int fft(CPLX buf[], int n);//����Ҷ�任��һ��
void shi2fu(CPLX buf[], double PPG_final[], int COUNT);//ʵ��ת����
void c_abs(CPLX f[], double out[], int n);//��������ȡģ
double findmax(double array[], int n);//���������ֵ
void initial_estimate(CPLX buf[], double PPG_final[], int COUNT, double srate, double Y_N_prev[]);//��ʼ���ʹ���

/*************************************** data_preprocess  *****************************************************/
void data_preprocess(CPLX buf[], double PPG1[], double srate);
void conjugate_complex(int n, CPLX in[], CPLX out[]);
void ifft(CPLX f[], int N);
void clean_up(CPLX buf[], double PPG1[], double srate, double PPG[]);

/*************************************** bpmFft3  *****************************************************/
void bpmFft3(double E[], double Y0, double srate, double N_prev, double BPM_N_prev[]);