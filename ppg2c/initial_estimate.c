#include"initial_estimate.h"

void _fft(CPLX buf[], CPLX out[], int n, int step)
{
    int i;
    CPLX t;
    double a, b, c, d, theta;

    if (step < n)
    {
        _fft(out, buf, n, step * 2);
        _fft(out + step, buf + step, n, step * 2);

        for (i = 0; i < n; i += 2 * step)
        {
            theta = -PI * i / n;
            a = cos(theta);
            b = sin(theta);
            c = out[i + step].real;
            d = out[i + step].image;
            t.real = a * c - b * d;
            t.image = b * c + a * d;
            buf[i / 2].real = out[i].real + t.real;
            buf[i / 2].image = out[i].image + t.image;
            buf[(i + n) / 2].real = out[i].real - t.real;
            buf[(i + n) / 2].image = out[i].image - t.image;
        }
    }
}

int fft(CPLX buf[], int n)
{
    int i;
    CPLX* out = (CPLX*)malloc(n * sizeof(CPLX));
    if (out == NULL)
    {
        return 0;
    }

    for (i = 0; i < n; i++)
    {
        out[i].real = buf[i].real;
        out[i].image = buf[i].image;
    }

    _fft(buf, out, n, 1);
    free(out);
    out = NULL;
    return 1;
}

void shi2fu(CPLX buf[], double PPG_final[], int COUNT)
{
    for (int i = 0; i < COUNT; i++)
    {
        buf[i].real = PPG_final[i];
        buf[i].image = 0;
    }
}

//求所有复数的模
void c_abs(CPLX f[], double out[], int n)
{
    int i = 0;
    double t;
    for (i = 0; i < n; i++)
    {
        t = f[i].real * f[i].real + f[i].image * f[i].image;
        out[i] = sqrt(t);
    }
}

double findmax(double array[], int n)
{
    double max = array[0];
    int i = 0;
    double posNum = 0;
    for (i = 0; i < n; i++)
    {
        if (max < array[i])
        {
            max = array[i];
            posNum = i + 1;
        }
    }

    return posNum;
}

void initial_estimate(CPLX buf[], double PPG_final[], int COUNT, double srate, double Y_N_prev[])
{
    double M, m;
    CPLX EX[4096];
    double EXout[4096];

    shi2fu(buf, PPG_final, COUNT);
    //对 buf 完成 COUNT点 的FFT，变换后仍为 buf
    fft(buf, COUNT);
    M = srate * 60 / COUNT;
    m = floor(242 / M);
    for (int i = 0; i < m; i++)
    {
        EX[i].real = buf[i].real;
        EX[i].image = buf[i].image;
    }
    c_abs(EX, EXout, m);//复数数组取模，数组 EXout 为取模后的结果
    //for (int i = 0; i < m; i++)
    //{
    //    printf("%d: ", i + 1);
    //    printf("%.4f\n", EXout[i]);
    //}
    Y_N_prev[1] = findmax(EXout, m);
    Y_N_prev[0] = (Y_N_prev[1] - 1) * M;

}