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

void shi2fu(CPLX buf[], double ppg[], int COUNT)
{
    for (int i = 0; i < COUNT; i++)
    {
        buf[i].real = ppg[i];
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

int findmax(double array[], int n)
{
    double max = array[0];
    int i = 0;
    int posNum = 0;
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