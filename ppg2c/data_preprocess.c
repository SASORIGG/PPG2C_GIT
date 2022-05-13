//#include "data_preprocess.h"
#include "initial_estimate.h"
#include <math.h>
//¸´ÊýµÄ½»»» 
void conjugate_complex(int n, CPLX in[], CPLX out[])
{
	int i = 0;
	for (i = 0; i < n; i++)
	{
		out[i].image = -in[i].image;
		out[i].real = in[i].real;
	}
}

//¸µÀïÒ¶Äæ±ä»»
void ifft(CPLX f[], int N )
{
	int i = 0;
	conjugate_complex(N, f, f);
	fft(f, N);
	conjugate_complex(N, f, f);
	for (i = 0; i < N; i++)
	{
		f[i].image = (f[i].image) / N;
		f[i].real = (f[i].real) / N;
	}
}

void clean_up(CPLX buf[], double PPG1[], double srate, double PPG[])
{
	double np[1024];
	shi2fu(buf, PPG1, 1024);
	fft(buf, 1024);
	//²âÊÔ
	//for (int i = 0; i < 1024; i++)
	//{
	//	printf("%d: ", i + 1);
	//	if (buf[i].image >= 0.0)
	//	{
	//		printf("%.4f + %.4fi\n", buf[i].real, buf[i].image);
	//	}
	//	else
	//	{
	//		printf("%.4f - %.4fi\n", buf[i].real, fabs(buf[i].image));
	//	}
	//}
	for (int i = 0; i < 1024; i++)
	{
		np[i] = i * 2 * PI / 1024;
		//printf("%d: ", i + 1);
		//printf("%f\n", np[i]);
	}
	for (int j = 0; j < 1024; j++)
	{
		if (np[j] < 2 * PI * 0.4 / srate | np[j] > 2 * PI * (1 - 0.4 / srate))
		{
			buf[j].real = 0;
			buf[j].image = 0;
		}
		else if (np[j] > 2 * PI * 3.5 / srate & np[j] < 2 * PI * (1 - 3.5 / srate))
		{
			buf[j].real = 0;
			buf[j].image = 0;
		}
	}
	//²âÊÔ
	//for (int i = 0; i < 1024; i++)
	//{
	//	printf("%d: ", i + 1);
	//	if (buf[i].image >= 0.0)
	//	{
	//		printf("%.4f + %.4fi\n", buf[i].real, buf[i].image);
	//	}
	//	else
	//	{
	//		printf("%.4f - %.4fi\n", buf[i].real, fabs(buf[i].image));
	//	}
	//}
	ifft(buf, 1024);
	//²âÊÔ
	for (int i = 0; i < 1000; i++)
	{
		PPG[i] = buf[i].real;
	}

}

void data_preprocess(CPLX buf[], double PPG1[], double srate)
{
	double PPG[1000];
	double PPG_sum = 0;
	double PPG_norm = 0;
	clean_up(buf, PPG1, srate, PPG);
	//²âÊÔ
	for (int i = 0; i < 1000; i++)
	{
		printf("%d: ", i + 1);
		printf("%f\n", PPG[i]);
		PPG_sum = PPG_sum + PPG[i] * PPG[i];
	}
	PPG_norm = sqrt(PPG_sum);
	printf("%f", PPG_norm);

}