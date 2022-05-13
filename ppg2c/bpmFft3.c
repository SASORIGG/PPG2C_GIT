#include "initial_estimate.h"
#include <math.h>

double findmax_value(double EXabs[], double m)
{
	double max = EXabs[0];
	int i = 0;
	double posNum = 0;
	for (i = 0; i < m; i++)
	{
		if (max < EXabs[i])
		{
			max = EXabs[i];
			posNum = i + 1;
		}
	}

	return max;
}

/*
*	EXabs[]：				在里面找峰值
*	EXabs_peak[]：			峰值大小写在里面
*	number：				峰值个数
*	EXabs_peak_location[]：	峰值频率对应的位置
*/
int findpeaks(double EXabs[], double EXabs_peak[], int number, double EXabs_peak_location[])
{
	int n = 0;
	int loc = 0;
	for (int i = 1; i < number; i++)
	{
		if (EXabs[i] > EXabs[i - 1] && EXabs[i] > EXabs[i + 1])
		{
			EXabs_peak[n++] = EXabs[i];
			EXabs_peak_location[loc++] = i+1;
			//printf("%f: ", EXabs_peak_location[loc - 1]);	//峰值频率所在的位置
			//printf("%f\n", EXabs_peak[n - 1]);				//峰值频率大小
		}
	}
	//for (int i = 0; i < n; i++)
	//{
	//	printf("%d:", i + 1);
	//	printf("%f\n", EXabs_peak[i]);
	//}
	return n;
}
/*
*	E[]：			去伪影后的PPG信号
*	Y0：			上一个时间窗口计算得到的心率
*	srate：			采样频率
*	N_prev：		上一个时间窗口心率估计点的位置
*	BPM_N_prev[2]：	第一个数存当前计算的心率值，第二个数存该心率对应的位置
*/
void bpmFft3(double E[], double Y0, double srate, double N_prev, double BPM_N_prev[])
{
	CPLX E_c[4096];	//复数形式的去伪影后的PPG信号
	CPLX EX[4096];	//存放 PPG 信号 FFT 后前 m 个数
	double EXabs[132];	//存放对 EX 取模后的结果
	double EXabs_peak[66];	//存放 EXabs 内的峰值大小
	double EXabs_peak_location[66];	//存放 EXabs 内的峰值位置
	double searchScale_peak[10];	//存放搜索范围内的峰值
	double searchScale_peak_location[10];	//存放搜索范围内峰值对应的位置
	double M, m;
	double I2;	//初始估计峰值位置
	double Y;	//初始估计心率
	double Y2;
	double DEL = 10;	//搜索区域Δs=10
	double T0 = 10;		//阈值
	double I;
	int EXabs_peak_num, searchScale_peak_num;	//范围内找到的峰值的个数
	double search_scale[21];
	int start, end;
	int n = 1;
	int n2 = 0;
	double maxValue;
	int p = 0;	//大于0.25倍的峰值频率个数

	//将去伪影后的PPG信号 E 转成复数形式 E_c，4096个采样点
	shi2fu(E_c, E, 4096);
	////测试
	//for (int i = 0; i < 1000; i++)
	//{
	//	printf("%d: ", i + 1);
	//	if (E_c[i].image >= 0.0)
	//	{
	//		printf("%.4f + %.4fi\n", E_c[i].real, E_c[i].image);
	//	}
	//	else
	//	{
	//		printf("%.4f - %.4fi\n", E_c[i].real, fabs(E_c[i].image));
	//	}
	//}

	//对 E_c 做4096点的 FFT，结果仍然存放在 E_c 里
	fft(E_c, 4096);
	////测试
	//for (int i = 0; i < 4096; i++)
	//{
	//	printf("%d: ", i + 1);
	//	if (E_c[i].image >= 0.0)
	//	{
	//		printf("%.4f + %.4fi\n", E_c[i].real, E_c[i].image);
	//	}
	//	else
	//	{
	//		printf("%.4f - %.4fi\n", E_c[i].real, fabs(E_c[i].image));
	//	}
	//}


	M = srate * 60 / 4096;
	m = floor(242 / M);
	//FFT 后取前 m(132) 个最关键的放到 EX 里
	for (int i = 0; i < m; i++)
	{
		EX[i].real = E_c[i].real;
		EX[i].image = E_c[i].image;
	}
	////测试
	//for (int i = 0; i < m; i++)
	//{
	//	printf("%d: ", i + 1);
	//	if (E_c[i].image >= 0.0)
	//	{
	//		printf("%.4f + %.4fi\n", EX[i].real, EX[i].image);
	//	}
	//	else
	//	{
	//		printf("%.4f - %.4fi\n", EX[i].real, fabs(EX[i].image));
	//	}
	//}

	//对 EX 前 m 个数取模，结果放在 EXabs 里
	c_abs(EX, EXabs, m);//复数数组取模，数组 EXabs 为取模后的结果
	// 测试
	//for (int i = 0; i < m; i++)
	//{
	//    printf("%d: ", i + 1);
	//    printf("%.4f\n", EXabs[i]);
	//}

	//在 m 个 EXabs 中找到最大值，得到初始估计心率位置 I2 ，和计算出初始心率 Y
	I2 = findmax(EXabs, m);
	BPM_N_prev[0] = (I2 - 1) * M;
	// 测试
	//printf("%f, %f", I2, Y);
	I = I2;
	//在 EXabs 找峰值
	//EXabs_peak_num 存放峰值个数，EXabs_peak[] 存放峰值大小，EXabs_peak_location[] 存放峰值位置
	EXabs_peak_num = findpeaks(EXabs, EXabs_peak, 131, EXabs_peak_location);
	//printf("\n%d\n", EXabs_peak_num);
	maxValue = findmax_value(EXabs, m);
	//printf("%\n%f\n", maxValue);
	for (int i = 0; i < EXabs_peak_num; i++)
	{
		if (EXabs_peak[i] > 0.25 * maxValue)
		{
			p++;
		}
	}
	if (p <= 1)
	{
		BPM_N_prev[0] = (I - 1) * M;
		BPM_N_prev[1] = I; 
		return;
	}
	//如果上一个时间窗口得到的心率大于零
	if (Y0 > 0)
	{
		I = N_prev;
		//printf("%f\n", I);
		start = (int)(I - 10);	//搜索范围的左边界
		end = (int)(I + 10);	//搜索范围的右边界
		//printf("%d	%d\n", start, end);

		for (int i = start - 1; i < end; i++)
		{
			search_scale[n2++] = EXabs[i];	//21个搜索范围
			//printf("%d:", n++);
			//printf("%f\n", EXabs[i]);
		}
		//searchScale_peak_num 存放峰值个数，searchScale_peak[] 存放峰值大小
		searchScale_peak_num = findpeaks(search_scale, searchScale_peak, 21, searchScale_peak_location);
		//将最大频率对应的位置改成绝对位置（不是搜索范围内的位置）
		for (int i = 0; i < searchScale_peak_num; i++)
		{
			searchScale_peak_location[i] += I - 10 - 1;	//加上左边界
		}//现在的 searchScale_peak_location 里存的是搜索范围内峰值的绝对位置（132里的）（对应 matlab 代码里的 ai）
		////测试
		//for (int i = 0; i < searchScale_peak_num; i++)
		//{
		//	printf("%f\n", searchScale_peak_location[i]);
		//}

		//去掉前面 searchScale_peak_location 里包含边界的峰值（即包含 start 和 end 的峰值）
		if (searchScale_peak_location[0] == (double)start)
		{
			for (int i = 0; i < searchScale_peak_num; i++)
			{
				searchScale_peak_location[i] = searchScale_peak_location[i + 1];
			}
			searchScale_peak_num = searchScale_peak_num - 1;
		}
		if (searchScale_peak_location[searchScale_peak_num - 1] == (double)end)
			searchScale_peak_num = searchScale_peak_num - 1;
		//测试
		//printf("\n");
		//for (int i = 0; i < searchScale_peak_num; i++)
		//{
		//	printf("%f：", searchScale_peak_location[i]);
		//	printf("%f\n", searchScale_peak[i]);
		//}
		//printf("\n%d\n", searchScale_peak_num);

		I = findmax(searchScale_peak, searchScale_peak_num); //第 I 个数
		I = searchScale_peak_location[(int)I - 1];	//对应数组 I 要 减1
		//printf("\n%f\n", I);
		Y2 = (I - 1) * M;
		if (abs(BPM_N_prev[0] - Y0) < T0)
		{
			BPM_N_prev[1] = I2;
		}
		else
		{
			BPM_N_prev[1] = I;
			BPM_N_prev[0] = Y2;
		}
	}
}