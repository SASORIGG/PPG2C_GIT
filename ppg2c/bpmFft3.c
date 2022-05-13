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
*	EXabs[]��				�������ҷ�ֵ
*	EXabs_peak[]��			��ֵ��Сд������
*	number��				��ֵ����
*	EXabs_peak_location[]��	��ֵƵ�ʶ�Ӧ��λ��
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
			//printf("%f: ", EXabs_peak_location[loc - 1]);	//��ֵƵ�����ڵ�λ��
			//printf("%f\n", EXabs_peak[n - 1]);				//��ֵƵ�ʴ�С
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
*	E[]��			ȥαӰ���PPG�ź�
*	Y0��			��һ��ʱ�䴰�ڼ���õ�������
*	srate��			����Ƶ��
*	N_prev��		��һ��ʱ�䴰�����ʹ��Ƶ��λ��
*	BPM_N_prev[2]��	��һ�����浱ǰ���������ֵ���ڶ�����������ʶ�Ӧ��λ��
*/
void bpmFft3(double E[], double Y0, double srate, double N_prev, double BPM_N_prev[])
{
	CPLX E_c[4096];	//������ʽ��ȥαӰ���PPG�ź�
	CPLX EX[4096];	//��� PPG �ź� FFT ��ǰ m ����
	double EXabs[132];	//��Ŷ� EX ȡģ��Ľ��
	double EXabs_peak[66];	//��� EXabs �ڵķ�ֵ��С
	double EXabs_peak_location[66];	//��� EXabs �ڵķ�ֵλ��
	double searchScale_peak[10];	//���������Χ�ڵķ�ֵ
	double searchScale_peak_location[10];	//���������Χ�ڷ�ֵ��Ӧ��λ��
	double M, m;
	double I2;	//��ʼ���Ʒ�ֵλ��
	double Y;	//��ʼ��������
	double Y2;
	double DEL = 10;	//��������s=10
	double T0 = 10;		//��ֵ
	double I;
	int EXabs_peak_num, searchScale_peak_num;	//��Χ���ҵ��ķ�ֵ�ĸ���
	double search_scale[21];
	int start, end;
	int n = 1;
	int n2 = 0;
	double maxValue;
	int p = 0;	//����0.25���ķ�ֵƵ�ʸ���

	//��ȥαӰ���PPG�ź� E ת�ɸ�����ʽ E_c��4096��������
	shi2fu(E_c, E, 4096);
	////����
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

	//�� E_c ��4096��� FFT�������Ȼ����� E_c ��
	fft(E_c, 4096);
	////����
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
	//FFT ��ȡǰ m(132) ����ؼ��ķŵ� EX ��
	for (int i = 0; i < m; i++)
	{
		EX[i].real = E_c[i].real;
		EX[i].image = E_c[i].image;
	}
	////����
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

	//�� EX ǰ m ����ȡģ��������� EXabs ��
	c_abs(EX, EXabs, m);//��������ȡģ������ EXabs Ϊȡģ��Ľ��
	// ����
	//for (int i = 0; i < m; i++)
	//{
	//    printf("%d: ", i + 1);
	//    printf("%.4f\n", EXabs[i]);
	//}

	//�� m �� EXabs ���ҵ����ֵ���õ���ʼ��������λ�� I2 ���ͼ������ʼ���� Y
	I2 = findmax(EXabs, m);
	BPM_N_prev[0] = (I2 - 1) * M;
	// ����
	//printf("%f, %f", I2, Y);
	I = I2;
	//�� EXabs �ҷ�ֵ
	//EXabs_peak_num ��ŷ�ֵ������EXabs_peak[] ��ŷ�ֵ��С��EXabs_peak_location[] ��ŷ�ֵλ��
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
	//�����һ��ʱ�䴰�ڵõ������ʴ�����
	if (Y0 > 0)
	{
		I = N_prev;
		//printf("%f\n", I);
		start = (int)(I - 10);	//������Χ����߽�
		end = (int)(I + 10);	//������Χ���ұ߽�
		//printf("%d	%d\n", start, end);

		for (int i = start - 1; i < end; i++)
		{
			search_scale[n2++] = EXabs[i];	//21��������Χ
			//printf("%d:", n++);
			//printf("%f\n", EXabs[i]);
		}
		//searchScale_peak_num ��ŷ�ֵ������searchScale_peak[] ��ŷ�ֵ��С
		searchScale_peak_num = findpeaks(search_scale, searchScale_peak, 21, searchScale_peak_location);
		//�����Ƶ�ʶ�Ӧ��λ�øĳɾ���λ�ã�����������Χ�ڵ�λ�ã�
		for (int i = 0; i < searchScale_peak_num; i++)
		{
			searchScale_peak_location[i] += I - 10 - 1;	//������߽�
		}//���ڵ� searchScale_peak_location ������������Χ�ڷ�ֵ�ľ���λ�ã�132��ģ�����Ӧ matlab ������� ai��
		////����
		//for (int i = 0; i < searchScale_peak_num; i++)
		//{
		//	printf("%f\n", searchScale_peak_location[i]);
		//}

		//ȥ��ǰ�� searchScale_peak_location ������߽�ķ�ֵ�������� start �� end �ķ�ֵ��
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
		//����
		//printf("\n");
		//for (int i = 0; i < searchScale_peak_num; i++)
		//{
		//	printf("%f��", searchScale_peak_location[i]);
		//	printf("%f\n", searchScale_peak[i]);
		//}
		//printf("\n%d\n", searchScale_peak_num);

		I = findmax(searchScale_peak, searchScale_peak_num); //�� I ����
		I = searchScale_peak_location[(int)I - 1];	//��Ӧ���� I Ҫ ��1
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