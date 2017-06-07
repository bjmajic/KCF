#ifndef SK_FFT_TOOL_H
#define SK_FFT_TOOL_H
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

namespace SK
{
	const float pi  = 3.141593f;
	const float EPS = 0.000001f;
	class complex
	{
	public:
		complex()
		{
			re = 0.0f;
			im = 0.0f;
		}
		complex(float real, float imag)
		{
			re = real;
			im = imag;
		}
		complex operator + (complex& c)
		{
			return complex(re + c.re, im + c.im);
		}
		complex operator - (complex& c)
		{
			return complex(re - c.re, im - c.im);
		}
		complex operator * (complex& c)
		{
			return complex((re*c.re) - (im*c.im), (re*c.im + im*c.re));
		}
		complex operator / (complex& c)
		{
			float tmp = c.re * c.re + c.im * c.im;
			return complex((re*c.re + im*c.im) / tmp, (im*c.re - re*c.im) / tmp);
		}
		//ȡ����
		void conj()
		{
			im = -im;
		}

		void half()
		{
			re = re * 0.5f;
			im = im * 0.5f;
		}
		void scale(int n)
		{
			re = re / n;
			im = im / n;
		}
		complex c_pow(complex& c, int n)
		{
			if (n >= 1)
			{
				complex retValue = c;
				for (int i = 1; i < n; i++)
				{
					retValue = retValue * c;
				}
				return retValue;
			}
			else
			{
				return complex(0, 0);
			}			
		}

		//��ֵ
		void SetValue(float real, float imag)
		{
			re = abs(real) > EPS ? real : 0.0f;
			im = abs(imag) > EPS ? imag : 0.0f;
		}

		void show()
		{
			//std::cout << ((re >= 0) ? " " : "-") << re << ((im >= 0) ? ("+") : ("-") << ("j")) << im;
			cout.setf(ios::fixed);
			cout.precision(6);
			std::cout << re << "\t" << im;
		}

	protected:
	private:
		float re; // ʵ��
		float im; // �鲿
	}; // class complex


	//��ȡWn
	void GetW(complex W[], int len)
	{
		for (int i = 0; i < len; i++)
		{
			float fAngle = (pi * 2 * i) / len;
			W[i].SetValue(cos(fAngle), -sin(fAngle));
		}
	}

	//��ȡmWn, ������任
	void GetMW(complex W[], int len)
	{
		for (int i = 0; i < len; i++)
		{
			float fAngle = (pi * 2 * i) / len;
			W[i].SetValue(cos(fAngle), sin(fAngle));
		}
	}

	// ��������
	unsigned int bit_rev(unsigned int v, unsigned int maxv)
	{
		unsigned int t = log(maxv + 1) / log(2);
		unsigned int ret = 0;
		unsigned int s = 1;
		for (unsigned int i = 0; i < t; ++i)
		{
			unsigned int r = v&(s << i); // ��ȡv��bitλ�ĵ�iλֵ��0 or 1
			ret |= r << (t - 1 - i) >> (i); // ������ t-1-i λ�� ������iλ
		}
		return ret;
	}
	void bit_reverse_copy(complex src[], complex des[], int len)
	{
		for (int i = 0; i < len; i++)
		{
			des[bit_rev(i, len - 1)] = src[i]; //		
		}
	}

	void fft1D(complex in[], complex out[], int len, bool bInverse = false)
	{
		bit_reverse_copy(in, out, len);
		int n = len;
		int nStage = log(n) / log(2); // ���м����׶ε����㣬��8��fft Ϊ������Ҫ3�λ�2����
		for (int s = 1; s <= nStage; s++)
		{
			int m = pow(2, s); //1�����ڼ���ÿ���׶ε�Wn����, 2������ÿ���׶β�ͬ��֮��ľ��룬���һ�׶εľ�����2���ڶ��׶εľ�����4
			complex wm(cos(2 * pi / m), -sin(2 * pi / m));
			for (int k = 0; k < n; k += m)
			{
				complex w(1, 0);
				for (int j = 0; j < m / 2; ++j) //��ǰ�׶ο�����ɶ��ٶԵ�������
				{
					complex vv = out[k + j + m / 2];
					complex t = w * vv;
					complex u = out[k + j];
					out[k + j] = u + t;
					out[k + j + m / 2] = u - t;
					w = w * wm;
				}
			}
		}

		if (true == bInverse)
		{
			for (int i = 0; i < n; ++i)
			{
				out[i].conj();
				out[i].scale(n);
			}
		}
	}
	void ifft1D(complex in[], complex out[], int len)
	{
		bit_reverse_copy(in, out, len);
		int n = len;
		int nStage = log(n) / log(2); // ���м����׶ε����㣬��8��fft Ϊ������Ҫ3�λ�2����
		for (int s = 1; s <= nStage; s++)
		{
			int m = pow(2, s); //1�����ڼ���ÿ���׶ε�Wn����, 2������ÿ���׶β�ͬ��֮��ľ��룬���һ�׶εľ�����2���ڶ��׶εľ�����4
			complex wm(cos(2 * pi / m), sin(2 * pi / m));
			for (int k = 0; k < n; k += m)
			{
				complex w(1, 0);
				for (int j = 0; j < m / 2; ++j) //��ǰ�׶ο�����ɶ��ٶԵ�������
				{
					complex vv = out[k + j + m / 2];
					complex t = w * vv;
					complex u = out[k + j];
					out[k + j] = u + t;
					out[k + j + m / 2] = u - t;
					w = w * wm;
				}
			}
		}
		for (int i = 0; i < n; ++i)
		{
			out[i].scale(n);
		}
	}



} // namespace SK
#endif // !SK_FFT_TOOL_H
