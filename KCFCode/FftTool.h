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
		//取共轭
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

		//设值
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
		float re; // 实部
		float im; // 虚部
	}; // class complex


	//获取Wn
	void GetW(complex W[], int len)
	{
		for (int i = 0; i < len; i++)
		{
			float fAngle = (pi * 2 * i) / len;
			W[i].SetValue(cos(fAngle), -sin(fAngle));
		}
	}

	//获取mWn, 用于逆变换
	void GetMW(complex W[], int len)
	{
		for (int i = 0; i < len; i++)
		{
			float fAngle = (pi * 2 * i) / len;
			W[i].SetValue(cos(fAngle), sin(fAngle));
		}
	}

	// 逆序重排
	unsigned int bit_rev(unsigned int v, unsigned int maxv)
	{
		unsigned int t = log(maxv + 1) / log(2);
		unsigned int ret = 0;
		unsigned int s = 1;
		for (unsigned int i = 0; i < t; ++i)
		{
			unsigned int r = v&(s << i); // 获取v的bit位的第i位值：0 or 1
			ret |= r << (t - 1 - i) >> (i); // 先左移 t-1-i 位， 再右移i位
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
		int nStage = log(n) / log(2); // 进行几个阶段的运算，以8点fft 为例，需要3次基2运算
		for (int s = 1; s <= nStage; s++)
		{
			int m = pow(2, s); //1、用于计算每个阶段的Wn增量, 2、控制每个阶段不同组之间的距离，如第一阶段的距离是2，第二阶段的距离是4
			complex wm(cos(2 * pi / m), -sin(2 * pi / m));
			for (int k = 0; k < n; k += m)
			{
				complex w(1, 0);
				for (int j = 0; j < m / 2; ++j) //当前阶段可以组成多少对蝶形运算
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
		int nStage = log(n) / log(2); // 进行几个阶段的运算，以8点fft 为例，需要3次基2运算
		for (int s = 1; s <= nStage; s++)
		{
			int m = pow(2, s); //1、用于计算每个阶段的Wn增量, 2、控制每个阶段不同组之间的距离，如第一阶段的距离是2，第二阶段的距离是4
			complex wm(cos(2 * pi / m), sin(2 * pi / m));
			for (int k = 0; k < n; k += m)
			{
				complex w(1, 0);
				for (int j = 0; j < m / 2; ++j) //当前阶段可以组成多少对蝶形运算
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
