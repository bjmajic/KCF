#include "FftTool.h"
#include "Eigen/Dense"
#include "FFT.h"


void main()
{
	SK::complex com[8];
	for (int i = 1; i <= 8; i++)
	{
		//cout << i << "----" << bit_rev(i, 7) << endl;
		com[i - 1] = SK::complex(i, 0);
		if (i == 8 || i == 7)
		{
			com[i - 1] = SK::complex(0, 0);
		}
	}
	cout << "********自定义复数结构*******" << endl;
	for (int i = 0; i < 8; i++)
	{
		com[i].show();
		cout << endl;
	}
	cout << endl << endl;

	cout << "******ffttool 正变换*******" << endl;
	SK::complex fftArr[8];
	SK::fft1D(com, fftArr, 8);
	for (int i = 0; i < 8; i++)
	{
		fftArr[i].show();
		cout << endl;
	}
	cout << endl << endl;


	cout << "******ffttool 逆变换*******" << endl;
	SK::complex ifft[8];
	SK::ifft1D(fftArr, ifft, 8);
	for (int i = 0; i < 8; i++)
	{
		ifft[i].show();
		cout << endl;
	}
	cout << endl << endl;

	cout << "******fftclass 正变换*******" << endl;
	vector<complexf> vecf;
	for (int i = 0; i < 6; i++)
	{
		complexf a(i + 1, i);
		vecf.push_back(a);
		cout << a << endl;
	}
	CFFT fft(6);
	vector<complexf> res;
	if (fft.bInitSuccess)
	{
		fft.fft1D(vecf, res);
		cout << "res size = " << res.size() << endl;
	}
	for (int i = 0; i < res.size(); i++)
	{
		cout << res[i] << endl;
	}

	cout << "******fftclass 逆变换*******" << endl;
	vector<complexf> invres;
	fft.ifft1D(res, invres);
	for (int i = 0; i < invres.size(); i++)
	{
		cout << invres[i] << endl;
	}
}