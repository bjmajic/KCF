#include "FFT.h"
#include <cmath>

CFFT::CFFT()
{
}


CFFT::CFFT(int n)
{
	try
	{
		if (n <= 0)
		{			
			throw "n must be greater than 0";
		}

		// 一般而言，n必须是2的幂次，这里为了兼容性，对非2的幂次进行补0操作
		unsigned int i = 2;
		m_nstage = 1;
		m_nZeroNum = 0;
		m_nLen = n;

		//int nBaseM = 2;
		m_vecDistance.push_back(i);
		m_vecOperationPairsNum.push_back(i / 2);
		//m_vecWm.push_back(complexf(cos(2 * pi / i), -sin(2 * pi / i)));  //complex wm(cos(2 * pi / m), -sin(2 * pi / m)); 此时m = 2
		m_vecWm.push_back(cos(2 * pi / i));
		m_vecWm.push_back(-sin(2 * pi / i));

		while (i < n)
		{
			i = i << 1;
			m_nstage++;

			//nBaseM = nBaseM << 1;
			m_vecDistance.push_back(i);
			m_vecOperationPairsNum.push_back(i / 2);
			//m_vecWm.push_back(complexf(cos(2 * pi / i), -sin(2 * pi / i)));
			m_vecWm.push_back(cos(2 * pi / i));
			m_vecWm.push_back(-sin(2 * pi / i));
		}

		m_nZeroNum = i - n;
		m_nLen = i;
		GetSortedIndex();

		if (m_vecWm.size() != m_vecOperationPairsNum.size())
		{
			throw "m_vecOperationPairsNum  m_vecWm size must be equal";
		}

		bInitSuccess = true;
	}
	catch (...)
	{
		bInitSuccess = false;
	}
}

CFFT::~CFFT()
{
}

void CFFT::fft1D(vector<complexf>& in, vector<complexf>& out)
{
	int inLen = in.size();
	if (inLen != m_nLen)
	{
		in.resize(m_nLen); //会自动赋值
		/*for (size_t i = inLen; i < m_nLen; i++)
		{
			in.push_back(complexf(0));
		}*/
	}
	out.resize(m_nLen);

	bit_reverse_copy(in, out);
	for (int s = 0; s < m_nstage; s++)
	{
		
		//int pm_vecOperationPairsNum = m_vecOperationPairsNum[s];
		int distance = m_vecDistance[s];
		int pairsNum = m_vecOperationPairsNum[s];
		//complexf cmf = m_vecWm[s];
		float cmf[2] = { m_vecWm[s * 2], m_vecWm[s * 2 + 1] };

		for (int k = 0; k < m_nLen; k += distance)
		{
			int pairsNumAddK = k + pairsNum;
			//complexf w(1, 0);
			float w[2] = { 1, 0 };

			for (int j = 0; j < pairsNum; ++j) //当前阶段可以组成多少对蝶形运算
			{
				int index1 = j + pairsNumAddK;
				int index2 = k + j;
				//complexf vv = out[j + pairsNumAddK];
				complexf vv(out[index1].real(), out[index1].imag());
				//complexf t = w * vv;
				complexf t((w.real()*vv.real() - w.imag()*vv.imag()), (w.real()*vv.imag() + w.imag()*vv.real()));
				//complexf u = out[k + j];
				complexf u(out[index2].real(), out[index2].imag());
				//out[k + j] = out[k + j] + t;
				out[index2].real(out[index2].real() + t.real());
				out[index2].imag(out[index2].imag() + t.imag());
				//out[j + pairsNumAddK] = u - t;
				out[index1].real(u.real() - t.real());
				out[index1].imag(u.imag() - t.imag());
				//w *= cmf;
				float new_r = w.real() * cmf.real() - w.imag() * cmf.imag();
				float new_i = w.real() * cmf.imag() + w.imag() * cmf.real();
				w.real(new_r);
				w.imag(new_i);
			}
		}
	}
}

void CFFT::ifft1D(vector<complexf>& in, vector<complexf>& out)
{
	//int inLen = in.size();
	//if (inLen != m_nLen)
	//{
	//	in.resize(m_nLen); //会自动赋值
	//}
	//out.resize(m_nLen);

	//bit_reverse_copy(in, out);
	//for (int s = 0; s < m_nstage; s++)
	//{
	//	for (int k = 0; k < m_nLen; k += m_vecDistance[s])
	//	{
	//		complexf w(1, 0);
	//		for (int j = 0; j < m_vecOperationPairsNum[s]; ++j) //当前阶段可以组成多少对蝶形运算
	//		{
	//			complexf vv = out[k + j + m_vecOperationPairsNum[s]];
	//			complexf t = w * vv;
	//			complexf u = out[k + j];
	//			out[k + j] = u + t;
	//			out[k + j + m_vecOperationPairsNum[s]] = u - t;
	//			w = w * conj(m_vecWm[s]);
	//		}
	//	}
	//}

	int inLen = in.size();
	if (inLen != m_nLen)
	{
		in.resize(m_nLen); //会自动赋值
	}
	out.resize(m_nLen);

	bit_reverse_copy(in, out);
	for (int s = 0; s < m_nstage; s++)
	{

		//int pm_vecOperationPairsNum = m_vecOperationPairsNum[s];
		int distance = m_vecDistance[s];
		int pairsNum = m_vecOperationPairsNum[s];
		complexf cmf = m_vecWm[s];

		for (int k = 0; k < m_nLen; k += distance)
		{
			int pairsNumAddK = k + pairsNum;
			complexf w(1, 0);
			for (int j = 0; j < pairsNum; ++j) //当前阶段可以组成多少对蝶形运算
			{
				int index1 = j + pairsNumAddK;
				int index2 = k + j;
				//complexf vv = out[j + pairsNumAddK];
				complexf vv(out[index1].real(), out[index1].imag());
				//complexf t = w * vv;
				complexf t((w.real()*vv.real() - w.imag()*vv.imag()), (w.real()*vv.imag() + w.imag()*vv.real()));
				//complexf u = out[k + j];
				complexf u(out[index2].real(), out[index2].imag());
				//out[k + j] = out[k + j] + t;
				out[index2].real(out[index2].real() + t.real());
				out[index2].imag(out[index2].imag() + t.imag());
				//out[j + pairsNumAddK] = u - t;
				out[index1].real(u.real() - t.real());
				out[index1].imag(u.imag() - t.imag());
				//w *= cmf;
				float new_r = w.real() * cmf.real() + w.imag() * cmf.imag();
				float new_i = w.imag() * cmf.real() - w.real() * cmf.imag();
				w.real(new_r);
				w.imag(new_i);
			}
		}
	}
	for (int i = 0; i < m_nLen; ++i)
	{
		out[i] /= m_nLen;
	}
}

unsigned int CFFT::bit_rev(unsigned int v, unsigned int maxv)
{
	unsigned int t = m_nstage/*log(maxv + 1) / log(2)*/;
	unsigned int ret = 0;
	unsigned int s = 1;
	for (unsigned int i = 0; i < t; ++i)
	{
		unsigned int r = v&(s << i); // 获取v的bit位的第i位值：0 or 1
		ret |= r << (t - 1 - i) >> (i); // 先左移 t-1-i 位， 再右移i位
	}
	return ret;
}

void CFFT::GetSortedIndex()
{
	int len = m_nLen - 1;
	for (int i = 0; i < m_nLen; i++)
	{
		m_nSoredIndex.push_back(bit_rev(i, len));
	}
}

void CFFT::bit_reverse_copy(vector<complexf>& src, vector<complexf>& des)
{
	for (int i = 0; i < m_nLen; i++)
	{
		des[m_nSoredIndex[i]] = src[i]; //		
	}
}
