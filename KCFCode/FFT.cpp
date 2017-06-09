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
	}
	out.resize(m_nLen);

	bit_reverse_copy(in, out);
	for (int s = 0; s < m_nstage; s++)
	{
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
				//complexf vv(out[index1].real(), out[index1].imag());
				float vv[2] = { out[index1].real(), out[index1].imag() };
				//complexf t = w * vv;
				//complexf t((w[0]*vv.real() - w[1]*vv.imag()), (w[0]*vv.imag() + w[1]*vv.real()));
				float t[2] = { w[0] * vv[0] - w[1] * vv[1], w[0] * vv[1] + w[1] * vv[0] };
				//complexf u = out[k + j];
				//complexf u(out[index2].real(), out[index2].imag());
				float out_index2_r = out[index2].real();
				float out_index2_i = out[index2].imag();
				float u[2] = { out_index2_r, out_index2_i };
				//out[k + j] = out[k + j] + t;
				out[index2].real(out_index2_r + t[0]);
				out[index2].imag(out_index2_i + t[1]);
				//out[j + pairsNumAddK] = u - t;
				out[index1].real(u[0] - t[0]);
				out[index1].imag(u[1] - t[1]);
				//w *= cmf;
				float new_r = w[0] * cmf[0] - w[1] * cmf[1];
				float new_i = w[0] * cmf[1] + w[1] * cmf[0];
				w[0] = new_r;
				w[1] = new_i;
			}
		}
	}
}

void CFFT::fft1D(vector<float>& in_real, vector<float>& in_imag, vector<float>& out_real, vector<float>& out_imag)
{
	int inLen = in_real.size();
	if (inLen != m_nLen)
	{
		in_real.resize(m_nLen); //会自动赋值
		in_imag.resize(m_nLen);
	}
	out_real.resize(m_nLen);
	out_imag.resize(m_nLen);

	bit_reverse_copy(in_real, in_imag, out_real, out_imag);

	for (int s = 0; s < m_nstage; s++)
	{
		int distance = m_vecDistance[s];
		int pairsNum = m_vecOperationPairsNum[s];
		//float cmf[2] = { m_vecWm[s * 2], m_vecWm[s * 2 + 1] };
		float cmf_0 = m_vecWm[s * 2];
		float cmf_1 = m_vecWm[s * 2 + 1];

		for (int k = 0; k < m_nLen; k += distance)
		{
			int pairsNumAddK = k + pairsNum;
			//complexf w(1, 0);
			//float w[2] = { 1, 0 };
			float w_0 = 1.0f;
			float w_1 = 0.0f;

			for (int j = 0; j < pairsNum; ++j) //当前阶段可以组成多少对蝶形运算
			{
				int index1 = j + pairsNumAddK;
				int index2 = k + j;
				//float vv[2] = { out[index1].real(), out[index1].imag() };
				float vv_0 = out_real[index1];
				float vv_1 = out_imag[index1];
				
				//float t[2] = { w[0] * vv[0] - w[1] * vv[1], w[0] * vv[1] + w[1] * vv[0] };
				float t_0 = w_0 * vv_0 - w_1 * vv_1;
				float t_1 = w_0 * vv_1 + w_1 * vv_0;
				
				//float out_index2_r = out[index2].real();
				//float out_index2_i = out[index2].imag();
				//float u[2] = { out_index2_r, out_index2_i };
				float u_0 = out_real[index2];
				float u_1 = out_imag[index2];
				//float u_0 = out_index2_r;
				//float u_1 = out_index2_i;
				
				//out[index2].real(out_index2_r + t[0]);
				//out[index2].imag(out_index2_i + t[1]);
				out_real[index2] = (u_0 + t_0);
				out_imag[index2] = (u_1 + t_1);
				
				//out[index1].real(u[0] - t[0]);
				//out[index1].imag(u[1] - t[1]);
				out_real[index1] = (u_0 - t_0);
				out_imag[index1] = (u_1 - t_1);

				//w *= cmf;
				float new_r = w_0 * cmf_0 - w_1 * cmf_1;
				float new_i = w_0 * cmf_1 + w_1 * cmf_0;
				w_0 = new_r;
				w_1 = new_i;
			}
		}
	}
}

void CFFT::ifft1D(vector<complexf>& in, vector<complexf>& out)
{
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
		//complexf cmf = m_vecWm[s];
		float cmf[2] = { m_vecWm[s * 2], m_vecWm[s * 2 + 1] };

		for (int k = 0; k < m_nLen; k += distance)
		{
			int pairsNumAddK = k + pairsNum;
			//complexf w(1, 0);
			float w[2] = { 1, 0 };
			float w_r = 1.0f;
			float w_i = 0.0f;

			for (int j = 0; j < pairsNum; ++j) //当前阶段可以组成多少对蝶形运算
			{
				int index1 = j + pairsNumAddK;
				int index2 = k + j;
				//complexf vv = out[j + pairsNumAddK];
				//complexf vv(out[index1].real(), out[index1].imag());
				float vv[2] = { out[index1].real(), out[index1].imag() };
				//complexf t = w * vv;
				//complexf t((w[0]*vv.real() - w[1]*vv.imag()), (w[0]*vv.imag() + w[1]*vv.real()));
				float t[2] = { w[0] * vv[0] - w[1] * vv[1], w[0] * vv[1] + w[1] * vv[0] };
				//complexf u = out[k + j];
				//complexf u(out[index2].real(), out[index2].imag());
				float out_index2_r = out[index2].real();
				float out_index2_i = out[index2].imag();
				float u[2] = { out_index2_r, out_index2_i };
				//out[k + j] = out[k + j] + t;
				out[index2].real(out_index2_r + t[0]);
				out[index2].imag(out_index2_i + t[1]);
				//out[j + pairsNumAddK] = u - t;
				out[index1].real(u[0] - t[0]);
				out[index1].imag(u[1] - t[1]);
				//w *= cmf;
				float new_r = w[0] * cmf[0] + w[1] * cmf[1];
				float new_i = w[1] * cmf[0] - w[0] * cmf[1];
				w[0] = new_r;
				w[1] = new_i;
			}
		}
	}
	for (int i = 0; i < m_nLen; ++i)
	{
		out[i] /= m_nLen;
	}
}

void CFFT::ifft1D(vector<float>& in_real, vector<float>& in_imag, vector<float>& out_real, vector<float>& out_imag)
{
	int inLen = in_real.size();
	if (inLen != m_nLen)
	{
		in_real.resize(m_nLen); //会自动赋值
		in_imag.resize(m_nLen);
	}
	out_real.resize(m_nLen);
	out_imag.resize(m_nLen);

	bit_reverse_copy(in_real, in_imag, out_real, out_imag);

	for (int s = 0; s < m_nstage; s++)
	{
		int distance = m_vecDistance[s];
		int pairsNum = m_vecOperationPairsNum[s];
		//float cmf[2] = { m_vecWm[s * 2], m_vecWm[s * 2 + 1] };
		float cmf_0 = m_vecWm[s * 2];
		float cmf_1 = m_vecWm[s * 2 + 1];

		for (int k = 0; k < m_nLen; k += distance)
		{
			int pairsNumAddK = k + pairsNum;
			//complexf w(1, 0);
			//float w[2] = { 1, 0 };
			float w_0 = 1.0f;
			float w_1 = 0.0f;

			for (int j = 0; j < pairsNum; ++j) //当前阶段可以组成多少对蝶形运算
			{
				int index1 = j + pairsNumAddK;
				int index2 = k + j;
				//float vv[2] = { out[index1].real(), out[index1].imag() };
				float vv_0 = out_real[index1];
				float vv_1 = out_imag[index1];

				//float t[2] = { w[0] * vv[0] - w[1] * vv[1], w[0] * vv[1] + w[1] * vv[0] };
				float t_0 = w_0 * vv_0 - w_1 * vv_1;
				float t_1 = w_0 * vv_1 + w_1 * vv_0;

				//float out_index2_r = out[index2].real();
				//float out_index2_i = out[index2].imag();
				//float u[2] = { out_index2_r, out_index2_i };
				float u_0 = out_real[index2];
				float u_1 = out_imag[index2];
				//float u_0 = out_index2_r;
				//float u_1 = out_index2_i;

				//out[index2].real(out_index2_r + t[0]);
				//out[index2].imag(out_index2_i + t[1]);
				out_real[index2] = (u_0 + t_0);
				out_imag[index2] = (u_1 + t_1);

				//out[index1].real(u[0] - t[0]);
				//out[index1].imag(u[1] - t[1]);
				out_real[index1] = (u_0 - t_0);
				out_imag[index1] = (u_1 - t_1);

				//w *= cmf;
				float new_r = w_0 * cmf_0 + w_1 * cmf_1;
				float new_i = w_1 * cmf_0 - w_0 * cmf_1;
				w_0 = new_r;
				w_1 = new_i;
			}
		}
	}
	for (int i = 0; i < m_nLen; ++i)
	{
		out_real[i] /= m_nLen;
		out_imag[i] /= m_nLen;
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

void CFFT::bit_reverse_copy(vector<float>& src_real, vector<float>& src_imag, vector<float>& des_real, vector<float>& des_imag)
{
	for (int i = 0; i < m_nLen; i++)
	{
		int j = m_nSoredIndex[i];
		des_real[j] = src_real[i]; //
		des_imag[j] = src_imag[i]; //
	}
}
