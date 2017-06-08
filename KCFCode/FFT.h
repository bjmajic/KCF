#pragma once
#include <complex>
#include <iostream>
#include <vector>
using namespace std;

typedef complex<float> complexf;
const float pi = 3.141593f;
class CFFT
{
public:
	CFFT();
	CFFT(int n); //������fft, n һ��Ϊ 2 ���ݴΣ� ������ǣ���Ҫ��0
	~CFFT();

public:
	bool bInitSuccess;
	void fft1D(vector<complexf>& in, vector<complexf>& out);
	void ifft1D(vector<complexf>& in, vector<complexf>& out);

	void fft1D(vector<float>& in_real, vector<float>& in_imag, vector<float>& out_real, vector<float>& out_imag);
	void ifft1D(vector<float>& in_real, vector<float>& in_imag, vector<float>& out_real, vector<float>& out_imag);

private:
	int m_nstage; //���м��ε������㣬��n=8��ʱ��m_nstage=3
	int  m_nZeroNum;  //��Ҫ��0�ĸ���
	int  m_nLen;      // Ԫ�ز�0��ĳ���
	vector<float> m_vecWm;  //ÿ���׶ε��η�֧�е�Wm����; ������ W = w * Wm (w = (1+j0))
	vector<int> m_vecOperationPairsNum; //ÿ���׶ε�������Ե�����
	vector<int> m_vecDistance; //ÿ���׶β�ͬ��ż��֮��ľ��룬���һ�׶εľ�����2���ڶ��׶εľ�����4, �����׶���8��...
	vector<int> m_nSoredIndex; // �洢�������к������, 001 --> 100 (4) 

	unsigned int bit_rev(unsigned int v, unsigned int maxv);
	void GetSortedIndex();
	void bit_reverse_copy(vector<complexf>& src, vector<complexf>& des);
	void bit_reverse_copy(vector<float>& src_real, vector<float>& src_imag, 
		vector<float>& des_real, vector<float>& des_imag);
};

