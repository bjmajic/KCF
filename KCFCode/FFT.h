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
	CFFT(int n); //做几点fft, n 一般为 2 的幂次， 如果不是，需要补0
	~CFFT();

public:
	bool bInitSuccess;
	void fft1D(vector<complexf>& in, vector<complexf>& out);
	void ifft1D(vector<complexf>& in, vector<complexf>& out);

	void fft1D(vector<float>& in_real, vector<float>& in_imag, vector<float>& out_real, vector<float>& out_imag);
	void ifft1D(vector<float>& in_real, vector<float>& in_imag, vector<float>& out_real, vector<float>& out_imag);

private:
	int m_nstage; //进行几次迭代运算，如n=8的时候，m_nstage=3
	int  m_nZeroNum;  //需要补0的个数
	int  m_nLen;      // 元素补0后的长度
	vector<float> m_vecWm;  //每个阶段蝶形分支中的Wm增量; 依次是 W = w * Wm (w = (1+j0))
	vector<int> m_vecOperationPairsNum; //每个阶段蝶形运算对的数量
	vector<int> m_vecDistance; //每个阶段不同奇偶组之间的距离，如第一阶段的距离是2，第二阶段的距离是4, 第三阶段是8，...
	vector<int> m_nSoredIndex; // 存储倒序排列后的索引, 001 --> 100 (4) 

	unsigned int bit_rev(unsigned int v, unsigned int maxv);
	void GetSortedIndex();
	void bit_reverse_copy(vector<complexf>& src, vector<complexf>& des);
	void bit_reverse_copy(vector<float>& src_real, vector<float>& src_imag, 
		vector<float>& des_real, vector<float>& des_imag);
};

