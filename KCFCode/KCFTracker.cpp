#include "KCFTracker.h"

#include <stdexcept>
#undef eigen_assert
#define eigen_assert(x) \
  if (!(x)) { throw (std::runtime_error("input size must equle output size")); }

namespace SK
{
	KCFTracker::KCFTracker(bool hog /*= false*/, bool fixed_window /*= true*/, bool multiscale /*= false*/)
	{
		lambda = 0.0001f;
		padding = 2.5f;
		output_sigma_factor = 0.125f;

		tmp_in_row_imag = nullptr;
		tmp_out_real = nullptr;
		tmp_out_imag = nullptr;

		if (hog)
		{
		}
		else
		{
			interp_factor = 0.075f;
			sigma = 0.2f;
			cell_size = 1;
			_hogfeatures = false;
		}

		if (multiscale)
		{
		}
		else if (fixed_window)
		{
			template_size = 96;
			scale_step = 1;
		}
		else
		{
			template_size = 1;
			scale_step = 1;
		}
	}

	KCFTracker::~KCFTracker()
	{
		if (tmp_in_row_imag != nullptr)
		{
			delete[] tmp_in_row_imag;
			tmp_in_row_imag = nullptr;
		}
		if (tmp_out_real != nullptr)
		{
			delete[] tmp_out_real;
			tmp_out_real = nullptr;
		}
		if (tmp_out_imag != nullptr)
		{
			delete[] tmp_out_imag;
			tmp_out_imag = nullptr;
		}
	}

	void KCFTracker::Init(skRect& roi, skMat image)
	{
		//m_roi = roi;
		RectCopy<int, float>(roi, m_roi);
		assert(roi.width > 0 && roi.height > 0);
		_tmpl = getFeatures(image, 1);

		//saveTxt("init_tmpl", _tmpl);

		fft_row = CFFT(_tmpl.cols());
		fft_col = CFFT(_tmpl.rows());

		tmp_in_row_imag = new float[_tmpl.cols()];  // 输入图像每一行的虚部
		memset(tmp_in_row_imag, 0, _tmpl.cols()*sizeof(float));
		
		tmp_out_real = new float[_tmpl.rows() * _tmpl.cols()];
		memset(tmp_out_real, 0, _tmpl.rows() * _tmpl.cols() * sizeof(float));

		tmp_out_imag = new float[_tmpl.rows() * _tmpl.cols()];
		memset(tmp_out_imag, 0, _tmpl.rows() * _tmpl.cols() * sizeof(float));

		/*createGaussianPeak(size_patch[0], size_patch[1], _prob);
		_alphaf = MatrixXcf(size_patch[0], size_patch[1]);
		_alphaf.fill(complexf(0));*/

		createGaussianPeak(size_patch[0], size_patch[1], _prob_real, _prob_imag);
		_alphaf_real = MatrixXf(size_patch[0], size_patch[1]);
		_alphaf_real.fill(0);
		_alphaf_imag = MatrixXf(size_patch[0], size_patch[1]);
		_alphaf_imag.fill(0);

		train(_tmpl, 1.0); // train with initial frame
	}

	skRect KCFTracker::update(skMat image, float& peakValue)
	{
		if (m_roi.x + m_roi.width <= 0) m_roi.x = -m_roi.width + 1;
		if (m_roi.y + m_roi.height <= 0) m_roi.y = -m_roi.height + 1;
		if (m_roi.x >= image.cols() - 1) m_roi.x = image.cols() - 2;
		if (m_roi.y >= image.rows() - 1) m_roi.y = image.rows() - 2;

		float cx = m_roi.x + m_roi.width / 2.0f;
		float cy = m_roi.y + m_roi.height / 2.0f;

		//float peak_value;
		skPoint2f res = detect(_tmpl, getFeatures(image, 0, 1.0f), peakValue);
		//cout << "peak_value =" << peak_value << endl;

		if (scale_step != 1)
		{
		}

		// Adjust by cell size and _scale
		m_roi.x = cx - m_roi.width / 2.0f + ((float)res.x * cell_size * _scale);
		m_roi.y = cy - m_roi.height / 2.0f + ((float)res.y * cell_size * _scale);

		if (m_roi.x >= image.cols() - 1) m_roi.x = image.cols() - 1;
		if (m_roi.y >= image.rows() - 1) m_roi.y = image.rows() - 1;
		if (m_roi.x + m_roi.width <= 0) m_roi.x = -m_roi.width + 2;
		if (m_roi.y + m_roi.height <= 0) m_roi.y = -m_roi.height + 2;

		assert(m_roi.width >= 0 && m_roi.height >= 0);
		MatrixXf x = getFeatures(image, 0);
		train(x, interp_factor);

		skRect roi;
		RectCopy<float, int>(m_roi, roi);
		return roi;
		//return Rect(1, 2, 3, 4);
	}

	Eigen::MatrixXf KCFTracker::getFeatures(const skMat& image, bool inithann, float scale_adjust /*= 1.0f*/)
	{
		skRect extracted_roi;

		//真实roi区域的中心（未pading之前）
		float cx = m_roi.x + m_roi.width * 0.5f; 
		float cy = m_roi.y + m_roi.height * 0.5f;

		if (inithann)
		{
			//汉宁窗操作只有第一帧初始化的时候执行
			int padded_w = static_cast<int>(m_roi.width) * padding;
			int padded_h = static_cast<int>(m_roi.height) * padding;

			if (template_size > 1)
			{
				//设法逼近给定的template size
				_scale = (padded_w >= padded_h) ? padded_w / (float)template_size :
					padded_h / (float)template_size;

				//保存缩放后的尺寸
				_tmpl_sz.nWidth = static_cast<int>(padded_w / _scale);
				_tmpl_sz.nHeight = static_cast<int>(padded_h / _scale);
			}
			else
			{
				//没有给定template size， 使用 ROI size
				_tmpl_sz.nWidth = padded_w;
				_tmpl_sz.nHeight = padded_h;
				_scale = 1;
			}
			if (_hogfeatures)
			{
			}
			else
			{
				// 确保像素是偶数
				//_tmpl_sz.nWidth = (_tmpl_sz.nWidth / 2) * 2;
				//_tmpl_sz.nHeight = (_tmpl_sz.nHeight / 2) * 2;
				int i = 1;
				int j = 1;
				while (i < _tmpl_sz.nWidth)
				{
					i = i << 1;
				}
				int diff1 = _tmpl_sz.nWidth - (i >> 1);
				int diff2 = i - _tmpl_sz.nWidth;
				_tmpl_sz.nWidth = diff1 < diff2 ? i>>1 : i;

				while (j < _tmpl_sz.nHeight)
				{
					j = j << 1;
				}
				diff1 = _tmpl_sz.nHeight - (j >> 1);
				diff2 = j - _tmpl_sz.nHeight;
				_tmpl_sz.nHeight = diff1 < diff2 ? j>>1 : j;

				cout << "nwidth = " << _tmpl_sz.nWidth << endl;
				cout << "nheight = " << _tmpl_sz.nHeight << endl;
				//_tmpl_sz.nWidth = 256;
				//_tmpl_sz.nHeight = 256;
			}
		}

		//最终的结果还是 约等于 pading后的尺寸（因为有一个确保偶数的操作和取整操作）
		//为什么费这么大劲，反反复复，原因不明
		extracted_roi.width = scale_adjust * _scale * _tmpl_sz.nWidth;
		extracted_roi.height = scale_adjust * _scale * _tmpl_sz.nHeight;

		extracted_roi.x = cx - extracted_roi.width * 0.5;
		extracted_roi.y = cy - extracted_roi.height * 0.5;  // 这样中心点没有移动

		Eigen::MatrixXf FeaturesMap;
		skMat z = SK::subwindow(image, extracted_roi);
		if (z.cols() != _tmpl_sz.nWidth || z.rows() != _tmpl_sz.nHeight) {
			SK::Resize(z, _tmpl_sz);
		}

		MatrixXf t = z.cast<float>();
		//saveTxt("init_z.txt", t);

		if (_hogfeatures)
		{
		}
		else
		{
			FeaturesMap = getGrayImage(z);
			FeaturesMap.array() -= 0.5f;
			// size_patch[] 个人理解就是特征: 特征只有一个度量属性（或者说只有一个通道）, 宽为z.cols， 高为z.rows
			size_patch[0] = z.rows();
			size_patch[1] = z.cols();
			size_patch[2] = 1;

			//saveTxt("befor_hanning.txt", FeaturesMap);
		}
		if (inithann) {
			CreateHanningMats();
		}
		//FeaturesMap = hann.mul(FeaturesMap); opencv 逐元素乘法
		FeaturesMap = hann.cwiseProduct(FeaturesMap);
		return FeaturesMap;
	}

	void KCFTracker::CreateHanningMats()
	{
		MatrixXf hann1t(1, size_patch[1]); //（行， 列）
		hann1t.fill(0);
		MatrixXf hann2t(size_patch[0], 1); // (行， 列)
		hann2t.fill(0);
		for (int i = 0; i < hann1t.cols(); i++)
		{
			hann1t(0, i) = 0.5f * (1 - std::cos(2 * 3.141593f * i / (hann1t.cols() - 1)));
		}
		for (int i = 0; i < hann2t.rows(); i++)
		{
			hann2t(i, 0) = 0.5f * (1 - std::cos(2 * 3.141593f * i / (hann2t.rows() - 1)));
		}

		MatrixXf hann2d = hann2t * hann1t;
		if (_hogfeatures)
		{
		}
		else
		{
			// Gray features
			hann = hann2d;
		}
	}

	float KCFTracker::subPixelPeak(float left, float center, float right)
	{
		float divisor = 2 * center - right - left;

		if (divisor == 0)
			return 0;

		return 0.5f * (right - left) / divisor;
	}

	Eigen::MatrixXf KCFTracker::gaussianCorrelation(MatrixXf x1, MatrixXf x2)
	{
		//cv::Mat c = cv::Mat(cv::Size(size_patch[1], size_patch[0]), CV_32F, cv::Scalar(0));
		MatrixXf c(size_patch[0], size_patch[1]);
		c.fill(0);
		if (_hogfeatures)
		{
		}
		else
		{
			MatrixXf x1_out_real;
			MatrixXf x2_out_real;

			MatrixXf x1_out_imag;
			MatrixXf x2_out_imag;

			//FFT2D(x1, x1_out);
			//FFT2D(x2, x2_out);
			FFT2D_F(x1, x1_out_real, x1_out_imag);
			FFT2D_F(x2, x2_out_real, x2_out_imag);

			//saveTxt("gaussianCorrelation_x1_out.txt", x1_out);
			//saveTxt("gaussianCorrelation_x2_out.txt", x2_out);


			//iFFT2D(x1_out.cwiseProduct(x2_out.conjugate()), c);
			//iFFT2D_F(x1_out.cwiseProduct(x2_out.conjugate()), c);
			iFFT2D_F(x1_out_real.cwiseProduct(x2_out_real) + x1_out_imag.cwiseProduct(x2_out_imag),
				     x1_out_imag.cwiseProduct(x2_out_real) - x1_out_real.cwiseProduct(x2_out_imag), c);

			//saveTxt("gaussianCorrelation_c.txt", c);

			rearrange(c);
		}

		//saveTxt("gaussianCorrelation_c_.txt", c);

		MatrixXf d;
		d = x1.cwiseProduct(x1).sum() + x2.cwiseProduct(x2).sum() - (c * 2).array();
		d = d / (size_patch[0] * size_patch[1] * size_patch[2]);
		d = (d.array() > 0).select(d, 0);
		
		MatrixXf k;
		//cv::exp((-d / (sigma * sigma)), k);
		
		k = (-d / (sigma * sigma)).array().exp();
		return k;
	}

	void KCFTracker::createGaussianPeak(int sizey, int sizex, MatrixXcf& resMatrix)
	{
		MatrixXf res(sizey, sizex);

		int syh = (sizey) / 2;
		int sxh = (sizex) / 2;

		float output_sigma = std::sqrt((float)sizex * sizey) / padding * output_sigma_factor; // 带宽
		float mult = -0.5f / (output_sigma * output_sigma);

		for (int i = 0; i < sizey; i++)
			for (int j = 0; j < sizex; j++)
			{
				int ih = i - syh;
				int jh = j - sxh;
				res(i, j) = std::exp(mult * (float)(ih * ih + jh * jh));
			}

		//FFT2D(res, resMatrix);
		FFT2D_F(res, resMatrix);
	}

	void KCFTracker::createGaussianPeak(int sizey, int sizex, MatrixXf& resMatrix_real, MatrixXf& resMatrix_imag)
	{
		MatrixXf res(sizey, sizex);

		int syh = (sizey) / 2;
		int sxh = (sizex) / 2;

		float output_sigma = std::sqrt((float)sizex * sizey) / padding * output_sigma_factor; // 带宽
		float mult = -0.5f / (output_sigma * output_sigma);

		for (int i = 0; i < sizey; i++)
			for (int j = 0; j < sizex; j++)
			{
				int ih = i - syh;
				int jh = j - sxh;
				res(i, j) = std::exp(mult * (float)(ih * ih + jh * jh));
			}
		FFT2D_F(res, resMatrix_real, resMatrix_imag);
	}

	void KCFTracker::FFT2D(const MatrixXf& in, MatrixXcf& out)
	{
		out.resize(in.rows(), in.cols());

		for (int i = 0; i < in.rows(); i++)
		{
			//每一行进行一维fft
			vector<complexf> in_row;
			vector<complexf> out_row;
			for (int j = 0; j < in.cols(); j++)
			{
				in_row.push_back(in(i, j));
			}
			fft_row.fft1D(in_row, out_row);
			for (size_t j = 0; j < out_row.size(); j++)
			{
				out(i, j) = out_row[j];
			}
		}

		for (int i = 0; i < out.cols(); i++)
		{
			//每一列进行一维fft
			vector<complexf> in_col;
			vector<complexf> out_col;
			for (int j = 0; j < out.rows(); j++)
			{
				in_col.push_back(out(j, i));
			}
			fft_col.fft1D(in_col, out_col);
			for (size_t j = 0; j < out_col.size(); j++)
			{
				out(j, i) = out_col[j];
			}
		}
	}

	void KCFTracker::iFFT2D(const MatrixXcf& in, MatrixXf& out)
	{
		out.resize(in.rows(), in.cols());
		MatrixXcf outTmp(in.rows(), in.cols());

		for (int i = 0; i < in.rows(); i++)
		{
			//每一行进行一维fft
			vector<complexf> in_row;
			vector<complexf> out_row;
			for (int j = 0; j < in.cols(); j++)
			{
				in_row.push_back(in(i, j));
			}
			fft_row.ifft1D(in_row, out_row);
			for (size_t j = 0; j < out_row.size(); j++)
			{
				outTmp(i, j) = out_row[j];
			}
		}

		for (int i = 0; i < out.cols(); i++)
		{
			//每一列进行一维fft
			vector<complexf> in_col;
			vector<complexf> out_col;
			for (int j = 0; j < out.rows(); j++)
			{
				in_col.push_back(outTmp(j, i));
			}
			fft_col.ifft1D(in_col, out_col);
			for (size_t j = 0; j < out_col.size(); j++)
			{
				out(j, i) = out_col[j].real();
			}
		}
	}

	void KCFTracker::FFT2D_F(const MatrixXf& in, MatrixXcf& out)
	{
		out.resize(in.rows(), in.cols());
		vector< vector<float> > outTmp_r(in.rows(), vector<float>(in.cols()));
		vector< vector<float> > outTmp_i(in.rows(), vector<float>(in.cols()));

		int colSize = in.cols();
		vector<float> in_row_real(colSize);
		vector<float> in_row_imag(colSize);
		vector<float> out_row_real(colSize);
		vector<float> out_row_imag(colSize);

		for (int i = 0; i < in.rows(); i++)
		{
			//每一行进行一维fft
			
			for (int j = 0; j < colSize; j++)
			{
				in_row_real[j] = in(i, j);
				in_row_imag[j] = 0;
			}
			fft_row.fft1D(in_row_real, in_row_imag, out_row_real, out_row_imag);

			for (size_t j = 0; j < out_row_real.size(); j++)
			{
				//out(i, j) = complexf(out_row_real[j], out_row_imag[j]);
				outTmp_r[i][j] = out_row_real[j];
				outTmp_i[i][j] = out_row_imag[j];
			}
		}

		int rowSize = out.rows();
		vector<float> in_col_real(rowSize);
		vector<float> in_col_imag(rowSize);
		vector<float> out_col_real(rowSize);
		vector<float> out_col_imag(rowSize);

		for (int i = 0; i < out.cols(); i++)
		{
			//每一列进行一维fft

			for (int j = 0; j < rowSize; j++)
			{
				/*in_col_real[j] = out(j, i).real();
				in_col_imag[j] = out(j, i).imag();*/
				in_col_real[j] = outTmp_r[j][i];
				in_col_imag[j] = outTmp_i[j][i];
			}
			fft_col.fft1D(in_col_real, in_col_imag, out_col_real, out_col_imag);
			for (size_t j = 0; j < out_col_real.size(); j++)
			{
				out(j, i) = complexf(out_col_real[j], out_col_imag[j]);
			}
		}
	}

	void KCFTracker::FFT2D_F(const MatrixXf& in, MatrixXf& out_real, MatrixXf& out_imag)
	{		
		int rowSize = in.rows();
		int colSize = in.cols();

		out_real.resize(rowSize, colSize);
		out_imag.resize(rowSize, colSize);

		//vector< vector<float> > outTmp_r(rowSize, vector<float>(colSize));
		//vector< vector<float> > outTmp_i(rowSize, vector<float>(colSize));

		vector<float> in_row_real(colSize);
		vector<float> in_row_imag(colSize);
		vector<float> out_row_real(colSize);
		vector<float> out_row_imag(colSize);

		for (int i = 0; i < rowSize; i++)
		{
			//每一行进行一维fft

			for (int j = 0; j < colSize; j++)
			{
				in_row_real[j] = in(i, j);
				in_row_imag[j] = 0;
			}
			fft_row.fft1D(in_row_real, in_row_imag, out_row_real, out_row_imag);

			for (size_t j = 0; j < out_row_real.size(); j++)
			{
				//out(i, j) = complexf(out_row_real[j], out_row_imag[j]);
				out_real(i, j) = out_row_real[j];
				out_imag(i, j) = out_row_imag[j];
			}
		}

		
		vector<float> in_col_real(rowSize);
		vector<float> in_col_imag(rowSize);
		vector<float> out_col_real(rowSize);
		vector<float> out_col_imag(rowSize);

		for (int i = 0; i < colSize; i++)
		{
			//每一列进行一维fft

			for (int j = 0; j < rowSize; j++)
			{
				in_col_real[j] = out_real(j, i);
				in_col_imag[j] = out_imag(j, i);
			}
			fft_col.fft1D(in_col_real, in_col_imag, out_col_real, out_col_imag);
			for (size_t j = 0; j < out_col_real.size(); j++)
			{
				//out(j, i) = complexf(out_col_real[j], out_col_imag[j]);
				out_real(j, i) = out_col_real[j];
				out_imag(j, i) = out_col_imag[j];
			}
		}
	}

	void KCFTracker::iFFT2D_F(const MatrixXcf& in, MatrixXf& out)
	{
		out.resize(in.rows(), in.cols());
		//MatrixXcf outTmp(in.rows(), in.cols());
		vector< vector<float> > outTmp_r(in.rows(), vector<float>(in.cols()));
		vector< vector<float> > outTmp_i(in.rows(), vector<float>(in.cols()));

		int colSize = in.cols();
		vector<float> in_row_real(colSize);
		vector<float> in_row_imag(colSize);
		vector<float> out_row_real(colSize);
		vector<float> out_row_imag(colSize);

		for (int i = 0; i < in.rows(); i++)
		{
			//每一行进行一维fft
		
			for (int j = 0; j < colSize; j++)
			{
				in_row_real[j] = in(i, j).real();
				in_row_imag[j] = in(i, j).imag();
			}
			fft_row.ifft1D(in_row_real, in_row_imag, out_row_real, out_row_imag);
			for (size_t j = 0; j < out_row_real.size(); j++)
			{
				//outTmp(i, j) = complexf(out_row_real[j], out_row_imag[j]);
				outTmp_r[i][j] = out_row_real[j];
				outTmp_i[i][j] = out_row_imag[j];
			}
		}

		int rowSize = out.rows();
		vector<float> in_col_real(rowSize);
		vector<float> in_col_imag(rowSize);
		vector<float> out_col_real(rowSize);
		vector<float> out_col_imag(rowSize);

		for (int i = 0; i < out.cols(); i++)
		{
			//每一列进行一维fft
			
			for (int j = 0; j < rowSize; j++)
			{
				//in_col_real[j] = (outTmp(j, i).real());
				//in_col_imag[j] = (outTmp(j, i).imag());
				in_col_real[j] = outTmp_r[j][i];
				in_col_imag[j] = outTmp_i[j][i];;
			}
			fft_col.ifft1D(in_col_real, in_col_imag, out_col_real, out_col_imag);
			for (size_t j = 0; j < out_col_real.size(); j++)
			{
				out(j, i) = out_col_real[j];
			}
		}
	}

	void KCFTracker::iFFT2D_F(const MatrixXf& in_real, const MatrixXf& in_imag, MatrixXf& out)
	{
		eigen_assert(in_real.rows() == in_imag.rows() && in_real.cols() == in_imag.cols());
		out.resize(in_real.rows(), in_real.cols());

		int rowSize = out.rows();
		int colSize = out.cols();
		vector< vector<float> > outTmp_r(rowSize, vector<float>(colSize));
		vector< vector<float> > outTmp_i(rowSize, vector<float>(colSize));
		

		vector<float> in_row_real(colSize);
		vector<float> in_row_imag(colSize);
		vector<float> out_row_real(colSize);
		vector<float> out_row_imag(colSize);

		for (int i = 0; i < rowSize; i++)
		{
			//每一行进行一维fft

			for (int j = 0; j < colSize; j++)
			{
				in_row_real[j] = in_real(i, j);
				in_row_imag[j] = in_imag(i, j);
			}
			fft_row.ifft1D(in_row_real, in_row_imag, out_row_real, out_row_imag);
			for (size_t j = 0; j < out_row_real.size(); j++)
			{
				//outTmp(i, j) = complexf(out_row_real[j], out_row_imag[j]);
				outTmp_r[i][j] = out_row_real[j];
				outTmp_i[i][j] = out_row_imag[j];
			}
		}

		
		vector<float> in_col_real(rowSize);
		vector<float> in_col_imag(rowSize);
		vector<float> out_col_real(rowSize);
		vector<float> out_col_imag(rowSize);

		for (int i = 0; i < out.cols(); i++)
		{
			//每一列进行一维fft

			for (int j = 0; j < rowSize; j++)
			{
				//in_col_real[j] = (outTmp(j, i).real());
				//in_col_imag[j] = (outTmp(j, i).imag());
				in_col_real[j] = outTmp_r[j][i];
				in_col_imag[j] = outTmp_i[j][i];;
			}
			fft_col.ifft1D(in_col_real, in_col_imag, out_col_real, out_col_imag);
			for (size_t j = 0; j < out_col_real.size(); j++)
			{
				out(j, i) = out_col_real[j];
			}
		}
	}

	void KCFTracker::FFT2D_F2(MatrixXf& in, MatrixXf& out_real, MatrixXf& out_imag)
	{
		int rowSize = in.rows();
		int colSize = in.cols();

		out_real.resize(rowSize, colSize);
		out_imag.resize(rowSize, colSize);

		MatrixXf tmp_in_row_imag(1, colSize);
		tmp_in_row_imag.fill(0.0f);

		//vector<float> in_row_real(colSize);
		//vector<float> in_row_imag(colSize);
		//vector<float> out_row_real(colSize);
		//vector<float> out_row_imag(colSize);

		int bufOff = 0;

		for (int i = 0; i < rowSize; i++)
		{
			bufOff = i * colSize;
			fft_row.fft1D(in.data() + bufOff, tmp_in_row_imag.data(), out_real.data() + bufOff, out_imag.data() + bufOff);
		}

		// 转置
		out_real.transpose();
		out_imag.transpose();
		int rowNum = out_real.rows();
		int colNum = out_real.cols();

		for (int i = 0; i < rowNum; i++)
		{
			//每一列进行一维fft
			bufOff = i * colSize;
			fft_col.fft1D(in.data() + bufOff, tmp_in_row_imag.data(), out_real.data() + bufOff, out_imag.data() + bufOff);			
		}
	}

	void KCFTracker::train(MatrixXf x, float train_interp_factor)
	{
		MatrixXf k = gaussianCorrelation(x, x);

		//saveTxt("train_k.txt", k);

		MatrixXf  k_fft_real;
		MatrixXf  k_fft_imag;

		FFT2D_F(k, k_fft_real, k_fft_imag);

		//saveTxt("train_kfft.txt", k_fft);

		//k_fft.array() += lambda;
		//MatrixXcf alphaf = _prob.cwiseQuotient(k_fft);

		k_fft_imag.array() += lambda;
		k_fft_real.array() += lambda;

		MatrixXf tmp = k_fft_real.array().square() + k_fft_imag.array().square();

		MatrixXf alphaf_real = (_prob_real.cwiseProduct(k_fft_real) + _prob_imag.cwiseProduct(k_fft_imag)).array()/tmp.array();
		MatrixXf alphaf_imag = (_prob_imag.cwiseProduct(k_fft_real) - _prob_real.cwiseProduct(k_fft_imag)).array()/tmp.array();

		_tmpl = (1 - train_interp_factor) * _tmpl + (train_interp_factor)* x;
		//_alphaf = (1 - train_interp_factor) * _alphaf + (train_interp_factor)* alphaf;
		_alphaf_real = (1 - train_interp_factor) * _alphaf_real + (train_interp_factor)* alphaf_real;
		_alphaf_imag = (1 - train_interp_factor) * _alphaf_imag + (train_interp_factor)* alphaf_imag;
	}

	SK::skPoint2f KCFTracker::detect(MatrixXf z, MatrixXf x, float &peak_value)
	{
		MatrixXf k = gaussianCorrelation(x, z);
		
		MatrixXf k_fft_real;
		MatrixXf k_fft_imag;
		
		FFT2D_F(k, k_fft_real, k_fft_imag);

		MatrixXf res;
		//iFFT2D(_alphaf.cwiseProduct(k_fft), res);
		//iFFT2D_F(_alphaf.cwiseProduct(k_fft), res);
		iFFT2D_F(_alphaf_real.cwiseProduct(k_fft_real) - _alphaf_imag.cwiseProduct(k_fft_imag), 
			     _alphaf_imag.cwiseProduct(k_fft_real) + _alphaf_real.cwiseProduct(k_fft_imag) ,res);

		//saveTxt("detect_res.txt", res);

		skPoint<int> pi;
		//float pv;
		int r = 0;
		int c = 0;
		peak_value = res.maxCoeff(&r, &c);


		//subpixel peak estimation, coordinates will be non-integer
		skPoint2f p((float)c, (float)r);

		if (pi.x > 0 && pi.x < res.cols() - 1) {
			p.x += subPixelPeak(res(pi.y, pi.x - 1), peak_value, res(pi.y, pi.x + 1));
		}

		if (pi.y > 0 && pi.y < res.rows() - 1) {
			p.y += subPixelPeak(res(pi.y - 1, pi.x), peak_value, res(pi.y + 1, pi.x));
		}

		p.x -= (res.cols()) / 2;
		p.y -= (res.rows()) / 2;

		return p;
	}

} // namespace SK
