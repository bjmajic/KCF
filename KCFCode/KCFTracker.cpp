#include "KCFTracker.h"
#include <cassert>
namespace SK
{
	KCFTracker::KCFTracker(bool hog /*= false*/, bool fixed_window /*= true*/, bool multiscale /*= false*/)
	{
		lambda = 0.0001f;
		padding = 2.5f;
		output_sigma_factor = 0.125f;
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

		createGaussianPeak(size_patch[0], size_patch[1], _prob);
		_alphaf = MatrixXcf(size_patch[0], size_patch[1]);
		_alphaf.fill(complexf(0));

		train(_tmpl, 1.0); // train with initial frame
	}

	skRect KCFTracker::update(skMat image)
	{
		if (m_roi.x + m_roi.width <= 0) m_roi.x = -m_roi.width + 1;
		if (m_roi.y + m_roi.height <= 0) m_roi.y = -m_roi.height + 1;
		if (m_roi.x >= image.cols() - 1) m_roi.x = image.cols() - 2;
		if (m_roi.y >= image.rows() - 1) m_roi.y = image.rows() - 2;

		float cx = m_roi.x + m_roi.width / 2.0f;
		float cy = m_roi.y + m_roi.height / 2.0f;

		float peak_value;
		skPoint2f res = detect(_tmpl, getFeatures(image, 0, 1.0f), peak_value);

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
				_tmpl_sz.nWidth = (_tmpl_sz.nWidth / 2) * 2;
				_tmpl_sz.nHeight = (_tmpl_sz.nHeight / 2) * 2;

				_tmpl_sz.nWidth = 256;
				_tmpl_sz.nHeight = 256;
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
			//灰度特征
			//MatrixXf c(x1.rows(), x1.cols());

			MatrixXcf x1_out;
			MatrixXcf x2_out;
			FFT2D(x1, x1_out);
			FFT2D(x2, x2_out);

			//saveTxt("gaussianCorrelation_x1_out.txt", x1_out);
			//saveTxt("gaussianCorrelation_x2_out.txt", x2_out);


			iFFT2D(x1_out.cwiseProduct(x2_out.conjugate()), c);

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
		//MatrixXcf res_fft;
		//FFT2D(res, res_fft);
		//return res_fft;
		//resMatrix.resize(res.rows(), res.cols());
		FFT2D(res, resMatrix);
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

		for (int i = 0; i < in.rows(); i++)
		{
			//每一行进行一维fft
			vector<float> in_row_real;
			vector<float> in_row_imag;
			vector<float> out_row_real;
			vector<float> out_row_imag;
			for (int j = 0; j < in.cols(); j++)
			{
				in_row_real.push_back(in(i, j));
				in_row_imag.push_back(0);
			}
			fft_row.fft1D(in_row_real, in_row_imag, out_row_real, out_row_imag);

			for (size_t j = 0; j < out_row_real.size(); j++)
			{
				out(i, j) = complexf(out_row_real[j], out_row_imag[j]);
			}
		}

		for (int i = 0; i < out.cols(); i++)
		{
			//每一列进行一维fft
			vector<float> in_col_real;
			vector<float> in_col_imag;
			vector<float> out_col_real;
			vector<float> out_col_imag;

			for (int j = 0; j < out.rows(); j++)
			{
				in_col_real.push_back(out(j, i).real());
				in_col_imag.push_back(out(j, i).imag());
			}
			fft_col.fft1D(in_col_real, in_col_imag, out_col_real, out_col_imag);
			for (size_t j = 0; j < out_col_real.size(); j++)
			{
				out(j, i) = complexf(out_col_real[j], out_col_imag[j]);
			}
		}
	}

	void KCFTracker::iFFT2D_F(const MatrixXcf& in, MatrixXf& out)
	{
		out.resize(in.rows(), in.cols());
		MatrixXcf outTmp(in.rows(), in.cols());

		for (int i = 0; i < in.rows(); i++)
		{
			//每一行进行一维fft
			vector<float> in_row_real;
			vector<float> in_row_imag;
			vector<float> out_row_real;
			vector<float> out_row_imag;

			for (int j = 0; j < in.cols(); j++)
			{
				in_row_real.push_back(in(i, j).real());
				in_row_imag.push_back(in(i, j).imag());
			}
			fft_row.ifft1D(in_row_real, in_row_imag, out_row_real, out_row_imag);
			for (size_t j = 0; j < out_row_real.size(); j++)
			{
				outTmp(i, j) = complexf(out_row_real[j], out_row_imag[j]);
			}
		}

		for (int i = 0; i < out.cols(); i++)
		{
			//每一列进行一维fft
			vector<float> in_col_real;
			vector<float> in_col_imag;
			vector<float> out_col_real;
			vector<float> out_col_imag;
			for (int j = 0; j < out.rows(); j++)
			{
				in_col_real.push_back(outTmp(j, i).real());
				in_col_imag.push_back(outTmp(j, i).imag());
			}
			fft_col.ifft1D(in_col_real, in_col_imag, out_col_real, out_col_imag);
			for (size_t j = 0; j < out_col_real.size(); j++)
			{
				out(j, i) = out_col_real[j];
			}
		}
	}

	void KCFTracker::train(MatrixXf x, float train_interp_factor)
	{
		MatrixXf k = gaussianCorrelation(x, x);

		//saveTxt("train_k.txt", k);

		MatrixXcf  k_fft;
		FFT2D(k, k_fft);

		//saveTxt("train_kfft.txt", k_fft);

		k_fft.array() += lambda;
		MatrixXcf alphaf = _prob.cwiseQuotient(k_fft);

		_tmpl = (1 - train_interp_factor) * _tmpl + (train_interp_factor)* x;
		_alphaf = (1 - train_interp_factor) * _alphaf + (train_interp_factor)* alphaf;

	}

	SK::skPoint2f KCFTracker::detect(MatrixXf z, MatrixXf x, float &peak_value)
	{
		MatrixXf k = gaussianCorrelation(x, z);
		//cv::Mat res = (real(fftd(complexMultiplication(_alphaf, fftd(k)), true)));
		MatrixXcf k_fft;
		FFT2D(k, k_fft);

		MatrixXf res;
		//iFFT2D(_alphaf.cwiseProduct(k_fft.conjugate()), res);
		iFFT2D(_alphaf.cwiseProduct(k_fft), res);

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
