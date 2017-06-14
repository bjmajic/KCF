#pragma once
#include "Tracker.h"
#include "FFT.h"
namespace SK
{
	class KCFTracker : public Tracker
	{
	public:
		KCFTracker(bool hog = false, bool fixed_window = true, bool multiscale = false);
		~KCFTracker();

		virtual void Init(skRect& roi, skMat image);
		virtual skRect update(skMat image, float& peakValue);

		float interp_factor;
		float sigma;                 // gaussian kernel bandwidth
		float lambda;
		float padding;               // extra area surrounding the target
		float output_sigma_factor;   // bandwidth of gaussian target
		int   template_size;         // template size
		float scale_step;            // scale step for multi-scale estimation
		int   cell_size;

	protected:
		//Eigen::MatrixXcf _alphaf; //频域的值
		//Eigen::MatrixXcf _prob;  // 频域的值
		skMatrix _alphaf_real; 
		skMatrix _alphaf_imag;

		skMatrix _prob_real;
		skMatrix _prob_imag;

		skMatrix _tmpl;
		skMatrix _num;
		skMatrix _den;
		//Eigen::MatrixXf _labCentroids;

		//Obtain sub-window from image, with replication-padding and extract features
		skMatrix getFeatures(const skMat& image, bool inithann, float scale_adjust = 1.0f);

		//Initialize Hanning window. Function called only in the first frame.
		void CreateHanningMats();

		//Calculate sub-pixel peak for one dimension
		float subPixelPeak(float left, float center, float right);

		// Evaluates a Gaussian kernel with bandwidth SIGMA for all relative shifts 
		// between input images X and Y, which must both be MxN. 
		// They must also be periodic (ie., pre-processed with a cosine window). 
		skMatrix gaussianCorrelation(skMatrix x1, skMatrix x2);

		// Create Gaussian Peak. Function called only in the first frame. 
		void createGaussianPeak(int sizey, int sizex, MatrixXcf& resMatrix);
		void createGaussianPeak(int sizey, int sizex, skMatrix& resMatrix_real, skMatrix& resMatrix_imag);

		void FFT2D(const skMatrix& in, MatrixXcf& out);
		void iFFT2D(const MatrixXcf& in, skMatrix& out);

		void FFT2D_F(const skMatrix& in, MatrixXcf& out);
		void iFFT2D_F(const MatrixXcf& in, skMatrix& out);

		void FFT2D_F(const skMatrix& in, skMatrix& out_real, skMatrix& out_imag);
		void iFFT2D_F(const skMatrix& in_real, const skMatrix& in_imag, skMatrix& out);

		void FFT2D_F2(skMatrix& in, skMatrix& out_real, skMatrix& out_imag);
		void iFFT2D_F2(skMatrix& in_real, skMatrix& in_imag, skMatrix& out);
		

		//cv::Point2f detect(cv::Mat z, cv::Mat x, float &peak_value); // Detect object in the current frame. 
		void train(skMatrix x, float train_interp_factor); // train tracker with a single image

		// Detect object in the current frame.
		skPoint2f detect(skMatrix z, skMatrix x, float &peak_value);

	private:
		CFFT fft_row;  // 对每一行进行fft
		CFFT fft_col;  // 对每一列进行fft

		bool _hogfeatures;
		int size_patch[3];
		float _scale;
		skSize<int> _tmpl_sz;
		skMatrix hann;
		
		//float *tmp_in_row_imag;
		//skMatrix tmp_out_real;
		//skMatrix tmp_out_imag;
	};

}


