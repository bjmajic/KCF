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
		Eigen::MatrixXf _alphaf_real; 
		Eigen::MatrixXf _alphaf_imag;

		Eigen::MatrixXf _prob_real;
		Eigen::MatrixXf _prob_imag;

		Eigen::MatrixXf _tmpl;
		Eigen::MatrixXf _num;
		Eigen::MatrixXf _den;
		//Eigen::MatrixXf _labCentroids;

		//Obtain sub-window from image, with replication-padding and extract features
		Eigen::MatrixXf getFeatures(const skMat& image, bool inithann, float scale_adjust = 1.0f);

		//Initialize Hanning window. Function called only in the first frame.
		void CreateHanningMats();

		//Calculate sub-pixel peak for one dimension
		float subPixelPeak(float left, float center, float right);

		// Evaluates a Gaussian kernel with bandwidth SIGMA for all relative shifts 
		// between input images X and Y, which must both be MxN. 
		// They must also be periodic (ie., pre-processed with a cosine window). 
		MatrixXf gaussianCorrelation(MatrixXf x1, MatrixXf x2);

		// Create Gaussian Peak. Function called only in the first frame. 
		void createGaussianPeak(int sizey, int sizex, MatrixXcf& resMatrix);
		void createGaussianPeak(int sizey, int sizex, MatrixXf& resMatrix_real, MatrixXf& resMatrix_imag);

		void FFT2D(const MatrixXf& in, MatrixXcf& out);
		void iFFT2D(const MatrixXcf& in, MatrixXf& out);

		void FFT2D_F(const MatrixXf& in, MatrixXcf& out);
		void iFFT2D_F(const MatrixXcf& in, MatrixXf& out);

		void FFT2D_F(const MatrixXf& in, MatrixXf& out_real, MatrixXf& out_imag);
		void iFFT2D_F(const MatrixXf& in_real, const MatrixXf& in_imag, MatrixXf& out);

		void FFT2D_F2(MatrixXf& in, MatrixXf& out_real, MatrixXf& out_imag);
		void iFFT2D_F2(MatrixXf& in_real, MatrixXf& in_imag, MatrixXf& out);
		

		//cv::Point2f detect(cv::Mat z, cv::Mat x, float &peak_value); // Detect object in the current frame. 
		void train(MatrixXf x, float train_interp_factor); // train tracker with a single image

		// Detect object in the current frame.
		skPoint2f detect(MatrixXf z, MatrixXf x, float &peak_value);

	private:
		CFFT fft_row;  // 对每一行进行fft
		CFFT fft_col;  // 对每一列进行fft

		bool _hogfeatures;
		int size_patch[3];
		float _scale;
		skSize<int> _tmpl_sz;
		MatrixXf hann;
		
		float *tmp_in_row_imag;
		float *tmp_out_real;
		float *tmp_out_imag;
	};

}


