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

		virtual void Init(Rect& roi, Mat image) override;
		virtual Rect update(Mat image) override;

		float interp_factor;
		float sigma;                 // gaussian kernel bandwidth
		float lambda;
		float padding;               // extra area surrounding the target
		float output_sigma_factor;   // bandwidth of gaussian target
		int   template_size;         // template size
		float scale_step;            // scale step for multi-scale estimation
		int   cell_size;

	protected:
		Eigen::MatrixXcf _alphaf; //频域的值
		Eigen::MatrixXcf _prob;  // 频域的值
		Eigen::MatrixXf _tmpl;
		Eigen::MatrixXf _num;
		Eigen::MatrixXf _den;
		//Eigen::MatrixXf _labCentroids;

		//Obtain sub-window from image, with replication-padding and extract features
		Eigen::MatrixXf getFeatures(const Mat& image, bool inithann, float scale_adjust = 1.0f);

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

		void FFT2D(const MatrixXf& in, MatrixXcf& out);
		void iFFT2D(const MatrixXcf& in, MatrixXf& out);
		

		//cv::Point2f detect(cv::Mat z, cv::Mat x, float &peak_value); // Detect object in the current frame. 
		void train(MatrixXf x, float train_interp_factor); // train tracker with a single image

		// Detect object in the current frame.
		Point2f detect(MatrixXf z, MatrixXf x, float &peak_value);

	private:
		CFFT fft_row;  // 对每一行进行fft
		CFFT fft_col;  // 对每一列进行fft

		bool _hogfeatures;
		int size_patch[3];
		float _scale;
		Size<int> _tmpl_sz;
		MatrixXf hann;
		
	};

}


