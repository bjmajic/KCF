#include "FftTool.h"
#include "Eigen/Dense"
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include "FFT.h"
#include "KCFTracker.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include "TrackerManager.h"

using namespace std;
//using namespace cv;

bool HOG = false;
bool FIXEDWINDOW = false;
bool MULTISCALE = false;
bool SILENT = true;
//bool LAB = false;

string TestImagePath = "D:\\KCF\\KCF\\FastCNN\\KCF_Test_Image\\";
string ImageType = ".JPG";

int xMin = 1140;
int yMin = 200;
int width = 80;
int height = 100;

string GetImagePath(int num)
{
	ostringstream stream;
	stream << num;
	string StringNum = stream.str();
	return TestImagePath + StringNum + ImageType;
}

void main()
{
	TrackerManager tm;
	cv::Mat  tmpMat;
	bool bFirst = true;
	bool stop1 = false;
	cv::VideoCapture capture("D:/KCF/1.mp4");//要读取的视频文
	if (!capture.isOpened())
	{
		cout << "con not open video 1.mp4" << endl;
		return;
	}
	double rate = capture.get(CV_CAP_PROP_FPS);

	int delay = 1000 / rate;
	while (!stop1)
	{
		if (!capture.read(tmpMat))
		{
			break;
		}
		//cv::Mat grayMat;
		//cv::cvtColor(tmpMat, grayMat, cv::COLOR_BGR2GRAY);

		//Eigen::Map<skMat> frame2(grayMat.ptr<uchar>(), grayMat.rows, grayMat.cols);
		if (bFirst)
		{
			tm.Init(tmpMat.data, tmpMat.cols, tmpMat.rows, 5);
			//tracker2.Init(SK::skRect(290, 184, 140, 130), frame2);
			//cv::imwrite("D:\\KCF\\1.jpg", framecv2);
			bFirst = false;

			tm.CreateTracker(290, 184, 140, 130);
			tm.CreateTracker(188,84,54,56);
			tm.CreateTracker(506,56,66,66);
		}

		int* pRect;
		float  fValue;
		int n;
		//SK::skRect result2 = tracker2.update(frame2, ret);
		pRect = tm.UpdateTracker(tmpMat, n);
		for (int i = 0; i < n; i++)
		{
			cv::Point P1, P2;
			P1.x = pRect[i * 4 + 0];
			P1.y = pRect[i * 4 + 1];
			P2.x = P1.x + pRect[i * 4 + 2];
			P2.y = P1.y + pRect[i * 4 + 3];
			cv::rectangle(tmpMat, P1, P2, cv::Scalar(128, 0, 128, 0.3), 2);
		}
		//cout << ret << endl;
	
		cv::imshow("Result", tmpMat);
		cv::waitKey(1);
	}

	return;

	//SK::KCFTracker tracker2(HOG, FIXEDWINDOW, MULTISCALE);
	//cv::Mat  framecv2;
	//bool bFirstFrame = true;
	//bool stop = false;
	//cv::VideoCapture capture2("D:/KCF/1.mp4");//要读取的视频文
	//if (!capture.isOpened())
	//{
	//	cout << "con not open video 1.mp4" << endl;
	//	return;
	//}
	//double rate2 = capture.get(CV_CAP_PROP_FPS);
	//
	//int delay2 = 1000 / rate2;
	//while (!stop)
	//{
	//	if (!capture.read(framecv2))
	//	{
	//		break;
	//	}
	//	cv::Mat grayMat;
	//	cv::cvtColor(framecv2, grayMat, cv::COLOR_BGR2GRAY);
	//	
	//	Eigen::Map<skMat> frame2(grayMat.ptr<uchar>(), grayMat.rows, grayMat.cols);
	//	if (bFirstFrame)
	//	{
	//		tracker2.Init(SK::skRect(290, 184, 140, 130), frame2);
	//		//cv::imwrite("D:\\KCF\\1.jpg", framecv2);
	//		bFirstFrame = false;
	//	}

	//	float ret;
	//	SK::skRect result2 = tracker2.update(frame2, ret);
	//	//cout << ret << endl;
	//	cv::Point P1, P2;
	//	P1.x = result2.x;
	//	P1.y = result2.y;
	//	P2.x = result2.x + result2.width;
	//	P2.y = result2.y + result2.height;
	//	cv::rectangle(framecv2, P1, P2, cv::Scalar(128, 0, 128, 0.3), 2);
	//	cv::imshow("Result", framecv2);
	//	cv::waitKey(1);
	//}
	//return;

	// just for single image
	SK::KCFTracker tracker(HOG, FIXEDWINDOW, MULTISCALE);
	cv::Mat  framecv;
	//cv::Rect result;
	//frame = cv::imread(GetImagePath(0));
	string filename = GetImagePath(0);
	int h = 0;
	int w = 0;

	unsigned char* p = nullptr;
	
	LoadJPG(filename.c_str(), p, w, h);

	if (p != nullptr)
	{
		skMat frame(h, w);
		BGR2GRAY(p, w, h, frame.data());
		
		for (int i = 0; i < frame.rows() / 2; i++)
		{
			frame.row(i).swap(frame.row(frame.rows() - 1 - i));
		}

		//MatrixXf t = frame.cast<float>();
		//saveTxt("frame.txt", t);

		tracker.Init(SK::skRect(xMin, yMin, width, height), frame);

		//return;

		memset(p, 0, w*h);

		for (int i = 0; i < 64; i++)
		{
			filename = GetImagePath(i);
			LoadJPG(filename.c_str(), p, w, h);
			skMat frame1(h, w);
			BGR2GRAY(p, w, h, frame1.data());
			for (int i = 0; i < frame1.rows() / 2; i++)
			{
				frame1.row(i).swap(frame1.row(frame1.rows() - 1 - i));
			}
			float retvalue;
			SK::skRect result = tracker.update(frame1, retvalue);
			//cout << result.x << "---" << result.y << "---" << result.width << "---" <<
				//result.height << endl;
			framecv = cv::imread(filename);
			cv::Point P1, P2;
			P1.x = result.x;
			P1.y = result.y;
			P2.x = result.x + result.width;
			P2.y = result.y + result.height;
			cv::rectangle(framecv, P1, P2, cv::Scalar(128, 0, 128, 0.3), 2);
			cv::imshow("Result", framecv);
			cv::waitKey(1);
		}
	}
	cout << "hahha" << endl;



	//SK::complex com[8];
	//for (int i = 1; i <= 8; i++)
	//{
	//	//cout << i << "----" << bit_rev(i, 7) << endl;
	//	com[i - 1] = SK::complex(i, 0);
	//	if (i == 8 || i == 7)
	//	{
	//		com[i - 1] = SK::complex(0, 0);
	//	}
	//}
	//cout << "********自定义复数结构*******" << endl;
	//for (int i = 0; i < 8; i++)
	//{
	//	com[i].show();
	//	cout << endl;
	//}
	//cout << endl << endl;

	//cout << "******ffttool 正变换*******" << endl;
	//SK::complex fftArr[8];
	//SK::fft1D(com, fftArr, 8);
	//for (int i = 0; i < 8; i++)
	//{
	//	fftArr[i].show();
	//	cout << endl;
	//}
	//cout << endl << endl;


	//cout << "******ffttool 逆变换*******" << endl;
	//SK::complex ifft[8];
	//SK::ifft1D(fftArr, ifft, 8);
	//for (int i = 0; i < 8; i++)
	//{
	//	ifft[i].show();
	//	cout << endl;
	//}
	//cout << endl << endl;

	//cout << "******fftclass 正变换*******" << endl;
	//vector<complexf> vecf;
	//for (int i = 0; i < 6; i++)
	//{
	//	complexf a(i + 1, i);
	//	vecf.push_back(a);
	//	cout << a << endl;
	//}
	//CFFT fft(6);
	//vector<complexf> res;
	//if (fft.bInitSuccess)
	//{
	//	fft.fft1D(vecf, res);
	//	cout << "res size = " << res.size() << endl;
	//}
	//for (int i = 0; i < res.size(); i++)
	//{
	//	cout << res[i] << endl;
	//}

	//cout << "******fftclass 逆变换*******" << endl;
	//vector<complexf> invres;
	//fft.ifft1D(res, invres);
	//for (int i = 0; i < invres.size(); i++)
	//{
	//	cout << invres[i] << endl;
	//}
}