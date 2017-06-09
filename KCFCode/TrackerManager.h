#pragma once
#include <opencv2/opencv.hpp>
#include "KCFTracker.h"
#include "Tracker.h"

using namespace SK;
class TrackerManager
{
public:
	TrackerManager();
	~TrackerManager();
	
	void Init(unsigned char* pImage, int width, int height, int num, bool bHog = false, bool bFixedWindow = false, bool bMultiscale = false);
	int CreateTracker(int x, int y, int w, int h);
	int* UpdateTracker(unsigned char* pImage, int width, int height, int& trackerNum);
	void DeleteTracker(Tracker* pTracker);

	// just for opencv test
	int* UpdateTracker(cv::Mat& mat, int& trackerNum);

	int* m_resRect;

	
private:
	bool m_bHog;
	bool m_bFixedWindow;
	bool m_bMultiscale;
	skMat m_frame;
	int  m_maxTracerNum;
	vector<Tracker*> m_vecTrackerPtr;

	//void Reset()
};

