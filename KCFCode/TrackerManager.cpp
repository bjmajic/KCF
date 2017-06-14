#include "TrackerManager.h"

TrackerManager::TrackerManager()
{
	
}

TrackerManager::~TrackerManager()
{
	for (size_t i = 0; i < m_vecTrackerPtr.size(); ++i)
	{
		Tracker* p = m_vecTrackerPtr[i];
		delete p;
		p = nullptr;
	}

	if (m_resRect != nullptr)
	{
		delete[] m_resRect;
		m_resRect = nullptr;
	}
}

void TrackerManager::Init(unsigned char* pImage, int width, int height, int num, bool bHog /*= false*/, bool bFixedWindow /*= false*/, bool bMultiscale /*= false*/)
{
	m_bHog = bHog;
	m_bFixedWindow = bFixedWindow;
	m_bMultiscale = bMultiscale;

	m_frame.resize(height*0.25, width*0.25);
	skMat grayMat(height, width);
	BGR2GRAY(pImage, width, height, grayMat.data());

	Zoom(grayMat.data(), width, height, m_frame.data(), width*0.25, height*0.25);
	m_maxTracerNum = num;

	m_resRect = new int[4 * m_maxTracerNum];
}

int TrackerManager::CreateTracker(int x, int y, int w, int h)
{
	if (m_vecTrackerPtr.size() >= m_maxTracerNum)
	{
		return -1; // tracter 数量已达上限，无法create
	}

	KCFTracker* tracker =  new KCFTracker(m_bHog, m_bFixedWindow, m_bMultiscale);
	if (tracker != nullptr)
	{
		tracker->Init(SK::skRect(x*0.25, y*0.25, w*0.25, h*0.25), m_frame);
		m_vecTrackerPtr.push_back(tracker);
		return 1; // 成功
	}
	return 0; // 分配内存失败
}

int* TrackerManager::UpdateTracker(unsigned char* pImage, int width, int height, int& trackerNum)
{
	skMat grayMat(height, width);
	BGR2GRAY(pImage, width, height, grayMat.data());

	Zoom(grayMat.data(), width, height, m_frame.data(), width*0.25, height*0.25);
	float pValue = 0.0f;
	skRect resRect;
	for (size_t i = 0; i < m_vecTrackerPtr.size(); ++i)
	{
		resRect = m_vecTrackerPtr[i]->update(m_frame, pValue);
		m_resRect[i * 4] = resRect.x/0.25;
		m_resRect[i * 4 + 1] = resRect.y/0.25;
		m_resRect[i * 4 + 2] = resRect.width/0.25;
		m_resRect[i * 4 + 3] = resRect.height/0.25;
	}
	trackerNum = m_vecTrackerPtr.size();
	return m_resRect;
}

int* TrackerManager::UpdateTracker(cv::Mat& mat, int& trackerNum)
{
	cv::Mat grayMat;
	cv::cvtColor(mat, grayMat, cv::COLOR_BGR2GRAY);
	cv::resize(grayMat, grayMat, cv::Size(), 0.25, 0.25);

	m_frame = Eigen::Map<skMat>(grayMat.ptr<uchar>(), grayMat.rows, grayMat.cols);
	float pValue = 0.0f;
	skRect resRect;
	for (size_t i = 0; i < m_vecTrackerPtr.size(); ++i)
	{
		resRect = m_vecTrackerPtr[i]->update(m_frame, pValue);
		m_resRect[i * 4] = resRect.x/0.25;
		m_resRect[i * 4 + 1] = resRect.y/0.25;
		m_resRect[i * 4 + 2] = resRect.width/0.25;
		m_resRect[i * 4 + 3] = resRect.height/0.25;

		//just for test
		//cout << "pValue = " << pValue << endl;
	}
	trackerNum = m_vecTrackerPtr.size();
	return m_resRect;
}

void TrackerManager::DeleteTracker(Tracker* pTracker)
{
	vector<Tracker*>::iterator it;
	for (it = m_vecTrackerPtr.begin(); it != m_vecTrackerPtr.end(); )
	{
		if (*it == pTracker)
		{
			Tracker* p = pTracker;
			it = m_vecTrackerPtr.erase(it);

			//释放内存
			delete p;
			p = nullptr;

			break;
		}
	}
}
