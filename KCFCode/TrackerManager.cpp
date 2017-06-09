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

	m_frame.resize(height, width);
	BGR2GRAY(pImage, width, height, m_frame.data());

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
		tracker->Init(SK::skRect(x, y, w, h), m_frame);
		m_vecTrackerPtr.push_back(tracker);
		return 1; // 成功
	}
	return 0; // 分配内存失败
}

int* TrackerManager::UpdateTracker(unsigned char* pImage, int width, int height, int& trackerNum)
{
	BGR2GRAY(pImage, width, height, m_frame.data());
	float pValue = 0.0f;
	skRect resRect;
	for (size_t i = 0; i < m_vecTrackerPtr.size(); ++i)
	{
		resRect = m_vecTrackerPtr[i]->update(m_frame, pValue);
		m_resRect[i * 4] = resRect.x;
		m_resRect[i * 4 + 1] = resRect.y;
		m_resRect[i * 4 + 2] = resRect.width;
		m_resRect[i * 4 + 3] = resRect.height;
	}
	trackerNum = m_vecTrackerPtr.size();
	return m_resRect;
}

int* TrackerManager::UpdateTracker(cv::Mat& mat, int& trackerNum)
{
	cv::Mat grayMat;
	cv::cvtColor(mat, grayMat, cv::COLOR_BGR2GRAY);

	m_frame = Eigen::Map<skMat>(grayMat.ptr<uchar>(), grayMat.rows, grayMat.cols);
	float pValue = 0.0f;
	skRect resRect;
	for (size_t i = 0; i < m_vecTrackerPtr.size(); ++i)
	{
		resRect = m_vecTrackerPtr[i]->update(m_frame, pValue);
		m_resRect[i * 4] = resRect.x;
		m_resRect[i * 4 + 1] = resRect.y;
		m_resRect[i * 4 + 2] = resRect.width;
		m_resRect[i * 4 + 3] = resRect.height;
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
