#pragma once
#include "RectTolls.h"
using namespace SK;
class Tracker
{
public:
	Tracker();
	virtual ~Tracker();
	virtual void Init(Rect& roi, Mat image) = 0;
	virtual Rect update(Mat image) = 0;
protected:
	Rect4f m_roi; 
};

