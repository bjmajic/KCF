#pragma once
#include "RectTolls.h"
//using namespace SK;
namespace SK
{
	class Tracker
	{
	public:
		Tracker();
		virtual ~Tracker();
		virtual void Init(SK::skRect& roi, skMat image) = 0;
		virtual SK::skRect update(skMat image, float& peakValue) = 0;
	protected:
		SK::skRect4f m_roi;
		//skMat        m_grayMat;
	};
}


