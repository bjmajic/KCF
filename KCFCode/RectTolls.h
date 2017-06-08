#ifndef SK_RECTTOOLS_HPP
#define SK_RECTTOOLS_HPP
#include "Eigen/Dense"
#include "ImageJpeg.h"
#include <iostream>
#include <fstream>
using namespace Eigen;
typedef Matrix<unsigned char, Dynamic, Dynamic, RowMajor> skMat;

// just for debug, record the res to txt

inline void saveTxt(string savedName, MatrixXf matdata)
{
	string basePath = "D:\\KCF\\";
	basePath = basePath + savedName;

	ofstream of1(basePath);
	of1.setf(ios::fixed);
	of1.precision(4);
	int height = matdata.rows();
	int width = matdata.cols();

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			of1 << matdata(i, j) << "\t";
		}
		of1 << "\n";
	}
}

inline void saveTxt(string savedName, MatrixXcf matdata)
{
	string basePath = "D:\\KCF\\";
	basePath = basePath + savedName;

	ofstream of1(basePath);
	//of1.setf(ios::fixed);
	//of1.precision(4);
	int height = matdata.rows();
	int width = matdata.cols();

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			of1 << matdata(i, j) << "\t";
		}
		of1 << "\n";
	}
}

inline void saveTxt(string savedName, skMat matdata)
{
	string basePath = "D:\\KCF\\";
	basePath = basePath + savedName;

	ofstream of1(basePath);
	//of1.setf(ios::fixed);
	//of1.precision(4);
	int height = matdata.rows();
	int width = matdata.cols();

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			of1 << matdata(i, j) << "\t";
		}
		of1 << "\n";
	}
}

namespace SK
{
	template<typename T>
	struct Rect_tag
	{
		T x;
		T y;
		T width;
		T height;
		Rect_tag()
		{
			x = 0;
			y = 0;
			width = 0;
			height = 0;
		}
		Rect_tag(T x_, T y_, T width_, T height_)
		{
			x = x_;
			y = y_;
			width = width_;
			height = height_;
		}

		friend bool operator != (Rect_tag<T>& rect1, Rect_tag<T>& rect2)
		{
			if (rect1.x == rect2.x && rect1.y == rect2.y && rect1.width == rect2.width &&
				rect1.height == rect2.height)
			{
				return false;
			}
			else
			{
				return true;
			}
		}

	};
	template<typename T1, typename T2>
	void RectCopy(Rect_tag<T1>& srcRect, Rect_tag<T2>& desRect)
	{
		desRect.x = static_cast<T2>(srcRect.x);
		desRect.y = static_cast<T2>(srcRect.y);
		desRect.width = static_cast<T2>(srcRect.width);
		desRect.height = static_cast<T2>(srcRect.height);
	}
	typedef Rect_tag<float> skRect4f;
	typedef Rect_tag<int> skRect;

	//定义两个辅助结构
	template <typename T>
	struct skSize
	{
		T nWidth;
		T nHeight;
		skSize(T w, T h)
		{
			nWidth = w;
			nHeight = h;
		}
		skSize()
		{
			nWidth = 0;
			nHeight = 0;
		}
	};
	
	template<typename T>
	struct skPoint
	{
		T x;
		T y;

		skPoint(T x_, T y_)
		{
			x = x_;
			y = y_;
		}

		skPoint()
		{
			x = 0; y = 0;
		}
	};

	typedef skPoint<float> skPoint2f;

	//获取rect的中心
	template<typename T>
	inline skPoint<T> center(const Rect_tag<T>& rect)
	{
		return skPoint<T>(rect.x + rect.width*0.5, rect.y + rect.height*0.5);
	}

	//获取rect右下角坐标（x, y)
	template<typename T>
	inline T x2(const Rect_tag<T> &rect)
	{
		return rect.x + rect.width;
	}
	template<typename T>
	inline T y2(const Rect_tag<T> &rect)
	{
		return rect.y + rect.height;
	}

	template <typename t>
	inline void resize(Rect_tag<t> &rect, float scalex, float scaley = 0)
	{
		if (!scaley)scaley = scalex;
		rect.x -= rect.width * (scalex - 1.f) / 2.f;
		rect.width *= scalex;

		rect.y -= rect.height * (scaley - 1.f) / 2.f;
		rect.height *= scaley;
	}
	
	// clip Rect to limit size
	template <typename t>
	inline void limit(Rect_tag<t> &rect, Rect_tag<t> limit)
	{
		if (rect.x + rect.width > limit.x + limit.width) rect.width = (limit.x + limit.width - rect.x);
		if (rect.y + rect.height > limit.y + limit.height) rect.height = (limit.y + limit.height - rect.y);
		if (rect.x < limit.x)
		{
			rect.width -= (limit.x - rect.x);
			rect.x = limit.x;
		}
		if (rect.y < limit.y)
		{
			rect.height -= (limit.y - rect.y);
			rect.y = limit.y;
		}
		if (rect.width < 0)rect.width = 0;
		if (rect.height < 0)rect.height = 0;
	}

	template <typename t>
	inline void limit(Rect_tag<t> &rect, t width, t height, t x = 0, t y = 0)
	{
		limit(rect, Rect_tag<t>(x, y, width, height));
	}


	// 获取边框的尺寸（其实就是上下左右四个边各填充多少像素）
	template <typename t>
	inline skRect getBorder(const Rect_tag<t> &original, Rect_tag<t> &limited)
	{
		skRect res;
		res.x = static_cast<int>(limited.x - original.x);
		res.y = static_cast<int>(limited.y - original.y);
		res.width = static_cast<int>(x2(original) - x2(limited));
		res.height = static_cast<int>(y2(original) - y2(limited));
		assert(res.x >= 0 && res.y >= 0 && res.width >= 0 && res.height >= 0);
		return res;
	}

	inline void copyMakeBorder(const skMat &src, skMat &dst, int top, int bottom, int left, int right, int boardtype = 0)
	{
		if (0 == boardtype) // 用0填充
		{
			//暂时只能想到以下这个方案，1、声明一个比较大的矩阵, 2、将原矩阵的内容赋值过来，3、填充边缘 4、赋值给dst
			skMat bigSize;
			bigSize.resize(src.rows() + top + bottom, src.cols() + left + right);
			bigSize.fill(0);
			bigSize.block(top, left, src.rows(), src.cols()) = src;
			dst = bigSize;
		}
	}

	inline skMat subwindow(const skMat &in, const skRect& window/*, int borderType = cv::BORDER_CONSTANT*/)
	{
		skRect cutWindow = window;
		limit<int>(cutWindow, in.cols(), in.rows());
		if (cutWindow.height <= 0 || cutWindow.width <= 0)assert(0); //return cv::Mat(window.height,window.width,in.type(),0) ;
		skRect border = getBorder(window, cutWindow);

		//MatrixXf res = in(cutWindow);
		skMat res = in.block(cutWindow.y, cutWindow.x, cutWindow.height, cutWindow.width);

		if (border != skRect(0, 0, 0, 0))
		{
			copyMakeBorder(res, res, border.y, border.height, border.x, border.width);
		}
		return res;
	}
	inline MatrixXf getGrayImage(const skMat& img)
	{
		MatrixXf resMatrix = img.cast<float>();
		//resMatrix *= 0.003922f;
		resMatrix /= 255;
		return resMatrix;
	}

	inline void Resize(skMat& z, const skSize<int>& newSize)
	{
		skMat newZ(z.rows(), z.cols());
		Zoom(z.data(), z.cols(), z.rows(), newZ.data(), newSize.nHeight, newSize.nHeight);
		z = newZ;
	}

	inline void rearrange(MatrixXf &img)
	{
		// img = img(cv::Rect(0, 0, img.cols & -2, img.rows & -2));
		int cx = img.cols() / 2;
		int cy = img.rows() / 2;

		//MatrixXf q0 = img.block(0, 0, cy, cx);
		//MatrixXf q1 = img.block(0, cx, cy, cx);
		//MatrixXf q2 = img.block(cy, 0, cy, cx);
		//MatrixXf q3 = img.block(cy, cx, cy, cx);

		img.block(0, 0, cy, cx).swap(img.block(cy, cx, cy, cx));
		img.block(0, cx, cy, cx).swap(img.block(cy, 0, cy, cx));

		//做这步的目的是讲频域的原点移动到图像的中心，详见论文附录A.1
		//saveTxt("rearrange_img.txt", img);
	}
}
#endif //SK_RECTTOOLS_HPP