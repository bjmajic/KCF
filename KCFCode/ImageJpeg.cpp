#include "ImageJpeg.h"
#include <cstdio>

bool LoadJPG(const char *filename, unsigned char *&outbuf, int &width, int &height)
{
	FILE *fimg;
	fimg = fopen(filename, "rb");
	if (fimg == NULL) return false;
	int jpglenght = 0;
	fseek(fimg, 0, SEEK_END);
	jpglenght = ftell(fimg);
	fseek(fimg, 0, SEEK_SET);
	unsigned char *jpgData = new unsigned char[jpglenght];
	fread(jpgData, 1, jpglenght, fimg);
	fclose(fimg);
	bool ret = jpeg2bmp_decompress2(jpgData, jpglenght, &outbuf, &width, &height);
	delete[]jpgData;
	return ret;

}

void Zoom(unsigned char *pSrcImg, long nSrcWidth, long nSrcHeight, unsigned char *pDstImg, long nDstWidth, long nDstHeight)
{
	int		n_x_d, n_y_d;
	int		n_x_s, n_y_s;
	float	lfXscl, lfYScl, lf_x_s, lf_y_s, lfNewB, lfNewG, lfNewR;
	float	lfWeight_x, lfWeight_y;

	if (nSrcWidth == nDstWidth && nSrcHeight == nDstHeight)
	{
		//memcpy(pDstImg, pSrcImg, nSrcWidth*nSrcHeight * 3);
		memcpy(pDstImg, pSrcImg, nSrcWidth*nSrcHeight);
		return;
	}

	lfXscl = float(nSrcWidth) / nDstWidth;
	lfYScl = float(nSrcHeight) / nDstHeight;

	nSrcWidth = nSrcWidth;
	nDstWidth = nDstWidth;

	for (n_y_d = 0; n_y_d < nDstHeight; n_y_d++)
	{
		for (n_x_d = 0; n_x_d < nDstWidth; n_x_d++)
		{
			lf_x_s = lfXscl * (n_x_d + 1) - 1;
			lf_y_s = lfYScl * (n_y_d + 1) - 1;
			n_x_s = int(lf_x_s);
			n_x_s = SMALLER(n_x_s, nSrcWidth - 2);
			n_y_s = int(lf_y_s);
			n_y_s = SMALLER(n_y_s, nSrcHeight - 2);
			lfWeight_x = lf_x_s - n_x_s;
			lfWeight_y = lf_y_s - n_y_s;

			lfNewB = (1 - lfWeight_y) * ((1 - lfWeight_x) * pSrcImg[(n_y_s * nSrcWidth + n_x_s)]
				+ lfWeight_x * pSrcImg[(n_y_s * nSrcWidth + n_x_s + 1)])
				+ lfWeight_y * ((1 - lfWeight_x) * pSrcImg[((n_y_s + 1) * nSrcWidth + n_x_s)]
				+ lfWeight_x * pSrcImg[((n_y_s + 1) * nSrcWidth + n_x_s + 1)]);
			pDstImg[(n_y_d * nDstWidth + n_x_d)] = unsigned char(lfNewB);

			/*lfNewB = (1 - lfWeight_y) * ((1 - lfWeight_x) * pSrcImg[(n_y_s * nSrcWidth + n_x_s) * 3 + 0]
				+ lfWeight_x * pSrcImg[(n_y_s * nSrcWidth + n_x_s + 1) * 3 + 0])
				+ lfWeight_y * ((1 - lfWeight_x) * pSrcImg[((n_y_s + 1) * nSrcWidth + n_x_s) * 3 + 0]
				+ lfWeight_x * pSrcImg[((n_y_s + 1) * nSrcWidth + n_x_s + 1) * 3 + 0]);
				pDstImg[(n_y_d * nDstWidth + n_x_d) * 3 + 0] = unsigned char(lfNewB);

				lfNewG = (1 - lfWeight_y) * ((1 - lfWeight_x) * pSrcImg[(n_y_s * nSrcWidth + n_x_s) * 3 + 1]
				+ lfWeight_x * pSrcImg[(n_y_s * nSrcWidth + n_x_s + 1) * 3 + 1])
				+ lfWeight_y * ((1 - lfWeight_x) * pSrcImg[((n_y_s + 1) * nSrcWidth + n_x_s) * 3 + 1]
				+ lfWeight_x * pSrcImg[((n_y_s + 1) * nSrcWidth + n_x_s + 1) * 3 + 1]);
				pDstImg[(n_y_d * nDstWidth + n_x_d) * 3 + 1] = unsigned char(lfNewG);

				lfNewR = (1 - lfWeight_y) * ((1 - lfWeight_x) * pSrcImg[(n_y_s * nSrcWidth + n_x_s) * 3 + 2]
				+ lfWeight_x * pSrcImg[(n_y_s * nSrcWidth + n_x_s + 1) * 3 + 2])
				+ lfWeight_y * ((1 - lfWeight_x) * pSrcImg[((n_y_s + 1) * nSrcWidth + n_x_s) * 3 + 2]
				+ lfWeight_x * pSrcImg[((n_y_s + 1) * nSrcWidth + n_x_s + 1) * 3 + 2]);
				pDstImg[(n_y_d * nDstWidth + n_x_d) * 3 + 2] = unsigned char(lfNewR);*/
		}
	}

}

void BGR2GRAY(unsigned char *scr, long width, long height, unsigned char *dst)
{
	int i, j;
	int b, g, r, t;
	int index1 = 0;
	int index2 = 0;
	for (j = 0; j < height; j++)
	{
		for (i = 0; i < width; i++)
		{
			b = (*(scr + index1))*cscGb;
			g = (*(scr + index1 + 1))*cscGg;
			r = (*(scr + index1 + 2))*cscGr;
			t = b + g + r;
			t = descale(t, shift);
			*(dst + index2) = fast_cast_8u(t);

			index1 += 3;
			index2 += 1;
		}
	}
}
