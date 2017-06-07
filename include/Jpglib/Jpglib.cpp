//JPG 格式转换库
//#include <iostream>
//#include <stdio.h>

extern "C"
{
	#include "jpeglib.h"
}
#include "Jpglib.h"
#include "jmyerror.h"

#pragma pack(push)
#pragma pack(2)        //两字节对齐，否则bmp_fileheader会占16Byte
typedef struct bmp_fileheader
{
	unsigned short    bfType;        //若不对齐，这个会占4Byte
	unsigned long     bfSize;
	unsigned short    bfReverved1;
	unsigned short    bfReverved2;
	unsigned long     bfOffBits;
}bmp_fileheader;
typedef struct bmp_infoheader
{
	unsigned long    biSize;
	unsigned long    biWidth;
	unsigned long    biHeight;
	unsigned short   biPlanes;
	unsigned short   biBitCount;
	unsigned long    biCompression;
	unsigned long    biSizeImage;
	unsigned long    biXPelsPerMeter;
	unsigned long    biYpelsPerMeter;
	unsigned long    biClrUsed;
	unsigned long    biClrImportant;
}bmp_infoheader;
#pragma pack(pop)


void write_bmp_header(j_decompress_ptr cinfo, FILE* output_file);	//写入bmp文件头
void write_bmp_data(j_decompress_ptr cinfo, unsigned char *src_buff, FILE* output_file);	//写入bmp文件数据

void write_bmp_header(j_decompress_ptr cinfo, FILE* output_file)
{
	struct bmp_fileheader bfh;
	struct bmp_infoheader bih;

	unsigned long width;
	unsigned long height;
	unsigned short depth;
	unsigned long headersize;
	unsigned long filesize;

	width = cinfo->output_width;
	height = cinfo->output_height;
	depth = cinfo->output_components;

	if (depth == 1)
	{
		headersize = 14 + 40 + 256 * 4;
		filesize = headersize + width*height;
	}

	if (depth == 3)
	{
		headersize = 14 + 40;
		filesize = headersize + width*height*depth;
	}

	memset(&bfh, 0, sizeof(struct bmp_fileheader));
	memset(&bih, 0, sizeof(struct bmp_infoheader));

	//写入比较关键的几个bmp头参数
	bfh.bfType = 0x4D42;
	bfh.bfSize = filesize;
	bfh.bfOffBits = headersize;

	bih.biSize = 40;
	bih.biWidth = width;
	bih.biHeight = height;
	bih.biPlanes = 1;
	bih.biBitCount = (unsigned short)depth * 8;
	bih.biSizeImage = width*height*depth;

	fwrite(&bfh, sizeof(struct bmp_fileheader), 1, output_file);
	fwrite(&bih, sizeof(struct bmp_infoheader), 1, output_file);

	if (depth == 1)        //灰度图像要添加调色板
	{
		unsigned char *platte;
		platte = (unsigned char*)malloc(256 * 4);
		unsigned char j = 0;
		int i;
		for (i = 0; i<1024; i+=4)
		{
			platte[i] = j;
			platte[i + 1] = j;
			platte[i + 2] = j;
			platte[i + 3] = 0;
			j++;
		}
		fwrite(platte, sizeof(unsigned char)* 1024, 1, output_file);
		free(platte);
	}
}


void write_bmp_data(j_decompress_ptr cinfo, unsigned char *src_buff, FILE* output_file)
{
	unsigned char *dst_width_buff;
	unsigned char *point;

	unsigned long width;
	unsigned long height;
	unsigned short depth;

	width = cinfo->output_width;
	height = cinfo->output_height;
	depth = cinfo->output_components;

	dst_width_buff = (unsigned char*)malloc(width*depth);
	memset(dst_width_buff, 0, sizeof(unsigned char)*width*depth);

	point = src_buff + width*depth*(height - 1);    //倒着写数据，bmp格式是倒的，jpg是正的
	unsigned long i;
	for (i = 0; i<height; i++)
	{
		unsigned long j;
		for (j = 0; j<width*depth; j += depth)
		{
			if (depth == 1)        //处理灰度图
			{
				dst_width_buff[j] = point[j];
			}

			if (depth == 3)        //处理彩色图
			{
				dst_width_buff[j + 2] = point[j + 0];
				dst_width_buff[j + 1] = point[j + 1];
				dst_width_buff[j + 0] = point[j + 2];
			}
		}
		point -= width*depth;
		fwrite(dst_width_buff, sizeof(unsigned char)*width*depth, 1, output_file);    //一次写一行
	}
	free(dst_width_buff);
}

void bmp2jpeg_compress(unsigned char *bgr, int width, int height, unsigned char **outbuf, unsigned long *outSize, int quality)
{
	int jpegWidth = width;//jpeg的宽度;
	int jpegHeight = height;//jpeg的高度;
	//开始进行jpg的数据写入
	struct jpeg_compress_struct toWriteCinfo;
	struct jpeg_error_mgr jerr;
	JSAMPROW row_pointer[1];
	int row_stride;


	toWriteCinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&toWriteCinfo);
	//确定要用于输出压缩的jpeg的数据空间
	jpeg_mem_dest(&toWriteCinfo, outbuf, outSize);


	toWriteCinfo.image_width = jpegWidth;
	toWriteCinfo.image_height = jpegHeight;
	toWriteCinfo.input_components = 3;
	toWriteCinfo.in_color_space = JCS_RGB;


	jpeg_set_defaults(&toWriteCinfo);
	jpeg_set_quality(&toWriteCinfo, quality, TRUE);
	jpeg_start_compress(&toWriteCinfo, TRUE);
	row_stride = jpegWidth*3;//如果图片为RGB，这个值要*3.灰度图像不用。

	int image_buffer_len = row_stride * jpegHeight; // DIB buffer 长度
	while (toWriteCinfo.next_scanline < toWriteCinfo.image_height)
	{
		//row_pointer[0] = &bgr[toWriteCinfo.next_scanline * row_stride];
		row_pointer[0] = &bgr[image_buffer_len - (toWriteCinfo.next_scanline + 1) * row_stride];
		(void)jpeg_write_scanlines(&toWriteCinfo, row_pointer, 1);
	}
	jpeg_finish_compress(&toWriteCinfo);
	jpeg_destroy_compress(&toWriteCinfo);
}


void jpg_to_bmp(const char* jpg_file, const char* bmp_file)
{
	FILE *input_file;
	FILE *output_file;
	input_file = fopen(jpg_file, "rb");
	output_file = fopen(bmp_file, "wb");
	if (NULL == input_file || NULL == output_file){
		printf("Func <jpg_to_bmp> open file ERR!!!\n\n");
		return;
	}

	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	JSAMPARRAY buffer;
	unsigned char *src_buff;
	unsigned char *point;

	cinfo.err = jpeg_std_error(&jerr);    //一下为libjpeg函数，具体参看相关文档
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, input_file);
	jpeg_read_header(&cinfo, TRUE);
	jpeg_start_decompress(&cinfo);

	unsigned long width = cinfo.output_width;
	unsigned long height = cinfo.output_height;
	unsigned short depth = cinfo.output_components;

	src_buff = (unsigned char*)malloc(width*height*depth);
	memset(src_buff, 0, sizeof(unsigned char)*width*height*depth);

	buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr)&cinfo, JPOOL_IMAGE, width*depth, 1);

	point = src_buff;
	while (cinfo.output_scanline<height)
	{
		jpeg_read_scanlines(&cinfo, buffer, 1);    //读取一行jpg图像数据到buffer
		memcpy(point, *buffer, width*depth);    //将buffer中的数据逐行给src_buff
		point += width*depth;            //一次改变一行
	}

	write_bmp_header(&cinfo, output_file);            //写bmp文件头
	write_bmp_data(&cinfo, src_buff, output_file);    //写bmp像素数据

	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	free(src_buff);

	fclose(input_file);
	fclose(output_file);
}

void bmp_to_jpg(const char* bmp_file, const char* jpg_file, int quality)
{
	FILE *input_file;
	FILE *output_file;
	input_file = fopen(bmp_file, "rb");
	output_file = fopen(jpg_file, "wb");
	if (NULL == input_file || NULL == output_file){
		printf("Func <jpg_to_bmp> open file ERR!!!\n\n");
		return;
	}

	//读取bmp图像信息
	bmp_fileheader file_header;
	bmp_infoheader info_header;
	fread(&file_header, sizeof(file_header), 1, input_file);
	fread(&info_header, sizeof(info_header), 1, input_file);

	if (24 != info_header.biBitCount){		//只处理24位情况
		printf("Read bmp file is not 24_bit_count!!!\n");
		fclose(input_file);
		fclose(output_file);
		return;
	}

	// 每行像素占用的字节数,每行要对齐4字节.
	int bmp_width = info_header.biWidth;
	int bmp_hight = info_header.biHeight;
	//int imageRowSize = (bmp_width * 24 + 31) / 32 * 4;
	int imageRowSize = ((bmp_width * 3 + 3) / 4) * 4;
	unsigned char* dibBuffer = (unsigned char*)malloc(imageRowSize * info_header.biWidth);
	if (NULL == dibBuffer){
		printf("New dibBuffer ERROR!\n");
		return;
	}
	memset(dibBuffer, 0, imageRowSize * bmp_hight);
	size_t rCount = fread(dibBuffer, imageRowSize*bmp_hight, 1, input_file);
	printf("Read count: %d\n", rCount);

	// DIB 中颜色的存放顺序是 BGR, 而 JPG 要求的顺序是 RGB, 所以要交换 R 和 B.
	// 由于有行对齐因素,所以逐行处理
	int row;
	for (row = 0; row < bmp_hight; ++row)
	{
		unsigned char* rowData = dibBuffer + imageRowSize * row;
		int col;
		for (col = 0; col < bmp_width * 3; col += 3)
		{
			unsigned char swap = rowData[col];
			rowData[col] = rowData[col + 2];
			rowData[col + 2] = swap;
		}
	}
	unsigned char *outbuf = NULL;
	unsigned long outSize = 0;
	//把位图数据压缩为 jpeg
	bmp2jpeg_compress(dibBuffer, bmp_width, bmp_hight, &outbuf, &outSize, quality);
	fwrite(outbuf, 1, outSize, output_file);

	free(outbuf);

	fclose(output_file);
	fclose(input_file);
	// 释放 DIB 数组
	free(dibBuffer);
}



/*
void bmp_to_jpg(const char* bmp_file, const char* jpg_file, int quality)
{
	FILE *input_file;
	FILE *output_file;
	input_file = fopen(bmp_file, "rb");
	output_file = fopen(jpg_file, "wb");
	if (NULL == input_file || NULL == output_file){
		printf("Func <jpg_to_bmp> open file ERR!!!\n\n");
		return;
	}

	//读取bmp图像信息
	bmp_fileheader file_header;
	bmp_infoheader info_header;
	fread(&file_header, sizeof(file_header), 1, input_file);
	fread(&info_header, sizeof(info_header), 1, input_file);

	if(24!=info_header.biBitCount){		//只处理24位情况
		printf("Read bmp file is not 24_bit_count!!!\n");
		fclose(input_file);
		fclose(output_file);
		return;
	}

	// 每行像素占用的字节数,每行要对齐4字节.
	int bmp_width = info_header.biWidth;
	int bmp_hight = info_header.biHeight;
	//int imageRowSize = (bmp_width * 24 + 31) / 32 * 4;
	int imageRowSize = ((bmp_width * 3 + 3)/4)*4;
	unsigned char* dibBuffer = new unsigned char[imageRowSize * info_header.biWidth];
	if (NULL == dibBuffer){
		printf("New dibBuffer ERROR!\n");
		return;
	}
	memset(dibBuffer, 0, imageRowSize * bmp_hight);
	size_t rCount=fread(dibBuffer, imageRowSize*bmp_hight, 1, input_file);
	printf("Read count: %d\n", rCount);

	// DIB 中颜色的存放顺序是 BGR, 而 JPG 要求的顺序是 RGB, 所以要交换 R 和 B.
	// 由于有行对齐因素,所以逐行处理
	for (int row = 0; row < bmp_hight; ++row)
	{
		unsigned char* rowData = dibBuffer + imageRowSize * row;
		for (int col = 0; col < bmp_width * 3; col += 3)
		{
			unsigned char swap = rowData[col];
			rowData[col] = rowData[col + 2];
			rowData[col + 2] = swap;
		}
	}

	//把位图数据压缩为 jpeg
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	//FILE * outfile;  // target file
	JSAMPROW row_pointer[1]; // pointer to JSAMPLE row[s]
	int row_stride;  // physical row width in image buffer
	int image_width = bmp_width;
	int image_height = bmp_hight;
	JSAMPLE* image_buffer = dibBuffer;     // DIB buffer
	int image_buffer_len = imageRowSize * image_height; // DIB buffer 长度


	{
		// Step 1: allocate and initialize JPEG compression object
		cinfo.err = jpeg_std_error(&jerr);

		// Now we can initialize the JPEG compression object.
		jpeg_create_compress(&cinfo);


		// Step 2: specify data destination (eg, a file)
		// Note: steps 2 and 3 can be done in either order.
		jpeg_stdio_dest(&cinfo, output_file);

		// Step 3: set parameters for compression
		// First we supply a description of the input image.
		// Four fields of the cinfo struct must be filled in:

		cinfo.image_width = image_width;  // image width and height, in pixels
		cinfo.image_height = image_height;
		cinfo.input_components = 3;  // # of color components per pixel  // 因为DIB数据是24位的,所以每个像素占用3个字节
		cinfo.in_color_space = JCS_RGB;  // colorspace of input image
		// Now use the library's routine to set default compression parameters.
		// (You must set at least cinfo.in_color_space before calling this,
		// since the defaults depend on the source color space.)

		jpeg_set_defaults(&cinfo);
		// Now you can set any non-default parameters you wish to.
		// Here we just illustrate the use of quality (quantization table) scaling:
		//
		jpeg_set_quality(&cinfo, quality, TRUE );// limit to baseline-JPEG values

		// Step 4: Start compressor
		// TRUE ensures that we will write a complete interchange-JPEG file.
		// Pass TRUE unless you are very sure of what you're doing.
		//
		jpeg_start_compress(&cinfo, TRUE);

		// Step 5: while (scan lines remain to be written)
		//           jpeg_write_scanlines(...);

		// Here we use the library's state variable cinfo.next_scanline as the
		// loop counter, so that we don't have to keep track ourselves.
		// To keep things simple, we pass one scanline per call; you can pass
		// more if you wish, though.

		row_stride = imageRowSize;
		while (cinfo.next_scanline < cinfo.image_height)
		{
			// jpeg_write_scanlines expects an array of pointers to scanlines.
			// Here the array is only one element long, but you could pass
			// more than one scanline at a time if that's more convenient.

			//row_pointer[0] = &image_buffer[cinfo.next_scanline * row_stride];
			row_pointer[0] = &image_buffer[image_buffer_len - (cinfo.next_scanline + 1) * row_stride];
			(void)jpeg_write_scanlines(&cinfo, row_pointer, 1);
		}
		// Step 6: Finish compression
		jpeg_finish_compress(&cinfo);

		// After finish_compress, we can close the output file.

		// Step 7: release JPEG compression object
		// This is an important step since it will release a good deal of memory.
		jpeg_destroy_compress(&cinfo);
	}

	fclose(output_file);
	fclose(input_file);
	// 释放 DIB 数组
	delete[]dibBuffer;
}


*/

bool jpeg2bmp_decompress(unsigned char *jpeg, int jpeg_size, unsigned char **outbuf, int *width, int *height)
{
	struct jpeg_decompress_struct cinfo;
	//struct jpeg_error_mgr jerr;
	struct my_error_mgr jerr;
	*width = 0;
	*height = 0;
	*outbuf = 0;

	//	JSAMPARRAY buffer;
	int row_stride = 0;
	unsigned char* tmp_buffer = NULL;
	int rgb_size;

	if (jpeg == NULL)
	{
		printf("no jpeg buffer here.\n");
		return false;
	}
	cinfo.err = jpeg_std_error(&jerr.pub);
	jerr.pub.error_exit = my_error_exit;
	if (setjmp(jerr.setjmp_buffer)) {
		jpeg_destroy_decompress(&cinfo);
		return false;
	}

	jpeg_create_decompress(&cinfo);
	jpeg_mem_src(&cinfo, jpeg, (unsigned long)(jpeg_size));

	jpeg_read_header(&cinfo, TRUE);
	//cinfo.out_color_space = JCS_RGB; //JCS_YCbCr;  // 设置输出格式 
	//cinfo.output_components = 3; //输出3个通道
	jpeg_start_decompress(&cinfo);

	row_stride = cinfo.output_width * cinfo.output_components;
	*width = cinfo.output_width;
	*height = cinfo.output_height;
	rgb_size = row_stride * cinfo.output_height; // 总大小 
	unsigned char *rgb_buffer = (unsigned char *)malloc(rgb_size + 54);
	memset(rgb_buffer, 0, sizeof(unsigned char)*rgb_size + 54);
	tmp_buffer = rgb_buffer;
	*outbuf = rgb_buffer;
	unsigned char * buffer[1];
	while (cinfo.output_scanline < cinfo.output_height) // 解压每一行  
	{
		buffer[0] = tmp_buffer;
		jpeg_read_scanlines(&cinfo, buffer, 1);
		tmp_buffer += row_stride;
	}

	for (int i = 0; i < (int)(cinfo.output_width*cinfo.output_height); i++)
	{
		unsigned char tmp;
		tmp = rgb_buffer[0];
		rgb_buffer[0] = rgb_buffer[2];
		rgb_buffer[2] = tmp;
		rgb_buffer += 3;
	}

	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	return true;
}


bool SaveRGB24JPG(void *img, int width, int height, const char* filename, int quality)
{
	FILE *output_file;
	output_file = fopen(filename, "wb");
	if (NULL == output_file) return false;

	// 每行像素占用的字节数,每行要对齐4字节.
	int bmp_width = width;
	int bmp_hight = height;
	int imageRowSize = bmp_width * 3;
	unsigned char* dibBuffer = (unsigned char*)malloc(imageRowSize * height);
	if (NULL == dibBuffer){
		fclose(output_file);
		return false;
	}

	//memcpy(dibBuffer, img , imageRowSize * bmp_hight);

	// DIB 中颜色的存放顺序是 BGR, 而 JPG 要求的顺序是 RGB, 所以要交换 R 和 B.
	// 由于有行对齐因素,所以逐行处理
	for (int row = 0; row < bmp_hight; ++row)
	{
		unsigned char* rowData = dibBuffer + imageRowSize * row;
		memcpy(rowData, (char*)img + (bmp_hight - row - 1) * imageRowSize, imageRowSize);
		for (int col = 0; col < bmp_width * 3; col += 3)
		{
			unsigned char swap = rowData[col];
			rowData[col] = rowData[col + 2];
			rowData[col + 2] = swap;
		}
	}

	//把位图数据压缩为 jpeg
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	//FILE * outfile;  /* target file */
	JSAMPROW row_pointer[1]; /* pointer to JSAMPLE row[s] */
	int row_stride;  /* physical row width in image buffer */
	int image_width = bmp_width;
	int image_height = bmp_hight;
	JSAMPLE* image_buffer = dibBuffer;     // DIB buffer
	int image_buffer_len = imageRowSize * image_height; // DIB buffer 长度



	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, output_file);

	cinfo.image_width = image_width;
	cinfo.image_height = image_height;
	cinfo.input_components = 3; // 因为DIB数据是24位的,所以每个像素占用3个字节
	cinfo.in_color_space = JCS_RGB;
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality, TRUE);
	jpeg_start_compress(&cinfo, TRUE);
	row_stride = imageRowSize;
	while (cinfo.next_scanline < cinfo.image_height)
	{
		row_pointer[0] = &image_buffer[cinfo.next_scanline * row_stride];
		(void)jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}
	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);

	fclose(output_file);
	free(dibBuffer);
	return true;
}

bool jpeg2bmp_decompress2(unsigned char *jpeg, int jpeg_size, unsigned char **outbuf, int *width, int *height)
{
	struct jpeg_decompress_struct cinfo;
	//struct jpeg_error_mgr jerr;
	struct my_error_mgr jerr;
	*width = 0;
	*height = 0;
	*outbuf = 0;

	//	JSAMPARRAY buffer;
	int row_stride = 0;
	unsigned char* tmp_buffer = NULL;
	int rgb_size;

	if (jpeg == NULL)
	{
		printf("no jpeg buffer here.\n");
		return false;
	}
	cinfo.err = jpeg_std_error(&jerr.pub);
	jerr.pub.error_exit = my_error_exit;
	if (setjmp(jerr.setjmp_buffer)) {
		jpeg_destroy_decompress(&cinfo);
		return false;
	}

	jpeg_create_decompress(&cinfo);
	jpeg_mem_src(&cinfo, jpeg, (unsigned long)(jpeg_size));

	jpeg_read_header(&cinfo, TRUE);
	//cinfo.out_color_space = JCS_RGB; //JCS_YCbCr;  // 设置输出格式 
	//cinfo.output_components = 3; //输出3个通道
	jpeg_start_decompress(&cinfo);
	
	row_stride = cinfo.output_width * cinfo.output_components;
	int rgbrow_stride = cinfo.output_width * 3;
	*width = cinfo.output_width;
	*height = cinfo.output_height;
	rgb_size = rgbrow_stride * cinfo.output_height; // 总大小 
	unsigned char *rgb_buffer = (unsigned char *)malloc(rgb_size + 54 + *height*3);
	memset(rgb_buffer, 0, sizeof(unsigned char)*rgb_size + 54 + *height * 3);
	tmp_buffer = rgb_buffer + rgb_size - rgbrow_stride;
	*outbuf = rgb_buffer;
	unsigned char * buffer[1];
	while (cinfo.output_scanline < cinfo.output_height) // 解压每一行  
	{
		buffer[0] = tmp_buffer;
		jpeg_read_scanlines(&cinfo, buffer, 1);
		tmp_buffer -= rgbrow_stride;
	}
	if (cinfo.output_components == 3)
	{
		for (int i = 0; i < (int)(cinfo.output_width*cinfo.output_height); i++)
		{
			unsigned char tmp;
			tmp = rgb_buffer[0];
			rgb_buffer[0] = rgb_buffer[2];
			rgb_buffer[2] = tmp;
			rgb_buffer += 3;
		}
	}
	else
	{
		unsigned char* buffer = rgb_buffer;
		for (int i = 0; i < (int)cinfo.output_height; i++)
		{
			buffer += row_stride;
			rgb_buffer += rgbrow_stride;
			for (int j = 0; j < (int)cinfo.output_width; j++)
			{
				rgb_buffer -= 3;
				buffer--;
				rgb_buffer[0] = buffer[0];
				rgb_buffer[1] = buffer[0];
				rgb_buffer[2] = buffer[0];
				
			}
			rgb_buffer += rgbrow_stride;
			buffer += rgbrow_stride;
		}
	}

	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	return true;
}