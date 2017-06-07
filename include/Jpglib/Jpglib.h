//JPG图像格式库
//读取jpg图片，输出bmp图片
void jpg_to_bmp(const char* jpg_file, const char* bmp_file);
//读取bmp图片，输出jpg图片
void bmp_to_jpg(const char* bmp_file, const char* jpg_file,int quality);
//在内存中转换数据格式, outbuf 需要外部进行内存释放free. 仅支持BGR倒立图像格式
void bmp2jpeg_compress(unsigned char *bgr, int width, int height, unsigned char **outbuf, unsigned long *outSize, int quality);
//在内存中转换数据格式, outbuf 需要外部进行内存释放free. 仅支持输出BGR正立图像格式
bool jpeg2bmp_decompress(unsigned char *jpeg, int jpeg_size, unsigned char **outbuf, int *width, int *height);
//在内存中转换数据格式, outbuf 需要外部进行内存释放free. 仅支持输出BGR倒立图像格式
bool jpeg2bmp_decompress2(unsigned char *jpeg, int jpeg_size, unsigned char **outbuf, int *width, int *height);
//将倒立BGR图像存储到文件中
bool SaveRGB24JPG(void *img, int width, int height, const char* filename, int quality);