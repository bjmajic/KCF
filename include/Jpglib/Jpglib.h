//JPGͼ���ʽ��
//��ȡjpgͼƬ�����bmpͼƬ
void jpg_to_bmp(const char* jpg_file, const char* bmp_file);
//��ȡbmpͼƬ�����jpgͼƬ
void bmp_to_jpg(const char* bmp_file, const char* jpg_file,int quality);
//���ڴ���ת�����ݸ�ʽ, outbuf ��Ҫ�ⲿ�����ڴ��ͷ�free. ��֧��BGR����ͼ���ʽ
void bmp2jpeg_compress(unsigned char *bgr, int width, int height, unsigned char **outbuf, unsigned long *outSize, int quality);
//���ڴ���ת�����ݸ�ʽ, outbuf ��Ҫ�ⲿ�����ڴ��ͷ�free. ��֧�����BGR����ͼ���ʽ
bool jpeg2bmp_decompress(unsigned char *jpeg, int jpeg_size, unsigned char **outbuf, int *width, int *height);
//���ڴ���ת�����ݸ�ʽ, outbuf ��Ҫ�ⲿ�����ڴ��ͷ�free. ��֧�����BGR����ͼ���ʽ
bool jpeg2bmp_decompress2(unsigned char *jpeg, int jpeg_size, unsigned char **outbuf, int *width, int *height);
//������BGRͼ��洢���ļ���
bool SaveRGB24JPG(void *img, int width, int height, const char* filename, int quality);