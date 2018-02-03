#ifndef IMAGE_INCLUDED
#define IMAGE_INCLUDED

class Image{
public:

	Image(int p_height, int p_width)
	{
		height = p_height;
		width = p_width;
		colors = new unsigned char[height*width*4];
	}
	Image(int p_height, int p_width, unsigned char * u_colors);
	~Image(){ delete colors; }

	//unsigned char * GenerateUnsignedCharPixels(const int channels);
	//void WriteImagePPM(const char * file_name);
	void WriteImagePNG(const char * file_name);

	unsigned char * colors;
	int height;
	int width;
	char name[100];
};

Image * ReadImagePNG(const char * file_name) ;

#endif// IMAGE_INCLUDED