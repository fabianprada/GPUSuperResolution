#include "image.h"
#include "PNG.h"

void Image::WriteImagePNG(const char * file_name)
{
 	PNGWriteColor(file_name, colors, width,height);
}


Image * ReadImagePNG(const char * file_name) 
{
	int p_width, p_height;
	unsigned char * p_colors = PNGReadColor(file_name, p_width, p_height);
	Image * image = new Image(p_height, p_width,p_colors);
	sprintf(image->name, "%s", file_name);
	return image;
}


Image::Image(int p_height, int p_width, unsigned char * p_colors)
{
	width = p_width;
	height = p_height;
	colors = p_colors;
}

