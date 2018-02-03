#include "casting.h"
#include <math.h>
#include <stdio.h>

unsigned char CastingTools::NormalFloat2Unsigned( float x)
{
	int i_color = static_cast<int>(floor(x*255.f));
	if (i_color < 0)
		i_color = 0;
	else if (i_color>255)
		i_color = 255;
	return static_cast<unsigned char>(i_color);
}

float * CastingTools::TransformUnsignedToNormalizedFloat(unsigned char * u_values, const unsigned int array_size, const unsigned int input_channels, const unsigned int output_channels)
{
	float * f_values =new float[array_size*output_channels];
	for(int i=0; i<array_size; i++){
				if(output_channels == input_channels){
					for(int k =0; k < input_channels; k++){
						f_values[i*input_channels + k] = ((float)u_values[i*input_channels + k])/255.f;
					}
				}
				else if(output_channels == 1 &&  input_channels == 3){
					f_values[i] = (((float)u_values[i*3])*0.3f + ((float)u_values[i*3 + 1])*0.59f + ((float)u_values[i*3 + 2])*0.11f) /255.f;
				}
				else if(output_channels == 1 &&  input_channels == 4){
					f_values[i] = (((float)u_values[i*4])*0.3f + ((float)u_values[i*4 + 1])*0.59f + ((float)u_values[i*4 + 2])*0.11f) /255.f;
				}
				else{
					printf("Undefined Dimension Transformation!!\n");
				}
	}
	return f_values;
}

unsigned char * CastingTools::TransformNormalizedFloatToUnsigned(unsigned char * f_values, const unsigned int array_size, const unsigned int input_channels, const unsigned int output_channels)
{
	unsigned char * u_values = new unsigned char[array_size*output_channels];
	for(int i=0; i<array_size;i++){

		if(output_channels == input_channels){
			for(int k =0; k < input_channels; k++){
			u_values[i*input_channels + k] = NormalFloat2Unsigned(f_values[i*input_channels + k]);
			}
		}
		else if(output_channels == 1 &&  input_channels == 3){
			u_values[i] = NormalFloat2Unsigned((f_values[i*3])*0.3f + (f_values[i*3 + 1])*0.59f + (f_values[i*3 + 2])*0.11f);
		}
		else if(output_channels == 1 &&  input_channels == 4){
			u_values[i] = NormalFloat2Unsigned((f_values[i*4])*0.3f + (f_values[i*4 + 1])*0.59f + (f_values[i*4 + 2])*0.11f);
		}
		else{
			printf("Undefined Dimension Transformation!!\n");
		}
	}
	return u_values;
}