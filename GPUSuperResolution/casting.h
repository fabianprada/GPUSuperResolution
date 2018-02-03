
namespace CastingTools{
	unsigned char NormalFloat2Unsigned( float x);
	float * TransformUnsignedToNormalizedFloat(unsigned char * u_values, const unsigned int array_size, const unsigned int input_channels, const unsigned int output_channels);
	unsigned char * TransformNormalizedFloatToUnsigned(unsigned char * f_values, const unsigned int array_size, const unsigned int input_channels, const unsigned int output_channels);
}