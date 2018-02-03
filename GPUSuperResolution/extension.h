#ifndef EXTENSION_INCLUDED
#define EXTENSION_INCLUDED

#include<iostream>

enum ExtensionMode { MIRROR_EXTENSION, ZERO_EXTENSION, REPEAT_EXTENSION, NEAREST_EXTENSION };
namespace ExtensionTools{
	int ExtensionIndex(int index, int n, ExtensionMode mode = MIRROR_EXTENSION );
	double ExtensionValue(double input, double x, ExtensionMode mode);
}

#endif //EXTENSION_INCLUDED