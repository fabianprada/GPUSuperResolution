#include "extension.h"

int ExtensionTools::ExtensionIndex(int index, int n, ExtensionMode mode)
{
	if (mode == MIRROR_EXTENSION){
		if (index < 0){
			index = -1 - index;
		}
			int ratio = index / n;
			int residue = index - ratio*n;
			if (ratio % 2 == 0)
				return residue;
			else
				return n - 1 - residue;
	}
	else if(mode == REPEAT_EXTENSION){
		if (index < 0){
			index = n + index;
			return ExtensionIndex(index,n,REPEAT_EXTENSION);
		}
			return (index % n);
	}
	else{
		printf("Unimplemented!!\n");
		return 0;
	}
}

double ExtensionTools::ExtensionValue(double input, double x, ExtensionMode mode)
{
	if (mode == MIRROR_EXTENSION){
		if (input < 0){
			input *= -1.f;
		}
		double ratio = floor(input / x);
		double residue = input - ratio*x;
		if (static_cast<int>(ratio) % 2 == 0)
			return residue;
		else
			return x - residue;
	}
	else{
		printf("Unimplemented!!\n");
		return 0.f;
	}
}



