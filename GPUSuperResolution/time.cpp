#include <stdlib.h>
#include <sys/timeb.h>
#ifndef WIN32
#include <sys/time.h>
#endif

double GetTime(void){
#ifdef WIN32
	struct _timeb t;
	_ftime(&t);
	return t.time+t.millitm/1000.0;
#else
	struct timeval t;
	gettimeofday(&t,NULL);
	return t.tv_sec+(double)t.tv_usec/1000000;
#endif
}
