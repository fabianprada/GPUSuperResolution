namespace FileIO{

void readInputdd(char * input, double *A, int num_entries);
void readInputff(char * input, float *A, int num_entries);
void readInputdf(char * input, float *A, int num_entries);
void readInputii(char * input, int *A, int num_entries);

void writeOutputdd(char * output, double *A, int num_entries);
void writeOutputff(char * output, float *A, int num_entries);
void writeOutputii(char * output, int *A, int num_entries);
void writeOutputui(char * output, unsigned char *A, int num_entries);

}