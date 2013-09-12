#include "common.h"

/* Arrange the N elements of ARRAY in random order.
   Only effective if N is much smaller than RAND_MAX;
   if this may not be the case, use a better random
   number generator. */
void shuffle(double *array, size_t n)
{
    if (n > 1)
    {
        size_t i;
        for (i = 0; i < n - 1; i++)
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          double t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}

int filexists(char *filename)
{
  FILE *fp;
  fp =fopen(filename,"r");
  if (fp) {
    fclose(fp);
    return 1;
  } else
    return 0;
}

char int2nt(int x) {
	char state = '-';
	if (x == 0)   {state='A';}
	else if (x==1){state='C';}
	else if (x==2){state='G';}
	else if (x==3){state='T';}
	return state;
}

int nt2int(char c) {
	int i = '-';
	if (c == 'A')   {i=0;}
	else if (c=='C'){i=1;}
	else if (c=='G'){i=2;}
	else if (c=='T'){i=3;}
	return i;
}


double drand() {   /* uniform distribution, (0..1] */
  return (rand()+1.0)/(RAND_MAX+1.0);
}

double random_normal() {  /* normal distribution, centered on 0, std dev 1 */
  return sqrt(-2*log(drand())) * cos(2*M_PI*drand());
}

int sample_from_cdf(double* p, int len){
	double sum = 0.0;
	for (int ii = 0; ii < len; ii++){
		sum += p[ii];
	}

	double randp = drand();
	double c = 0.0;
	for (int ii = 0; ii < len; ii++){
		c += p[ii];
		if (c > randp){
			return ii;
		}
	}
	return len-1;
}

int Filexists(char *filename)
{
  FILE *fp;
  fp =fopen(filename,"r");
  if (fp) {
    fclose(fp);
    return 1;
  } else
    return 0;
}
