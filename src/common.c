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

/* This returns a random value sampled from the cosine-approximated
 * normal distribution.
 */
double random_normal() {  /* normal distribution, centered on 0, std dev 1 */
  return sqrt(-2*log(drand())) * cos(2*M_PI*drand());
}

double max( double* x, int len){
	double ret = 0.0;
	for (int ii = 0; ii < len; ii++){
		if (x[ii] > ret){ ret = x[ii]; }
	}
	return ret;
}
/* Returns the minimum fitness in an array of fitnesses, x */
double min( double* x, int len){
	double ret = 1.0;
	for (int ii = 0; ii < len; ii++){
		if (x[ii] < ret){ ret = x[ii]; }
	}
	return ret;
}
double mean( double* x, int len){
	double sum = 0.0;
	for (int ii = 0; ii < len; ii++){
		sum += x[ii];
	}
	return sum / len;
}
double stdev( double* x, int len){
	double m = mean(x, len);
	double sum = 1.0;
	for (int ii = 0; ii < len; ii++){
		int qq = (x[ii] - m);
		sum += qq*qq;
	}
	return sqrt(sum/len);
}
double sderr( double* x, int len){
	return stdev(x,len) / sqrt(len);
}

/* Given an array of doubles, and a length of that array, this method
 * returns a random index from that array, where the probability of each
 * index is weighted by its relative value.
 */
int sample_from_cdf(double* p, int len){
	double sum = 0.0;
	for (int ii = 0; ii < len; ii++){
		sum += p[ii];
	}

	double randp = drand() * sum;

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


/* This method comes from:
 * http://www.programmingsimplified.com/c/source-code/c-program-for-pattern-matching
 */
int match(char *a, char *b)
{
   int c;
   int position = 0;
   char *x, *y;

   x = a;
   y = b;

   while(*a)
   {
      while(*x==*y)
      {
         x++;
         y++;
         if(*x=='\0'||*y=='\0')
            break;
      }
      if(*y=='\0')
         break;

      a++;
      position++;
      x = a;
      y = b;
   }
   if(*a)
      return position;
   else
      return -1;
}

/* The comparison function for qsort */
double cmpfunc (double a, double b) {
   return ( a - b );
}

/* This qsort algorithm comes from http://www.ontko.com/pub/rayo/cs40/mergesort.c
 */
void quicksort(double *data, int begin, int end) {
  /***********************************************************************/
  /* Purpose: sorts a list of numbers in increasing order                */
  /* Input:                                                              */
  /*   data: an unsorted array of numbers                                */
  /*   begin: the beginning index of "data"                              */
  /*   end: the final index of "data"                                    */
  /* Output:                                                             */
  /*   upon termination of the subroutine, "data" contains the sorted    */
  /*     list of numbers                                                 */
  /***********************************************************************/

  int leftarrow, rightarrow;
  double temp, pivot;

  leftarrow = begin;
  rightarrow = end;
  pivot = data[(begin+end)/2];
  while (1)
  {
    while (data[rightarrow] > pivot)
      rightarrow = rightarrow - 1;

    while (data[leftarrow] < pivot)
      leftarrow = leftarrow + 1;

    if (leftarrow <= rightarrow)
    {
      temp = data[leftarrow];
      data[leftarrow] = data[rightarrow];
      data[rightarrow] = temp;
      leftarrow = leftarrow + 1;
      rightarrow = rightarrow - 1;
    }

    if (rightarrow < leftarrow)
      break;
  }

  if (begin < rightarrow)
    quicksort(data,begin,rightarrow);
  if (leftarrow < end)
    quicksort(data,leftarrow,end);

  return;
}


/*
 * Sample a random number from a normal distribution, with mean mu and stdev sigma.
 *
 * Code from here: http://phoxis.org/2013/05/04/generating-random-numbers-from-normal-distribution-in-c/
 */
double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}
