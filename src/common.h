/* The order of these include statements matters. */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "configuration.h"
#include "psam.h"
#include "ptable.h"
#include "settings.h"
#include "gene.h"

#include "genome.h"
#include "population.h"
#include "landscape.h"
#include "ga.h"
#include "fitness.h"
#include "mutate.h"
#include "argp.h"

#include "fout.h"
void shuffle(double *array, size_t n);
int filexists(char *filename);
char int2nt(int x);
int nt2int(char c);

double drand();
double random_normal();
int sample_from_cdf(double* p, int len);
int Filexists(char *filename);
