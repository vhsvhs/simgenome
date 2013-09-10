/* The order of these include statements matters. */
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <time.h>

#include "configuration.h"
#include "psam.h"
#include "ptable.h"

#include "gene.h"
#include "settings.h"
#include "genome.h"
#include "population.h"

#include "mutate.h"
#include "fout.h"

void shuffle(double *array, size_t n);
int filexists(char *filename);
char int2nt(int x);
int nt2int(char c);
