#include "common.h"


/* Returns a random delta-delta-G value */
double get_random_ddg() {
	return (double)rand() / (double)RAND_MAX * (MAX_DDG - MIN_DDG) / 2.0;
}

/* * Mutates a single random site in the psam */
void mutate_psam(psam *p, settings *ss) {
	int rand_ii = (double)rand() / (double)RAND_MAX * (p->nsites*p->nstates);
	p->data[rand_ii] = get_random_ddg();\
	if (ss->verbosity > 100) {
		printf("(mutate 14) rand_ii = %d %f\n", rand_ii, p->data[rand_ii]);
	}
}
