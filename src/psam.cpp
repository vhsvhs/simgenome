#include "common.h"

/* Initializes memory for a PSAM */
psam* make_psam(int nsites, int nstates) {
	//printf("(psam 5)\n");
	psam *p = (psam *)malloc(1*sizeof(psam));
	p->nstates = nstates;
	p->nsites = nsites;
	//printf("(psam 8)\n");
	p->data = (double*)malloc(nstates*nsites*sizeof(double));
	return p;
}

void shuffle_psam(psam *p) {
	shuffle(p->data, p->nstates*p->nsites);
}

/* Assigns a random affinity for all states and all sites */
void rand_init_psam(psam *p) {
	for( int ii = 0; ii < p->nstates*p->nsites; ii++) {
		p->data[ii] = get_random_ddg();
		if (p->data[ii] < 0.5) { p->data[ii] *= -1; }
	}
}

void print_psam(psam *p) {
	for (int ii=0; ii < p->nsites; ii++){ // ii = site
		for (int jj=0; jj < p->nstates; jj++){ // jj = state
			printf("\n. site %d state %d ddG %f", ii, jj, p->data[ii*p->nstates + jj]);
		}
	}
	printf("\n.");
}

/* Returns the summed delta-delta-G affinity for the entire sequence. */
double get_affinity(psam *p, int *seq, int len) {
	double sum = 0.0;
	for (int ii=0; ii < len; ii++) {
		sum += p->data[ ii*p->nstates + seq[ii] ];
	}
	return sum;
}

