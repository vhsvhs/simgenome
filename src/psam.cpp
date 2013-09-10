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

void shuffle_sites(psam *p) {
	shuffle(p->data, p->nstates*p->nsites);
}

void rand_init(psam *p) {
	for( int ii = 0; ii < p->nstates*p->nsites; ii++) {
		p->data[ii] = (double)rand() / (double)RAND_MAX * (MAX_DDG - MIN_DDG) / 2.0;
		if (p->data[ii] < 0.5) { p->data[ii] *= -1; }
		//printf("\n(psam 21) %d\n", ii);
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
