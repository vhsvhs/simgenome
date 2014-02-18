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

void free_psam(psam* p){
	free(p->data);
	free(p);
}

/* Copies a psam, assuming the memory has already been allocated. */
void copy_psam(psam *to, psam *from ) {
	/* If the PSAMs have different dimensions, then re-allocate memory */
	if (to->nstates*to->nsites != from->nstates*from->nsites){
		free(to->data);
		to->data = (double*)malloc(from->nstates*from->nsites*sizeof(double));
	}
	to->nstates = from->nstates;
	to->nsites = from->nsites;
	for(int ii=0; ii<from->nstates*from->nsites; ii++){
		to->data[ii] = from->data[ii];
	}
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
	for (int ii=0; ii < p->nsites; ii++){ // sites
		for (int jj=0; jj < p->nstates; jj++){ //states
			char state = int2nt( jj );
			double value = p->data[ii*p->nstates + jj];
			printf("%c (%f)\t", state, value );
		}
		printf("\n");
	}
}

/* Returns the summed delta-delta-G affinity for the entire sequence. */
double get_affinity(psam *p, int *seq, int seqlen, int startsite) {
	//printf("psam 57 %d %d\n", seqlen, startsite);
	//print_psam(p);

	double sum = 0.0;
	if (seqlen-startsite > p->nsites){
		seqlen = startsite + p->nsites;
	}
	for (int ii=startsite; ii < seqlen; ii++) {
		sum += p->data[ (ii-startsite)*p->nstates + seq[ii] ];
	}

	if (sum < 0){
		sum = 0.0;
	}

	return sum;
}

