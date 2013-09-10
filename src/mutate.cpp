#include "common.h"


/* Returns a random delta-delta-G value */
double get_random_ddg() {
	return (double)rand() / (double)RAND_MAX * (MAX_DDG - MIN_DDG) / 2.0;
}

/* Mutates a single random site in the psam */
void mutate_psam(psam *p, settings *ss) {
	int rand_ii = (double)rand() / (double)RAND_MAX * (p->nsites*p->nstates);
	p->data[rand_ii] = get_random_ddg();\
	if (ss->verbosity > 100) {
		printf("(mutate 14) rand_ii = %d %f\n", rand_ii, p->data[rand_ii]);
	}
}

/* Inserts a point mutation into the URS of gene g */
void mutate_urs(t_gene *g, settings *ss) {
	int rand_site = (double)rand() / (double)RAND_MAX * (g->urslen);
	int rand_state = (int)( (double)rand() / (double)RAND_MAX * N_STATES );
	if (rand_state == g->urs[rand_site]) {
		rand_state += 1;
		rand_state = rand_state%N_STATES;
	}
	g->urs[ rand_site ] = rand_state;
}

void mutate_gamma(t_gene *g, settings *ss) {
	int rand_site = (int)( (float)rand()/(float)RAND_MAX * g->gammalen );
	double rand_value = ( (float)rand()/(float)RAND_MAX * 11 ) -1;
	g->gamma[ rand_site ] = rand_value;
}
