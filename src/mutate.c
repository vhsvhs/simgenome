#include "common.h"


/* Returns a random delta-delta-G value */
double get_random_ddg() {
	return (double)rand() / (double)RAND_MAX * (MAX_DDG - MIN_DDG) / 2.0;
}

/*
 * Mutates every individual in the population,
 * using the mutational settings defined in ss.
 */
void mutate(t_pop* pop, settings* ss){
	if (ss->verbosity > 2){
		printf("\n. The population is mutating. . .\n");
	}

	if (!ss->do_mutation){
		return;
	}

	for (int ii = 0; ii < pop->ngenomes; ii++) {
		/*
		 * n mutations to upstream regulatory sequence
		 */
		int n = count_urslen( pop->genomes[ii] ) * ss->urs_mu_rate;

		if (n < 0){
			n = 0;
		}
		if (ss->verbosity > 10){
			printf("\n. I'm making %d URS point mutations to ID %d.", n, ii );
		}
		for (int jj = 0; jj < n; jj++){
			int rand_gene = rand()%pop->genomes[ii]->ngenes;
			mutate_urs(pop->genomes[ii]->genes[rand_gene], ss);
		}

		/*
		 * n mutations to PSAMs
		 */
		n = count_psamlen( pop->genomes[ii] ) * ss->psam_mu_rate;
		if (n < 0){
			n = 0;
		}
		if (ss->verbosity > 10){
			printf("\n. I'm making %d PSAM point mutations to ID %d.", n, ii );
		}
		for (int jj = 0; jj < n; jj++){
			int rand_gene = rand()%pop->genomes[ii]->ntfs;
			mutate_psam(pop->genomes[ii]->genes[rand_gene]->dbd, ss);
		}


		/*
		 * to-do: mutate co-factor interactions
		 */

		// Rebuilt the co-factor affinity matrix based
		// on the new mutant affinity values.
		for (int jj=0; jj < pop->genomes[ii]->ntfs; jj++){
			calc_gamma( pop->genomes[ii]->genes[jj],
					pop->genomes[ii]->ntfs,
					ss->maxgd);
		}


		/*
		 * to-do: indels
		 */

	}

}

/* Mutates a single random site in the psam */
void mutate_psam(psam *p, settings *ss) {
	int rand_ii = (double)rand() / (double)RAND_MAX * (p->nsites*p->nstates);
	double old = p->data[rand_ii];
	p->data[rand_ii] = get_random_ddg();
	if (ss->verbosity > 20) {
		printf("\n. Mutating PSAM rand_ii %d old %f new %f\n", rand_ii, old, p->data[rand_ii]);
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
	int old = g->urs[ rand_site ];
	g->urs[ rand_site ] = rand_state;
	if (ss->verbosity > 20){
		printf("\n. Mutating URS rand_site %d rand_state %d old %d new %d\n",
				rand_site, rand_state, old, g->urs[ rand_site ]);
	}
}

void mutate_gamma(t_gene *g, settings *ss) {
	int rand_site = (int)( (float)rand()/(float)RAND_MAX * g->gammalen );
	double rand_value = ( (float)rand()/(float)RAND_MAX * 11 ) -1;
	g->gamma[ rand_site ] = rand_value;
}
