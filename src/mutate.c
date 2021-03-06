#include "common.h"


/* Returns a random delta-delta-G value */
double get_random_ddg() {
	double rand_ii = drand();
	if (rand_ii > 0.5){
		return (double)rand() / (double)RAND_MAX * (MAX_DDG) / 2.0;
	}
	else{
		return (double)rand() / (double)RAND_MAX * (MIN_DDG) / 2.0;
	}
}

/*
 * Mutates every individual in the population,
 * using the mutational settings defined in ss.
 */
void mutate(t_pop* pop, settings* ss){

	if (ss->do_mutation == false){
		if (ss->verbosity > 5){
			printf("\n\n. Mutation is disabled, skipping it.");
		}
		return;
	}
	if (ss->elite_proportion == 1.0){
		if (ss->verbosity > 3){
			printf("\n\n. All individuals in the population are 'elite',\n");
			printf("  and therefore no mutations will be made.\n");
			printf("  You can enable mutations by setting the value of\n");
			printf("  the parameter --elite_prop to be less than 1.\n");
		}
		return;
	}

	FILE *fp;
	if (ss->verbosity > 2){
		printf("\n\n. The population is mutating. . .\n");

		/* Open the mating log file */
		char *gs;
		gs = (char*)malloc(10*sizeof(char));
		sprintf(gs, "%d", ss->gen_counter);
		char* p = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
		strcat( strcat( strcat( strcat(p, ss->outdir), "/MUTATIONS/mu.gen"), gs), ".txt");
		free(gs);

		fp = fopen(p, "w");
		if (fp == NULL) {
		  fprintf(stderr, "Error: can't open the mating log file %s!\n", p);
		}
		free(p);
	}

	for (int ii = 0; ii < pop->ngenomes; ii++) {

		/* No mutations to elite individuals. */
		if (pop->genomes[ii]->is_elite){
			continue;
		}


		/*
		 * n mutations to upstream regulatory sequence
		 */
		int n = count_urslen( pop->genomes[ii] ) * randn( ss->urs_mu_rate, ss->mu_stdev ); //randn draws from a normal distribution

		if (n < 0){
			n = 0;
		}
		for (int jj = 0; jj < n; jj++){ // for gene jj
			int rand_gene = rand()%pop->genomes[ii]->ngenes;

			if (ss->mu_reporters_only && rand_gene < pop->genomes[ii]->ntfs){ continue; } // skip TF mutations
			if (ss->mu_regulators_only && rand_gene >= pop->genomes[ii]->ntfs){ continue; } // skip TF mutations

			int rs = mutate_urs(pop->genomes[ii]->genes[rand_gene], ss);

			if (ss->verbosity > 2){
				fprintf(fp, "SNP: ID %d URS %d site %d\n", ii, rand_gene, rs);
			}
		}
		if (ss->verbosity > 3){
			printf("\n. +%d SNPs to ID %d.", n, ii );
		}

		//
		// to-do: mutate URS lengths
		//

		/*
		 * n mutations to PSAMs
		 */
		n = count_psamlen( pop->genomes[ii] ) * randn(ss->psam_mu_rate, ss->mu_stdev); // randn draws from a normal distribution
		if (n < 0){
			n = 0;
		}
		if (ss->mu_reporters_only){ n=0; } // if we're only mutating reporter genes, then skip this next part.
		for (int jj = 0; jj < n; jj++){
			int rand_gene = rand()%pop->genomes[ii]->ntfs;

			int rs = mutate_psam(pop->genomes[ii]->genes[rand_gene]->dbd, ss);

			if (ss->verbosity > 2){
				fprintf(fp, "PSAM mutant: ID %d gene %d index %d\n", ii, rand_gene, rs);
			}
		}
		if (ss->verbosity > 3){
			printf("\n. I'm making %d PSAM point mutations to ID %d.", n, ii );
		}


		/*
		 * to-do: mutate co-factor interactions
		 */


		/*
		 * Mutate PSAM lengths
		 * */
		n = count_psamlen( pop->genomes[ii] ) * randn(ss->psamlenmu, ss->psamlensd);
		if (ss->mu_reporters_only){ n=0; }
		for (int jj = 0; jj < n; jj++){
			int rand_gene = rand()%pop->genomes[ii]->ntfs;
			int before_len = pop->genomes[ii]->genes[rand_gene]->dbd->nsites;
			mutate_psamlength(pop->genomes[ii]->genes[rand_gene]->dbd, ss);
			if (ss->verbosity > 3){
				fprintf(fp, "PSAM indel: ID %d gene %d old_len = %d, new_len = %d\n",
						ii, rand_gene, before_len, pop->genomes[ii]->genes[rand_gene]->dbd->nsites);
			}
		}
		if (ss->verbosity > 3){
			printf("\n. I'm making %d PSAM indels to ID %d.", n, ii );
		}

		/* Finally, rebuild the co-factor affinity matrix based
		 * on the new mutant affinity values.
		 */
		for (int jj=0; jj < pop->genomes[ii]->ntfs; jj++){
			calc_gamma( pop->genomes[ii]->genes[jj],
					pop->genomes[ii]->ntfs,
					ss->maxgd);
		}
	}

	if (ss->verbosity > 2){
		fclose(fp);
	}
}

/* Mutates a single random site in the psam.
 * Returns the index in psam->data that was mutated. */
int mutate_psam(psam *p, settings *ss) {
	int rand_ii = (double)rand() / (double)RAND_MAX * (p->nsites*p->nstates);
	double old = p->data[rand_ii];
	p->data[rand_ii] = get_random_ddg();
	if (ss->verbosity > 20) {
		printf("\n. Mutating PSAM rand_ii %d old %f new %f\n", rand_ii, old, p->data[rand_ii]);
	}
	return rand_ii;
}

void mutate_psamlength(psam *p, settings *ss){
	/* Insertion or deletion? */
	bool insert = true;
	if (p->nsites > 1){
		double rand_ii = drand();
		if (rand_ii > 0.5){
			insert = false;
		}
	}

	int erg = p->nstates * p->nsites;
	int size = 0.0;
	if (insert){
		size =  p->nstates*(p->nsites+1);
	}
	else{ //delete
		size = p->nstates*(p->nsites-1);
	}
	double* newdata = (double *)malloc( size * sizeof(double));

	if (insert){
		for (int ii = 0; ii < erg; ii++){
			newdata[ii] = p->data[ii];
		}
		for (int ii = erg; ii < erg + p->nstates; ii++){
			newdata[ii] = get_random_ddg();
		}
	}
	else{ // delete
		for (int ii = 0; ii < size; ii++){
			newdata[ii] = p->data[ii];
		}
	}
	free(p->data);
	p->data =(double *)malloc( size * sizeof(double));
	for (int ii = 0; ii < size; ii++){
		p->data[ii] = newdata[ii];
	}

	if (insert){ p->nsites++; }
	else{ p->nsites--; }
}

/* Inserts a point mutation into the URS of gene g.
 * Returns the site number of the URS that was mutated.
 * */
int mutate_urs(t_gene *g, settings *ss) {
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
	return rand_site;
}

void mutate_gamma(t_gene *g, settings *ss) {
	int rand_site = (int)( (float)rand()/(float)RAND_MAX * g->gammalen );
	double rand_value = ( (float)rand()/(float)RAND_MAX * 11 ) -1;
	g->gamma[ rand_site ] = rand_value;
}
