#include "common.h"

t_pop* make_population_random(int popsize, int ngenes, settings *ss){
	t_pop *pop;
	if (ngenes == -1) {
		ngenes = NGENES_DEFAULT;
	}
	pop = (t_pop*)malloc(1*sizeof(t_pop));
	pop->genomes = (t_genome**)malloc(popsize*sizeof(t_genome));
	for(int ii=0; ii<popsize; ii++){
		t_genome* gn;
		gn = make_genome_random(ngenes);
		pop->genomes[ii] = gn;
	}
	return pop;
}

/* Builds a population that is 'popsize' copies of the genome 'gn'. */
t_pop* make_population(t_genome* gn, settings *ss){
	t_pop *pop;
	pop = (t_pop*)malloc(1*sizeof(t_pop));
	pop->genomes = (t_genome**)malloc(ss->popsize*sizeof(t_genome));
	pop->ngenomes = ss->popsize;
	/* Copy the gn into each individual */
	for(int ii=0; ii< ss->popsize; ii++){
		pop->genomes[ii] = (t_genome*)malloc(1*sizeof(t_genome));
		pop->genomes[ii] = copy_genome(gn);
		// The copy's ID won't be same as the input genome:
		pop->genomes[ii]->id = ii;
	}
	return pop;
}

/* Builds a population that is lacking genome objects */
t_pop* make_population_basic(t_genome* gn, settings *ss){
	t_pop *pop;
	pop = (t_pop*)malloc(1*sizeof(t_pop));
	pop->genomes = (t_genome**)malloc(ss->popsize*sizeof(t_genome));
	pop->ngenomes = ss->popsize;
	return pop;
}

void free_pop(t_pop* pop){
	for(int ii=0; ii < pop->ngenomes; ii++){
		free_genome(pop->genomes[ii]);
		free(pop->genomes[ii]);
	}
	free(pop->genomes);
}

/* f is an array of fitness scores */
t_pop* reproduce(t_pop* pop, settings* ss, double* f) {
	if (ss->verbosity > 2){
		printf("\n. The population is selectively reproducing. . .\n");
	}

	/* Open the mating log file */
	char *gs;
	gs = (char*)malloc(10*sizeof(char));
	sprintf(gs, "%d", ss->gen_counter);
	char* p = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat( strcat( strcat( strcat(p, ss->outdir), "/MATING/mating.gen"), gs), ".txt");
	free(gs);
	FILE *fp;
	fp = fopen(p, "w");
	if (fp == NULL) {
	  fprintf(stderr, "Error: can't open the mating log file %s!\n", p);
	}
	free(p);

	/* Build the population at t+1 */
	t_pop *newpop;
	newpop = (t_pop*)malloc(1*sizeof(t_pop));
	newpop->genomes = (t_genome**)malloc(pop->ngenomes*sizeof(t_genome));
	newpop->ngenomes = pop->ngenomes;

	for (int ii = 0; ii < pop->ngenomes; ii++){

		/* Elite individuals are cloned. */
		if (pop->genomes[ii]->is_elite == true){
			newpop->genomes[ii] = copy_genome( pop->genomes[ii]);
			if (ss->verbosity > 2){
				printf("\n\t. Elite individual %d will be cloned.", ii);
			}
			fprintf(fp, "generation %d cloned id %id into child %d\n",
						ss->gen_counter,
						ii, ii);
		}
		else{
			/* Pick random parents */
			int parent1 = sample_from_cdf(f, pop->ngenomes);
			int parent2 = sample_from_cdf(f, pop->ngenomes);
			/* Mate those parents */
			if (ss->verbosity > 2){
				printf("\n\t. ID %d mating with ID %d", parent1, parent2);
			}

			fprintf(fp, "generation %d id %d x id %d = child %d\n",
						ss->gen_counter, parent1, parent2, ii);
			newpop->genomes[ii] = mate(pop->genomes[parent1],
					pop->genomes[parent2]);
		}
		newpop->genomes[ii]->id = ii;
		build_lifespan(newpop->genomes[ii], ss->maxtime);
	}

	/* Free the old population */
	free_pop(pop);
	fclose(fp);
	return newpop;
}


/* Mates two genomes, and returns their child as a new object.*/
t_genome* mate(t_genome* par1, t_genome* par2){
	t_genome *f1; /* The F1 cross of parent par1 with parent par2 */
	f1 = (t_genome*)malloc(1*sizeof(t_genome));
	f1->genes = (t_gene**)malloc(par1->ngenes*sizeof(t_gene));

	/* For each gene */
	for (int gg = 0; gg < par1->ngenes; gg++){
		/* Flip a coin for the child to get par1's copy of the gene
		 * versus par2's copy.
		 */
		double flip = (float)rand() / (float)RAND_MAX;
		int psamlen = 0;
		int urslen = 0;
		if (flip > 0.5){
			if (par1->genes[gg]->has_dbd){
				psamlen = par1->genes[gg]->dbd->nsites;
				urslen = par1->genes[gg]->urslen;
			}
			f1->genes[gg] = make_gene(psamlen, urslen);
			f1->genes[gg] = copy_gene(par1->genes[gg]);
		}
		else{
			if (par2->genes[gg]->has_dbd){
				psamlen = par2->genes[gg]->dbd->nsites;
				urslen = par2->genes[gg]->urslen;
			}
			f1->genes[gg] = make_gene(psamlen, urslen);
			f1->genes[gg] = copy_gene(par2->genes[gg]);
		}
	}

	f1->ngenes = par1->ngenes;
	f1->ntfs = par1->ntfs;
	f1->is_elite = false;

	return f1;
}

void print_population(t_pop* pop, settings* ss){
	for(int ii=0; ii< pop->ngenomes; ii++){
		printf("\n Individual %d:\n", ii);
		for(int jj=0; jj < pop->genomes[ii]->ngenes; jj++){
			printf("\n Gene %d:\n", jj);
			print_urs( pop->genomes[ii]->genes[jj]->urs, pop->genomes[ii]->genes[jj]->urslen);
			print_psam( pop->genomes[ii]->genes[jj]->dbd );
		}
	}
}

