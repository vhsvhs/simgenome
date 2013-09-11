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
t_pop* make_population(int popsize, t_genome* gn, settings *ss){
	t_pop *pop;
	pop = (t_pop*)malloc(1*sizeof(t_pop));
	pop->genomes = (t_genome**)malloc(popsize*sizeof(t_genome));
	pop->ngenomes = popsize;
	for(int ii=0; ii<popsize; ii++){
		pop->genomes[ii] = make_genome(gn->ngenes, gn->genes, ss);
		pop->genomes[ii]->id = ii;
	}
	return pop;
}

void free_pop(t_pop* pop){
	for(int ii=0; ii < pop->ngenomes; ii++){
		free_genome(pop->genomes[ii]);
		free(pop->genomes[ii]);
	}
	free(pop->genomes);
}

void reproduce(t_pop* pop, settings* ss) {
	int x;
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
