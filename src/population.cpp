#include "common.h"

t_pop* make_population(int popsize, int ngenes, settings *ss){
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
