#include "common.h"

t_genome* make_genome_basic(int ngenes){
//	t_genome *gn;
//	gn = (t_genome *)malloc(1*sizeof(t_genome));
//	gn->is_elite = false;
//
//	gn->genes = (t_gene*)malloc(ngenes*sizeof(t_gene));
//	for (int ii=0; ii<ngenes; ii++){
//		gn->genes[ii] = make_gene();
//		gn->genes[ii]->id = ii;
//	}
}

t_genome* make_genome_default(int ngenes, settings *ss){
	t_genome *gn;
	gn = (t_genome*)malloc(1*sizeof(t_genome));
	gn->ngenes = ngenes;
	gn->genes = (t_gene**)malloc(ngenes*sizeof(t_gene));
	for (int ii=0; ii<ngenes; ii++){
		t_gene *this_gene;
		this_gene = make_gene(PSAMLEN_DEFAULT, URSLEN_DEFAULT);
		gn->genes[ii] = (t_gene*)this_gene;
		gn->genes[ii]->id = ii;
	}
	return gn;
}

/* Makes a genome with all random gene sequences and affinities */
t_genome* make_genome_random(int ngenes){
	t_genome *gn;
	gn = (t_genome *)malloc(1*sizeof(t_genome));
	gn->is_elite = false;
	gn->ngenes = ngenes;
	gn->genes = (t_gene**)malloc(ngenes*sizeof(t_gene));
	for (int ii=0; ii<ngenes; ii++){
		t_gene *this_gene;
		this_gene = make_gene(PSAMLEN_DEFAULT, URSLEN_DEFAULT);
		gn->genes[ii] = this_gene;
		gn->genes[ii]->id = ii;
		gn->genes[ii]->urs = get_random_seq(URSLEN_DEFAULT);
		gn->genes[ii]->urslen = URSLEN_DEFAULT;

		double rand_f = (float)rand()/(float)RAND_MAX;
		if (rand_f > 0.5){
			gn->genes[ii]->has_dbd = true;
			rand_init_psam( gn->genes[ii]->dbd );
		}
		else {
			gn->genes[ii]->has_dbd = false;
		}

		rand_f = (float)rand()/(float)RAND_MAX;
		if (rand_f > 0.5){
			gn->genes[ii]->is_repressor = true;
		}
		else {
			gn->genes[ii]->is_repressor = false;
		}
	}
	return gn;
}

