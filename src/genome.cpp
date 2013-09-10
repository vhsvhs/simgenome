#include "common.h"


/* Build a genome, given a pre-built set of genes. */
t_genome* make_genome(int ngenes, t_gene** ingenes, settings *ss) {
	t_genome *gn;
	gn = (t_genome*)malloc(1*sizeof(t_genome));
	gn->ngenes = ngenes;
	gn->genes = (t_gene**)malloc(ngenes*sizeof(t_gene));
	for (int ii=0; ii<ngenes; ii++){
		t_gene *this_gene;
		gn->genes[ii] = make_gene(ingenes[ii]->dbd->nsites, ingenes[ii]->urslen);
		copy_gene(gn->genes[ii], ingenes[ii]);
	}
	return gn;
}

/* Copies a genome, assuming that memory has already been allocated for the 'to' genome. */
void copy_genome(t_genome* to, t_genome* from){
	to->id = from->id;
	for (int ii=0; ii < to->ngenes; ii++){
		free_gene( to->genes[ii] );
	}
	free(to->genes);
	to->genes = (t_gene**)malloc(from->ngenes*sizeof(t_gene));
	for (int ii=0; ii< from->ngenes; ii++){
		to->genes[ii] = make_gene( from->genes[ii]->dbd->nsites,from->genes[ii]->urslen);
		copy_gene(to->genes[ii],from->genes[ii]);
	}
	to->ngenes = from->ngenes;
	to->is_elite = from->is_elite;
}

void free_genome(t_genome* gn){
	for(int ii=0; ii<gn->ngenes; ii++){
		free_gene( gn->genes[ii] );
		free(gn->genes[ii]);
	}
	free(gn->genes);
	free(gn->gene_expr);
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
			gn->genes[ii]->reg_mode = 1;
		}
		else {
			gn->genes[ii]->reg_mode = 0;
		}
	}
	return gn;
}

