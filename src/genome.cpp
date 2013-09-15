#include "common.h"


/* Build a genome, given a pre-built set of genes.
 * The input genes are copied (new) into the new genome.
 *
 * */
t_genome* make_genome(int ngenes, t_gene** ingenes, settings *ss) {
	t_genome *gn;
	gn = (t_genome*)malloc(1*sizeof(t_genome));
	gn->ngenes = ngenes;
	gn->genes = (t_gene**)malloc(ngenes*sizeof(t_gene));
	gn->ntfs = 0;
	for (int ii=0; ii<ngenes; ii++){
		int psamlen = 0;
		if (ingenes[ii]->has_dbd == true){
			psamlen = ingenes[ii]->dbd->nsites;
		}
		gn->genes[ii] = copy_gene(ingenes[ii]);
		if (gn->genes[ii]->has_dbd == true){
			gn->ntfs++;
		}
	}
	return gn;
}



/* Copies a genome. Creates a new genome object from scratch.
 * The new genome will be identical to the original genome,
 * except that it's expression vector will be empty.
 */
t_genome* copy_genome(t_genome *org){
	t_genome *gn;
	gn = (t_genome*)malloc(1*sizeof(t_genome));
	gn->id = org->id;
	gn->ngenes = org->ngenes;
	gn->ntfs = org->ntfs;
	gn->is_elite = org->is_elite;
	gn->expr_timeslices = org->expr_timeslices;

	gn->genes = (t_gene**)malloc(org->ngenes*sizeof(t_gene));
	for (int ii=0; ii<org->ngenes; ii++){
		gn->genes[ii] = copy_gene(org->genes[ii]);
	}
	return gn;
}

void free_genome(t_genome* gn){
	for(int ii=0; ii<gn->ngenes; ii++){
		free_gene( gn->genes[ii] );
		free(gn->genes[ii]);
	}
	free(gn->genes);
	free(gn->gene_expr);
	free(gn->r);
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

/*
 * Builds gene expression array as a new object.
 *
 * Build timeslice-related genomic data.
 * Do this only after you've built the landscape
 * and, thus, know how many timeslices are required
 * per individual.
 */
void build_lifespan(t_genome* g, int t){
	g->gene_expr = (double *)malloc(g->ngenes*t*sizeof(double));
	/* gene_expr[gene->id * ntimeslices + timeslice] = value ranging from 0.0 to 1.0 */
	g->expr_timeslices = t;
	reset_lifespan(g);
}

/*
 * Reset the values in expr_timeslices to 0.0 */
void reset_lifespan(t_genome* g){
	for (int ii = 0; ii < g->ngenes*g->expr_timeslices; ii++){
		g->gene_expr[ii] = 0.0;
	}
}

/* How many nucleotides are in this genome?
 */
int count_urslen(t_genome* g){
	int c = 0;
	for (int ii = 0; ii < g->ngenes; ii++){
		c += g->genes[ii]->urslen;
	}
	return c;
}

/* What's the sum length of all PSAMs in this genome?
 */
int count_psamlen(t_genome* g){
	int c= 0;
	for (int ii = 0; ii < g->ngenes; ii++){
		if (g->genes[ii]->has_dbd){
			c += g->genes[ii]->dbd->nsites;
		}
	}
	return c;
}
