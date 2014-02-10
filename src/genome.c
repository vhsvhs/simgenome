#include "common.h"


/* Build a genome, given a pre-built set of genes.
 * The input genes are copied (new) into the new genome.
 *
 * */
t_genome* make_genome(int ngenes, t_gene** ingenes, settings *ss) {
	//printf("\n. genome 9 - entered make_genome");
	fflush(stdout);
	t_genome *gn;
	gn = (t_genome*)malloc(1*sizeof(t_genome));
	//printf("\n. genome 11");
	gn->ngenes = ngenes;
	//printf("\n. genome 12");
	gn->genes = (t_gene**)malloc(ngenes*sizeof(t_gene));
	//printf("\n. genome 13");
	gn->ntfs = 0;
	//printf("\n. genome 14");
	gn->is_elite = false;
	if (ingenes != NULL){
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
}

/* Makes a genome with all random gene sequences and affinities */
t_genome* make_genome_random(settings* ss){
	t_gene** genes = (t_gene**)malloc(ss->ngenes*sizeof(t_gene));

	/* Build the genes, using the URSes and PSAMs we previously found. */
	for (int ii=0; ii < ss->ngenes; ii++){ // ii = gene
		int psamlen = 0;
		if (ii < ss->nreg){
			psamlen = PSAMLEN_DEFAULT;
		}

		genes[ii] = make_gene(psamlen, ss->urslen);
		genes[ii]->id = ii;
		for (int jj=0; jj < ss->urslen; jj++){ // jj = site
			genes[ii]->urs[jj] = get_random_state();
		}
		if (ii < ss->nreg){
			for (int jj=0; jj<N_STATES*psamlen; jj++){
				rand_init_psam( genes[ii]->dbd );
				genes[ii]->has_dbd = true;
				genes[ii]->reg_mode = ii%2;
			}
		}
		if (ss->verbosity > 10){
			printf("\n+ Gene %d, URS (%d sites)", ii, genes[ii]->urslen);
			if (genes[ii]->has_dbd == true){
				printf(", PSAM (%d sites)", genes[ii]->dbd->nsites);
			}
			printf("\n");
		}
	}

	/* Setup cofactor variables */
	for (int ii=0; ii<ss->nreg; ii++){
		build_coop(genes[ii], ss->nreg, ss->maxgd);
		init_coop(genes[ii]);
		calc_gamma(genes[ii], ss->nreg, ss->maxgd);
	}

	t_genome *gn = make_genome(ss->ngenes, genes, ss);
	return gn;
}

/*
 * Builds gene expres	sion array as a new object.
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
	for (int ii = 0; ii < g->ntfs; ii++){ //to-do
		if (g->genes[ii]->has_dbd){
			c += g->genes[ii]->dbd->nsites;
		}
	}
	return c;
}
