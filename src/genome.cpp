#include "common.h"


/* Build a genome, given a pre-built set of genes.
 * The input genes are copied (new) into the new genome. */
t_genome* make_genome(int ngenes, t_gene** ingenes, settings *ss) {
	//printf("genome 7\n");
	t_genome *gn;
	gn = (t_genome*)malloc(1*sizeof(t_genome));
	gn->ngenes = ngenes;
	gn->genes = (t_gene**)malloc(ngenes*sizeof(t_gene));
	gn->ntfs = 0;
	//printf("genome 13\n");
	for (int ii=0; ii<ngenes; ii++){
		//printf("genome 15 gene %d\n", ii);
		int psamlen = 0;
		if (ingenes[ii]->has_dbd == true){
			psamlen = ingenes[ii]->dbd->nsites;
		}
		gn->genes[ii] = make_gene(psamlen, ingenes[ii]->urslen);
		//printf("genome 17\n");
		copy_gene(gn->genes[ii], ingenes[ii]);
		if (gn->genes[ii]->has_dbd == true){
			gn->ntfs++;
		}
	}
	//printf("genome 20\n");
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
	to->ntfs = from->ntfs;
}

/* An alternative to copy_genome.  This version creates a new
 * genome object from scratch.
 */
t_genome* dup_genome(t_genome *org){
	t_genome *gn;
	gn = (t_genome*)malloc(1*sizeof(t_genome));
	gn->ngenes = org->ngenes;
	gn->genes = (t_gene**)malloc(org->ngenes*sizeof(t_gene));
	gn->ntfs = org->ntfs;
	gn->is_elite = org->is_elite;
	gn->expr_timeslices = org->expr_timeslices;
	for (int ii=0; ii<org->ngenes; ii++){
		gn->genes[ii] = make_gene(org->genes[ii]->dbd->nsites, org->genes[ii]->urslen);
		copy_gene(gn->genes[ii], org->genes[ii]);
	}
	init_lifespan(gn, org->expr_timeslices);
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

t_genome* make_genome_default(int ngenes, settings *ss){
	t_genome *gn;
	gn = (t_genome*)malloc(1*sizeof(t_genome));
	gn->ngenes = ngenes;
	gn->genes = (t_gene**)malloc(ngenes*sizeof(t_gene));
	gn->ntfs = 0;
	for (int ii=0; ii<ngenes; ii++){
		t_gene *this_gene;
		this_gene = make_gene(PSAMLEN_DEFAULT, URSLEN_DEFAULT);
		gn->genes[ii] = (t_gene*)this_gene;
		gn->genes[ii]->id = ii;
		if (gn->genes[ii]->has_dbd){
			gn->ntfs++;
		}
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

/*
 * Builds gene expression array as a new object.
 *
 * Build timeslice-related genomic data.
 * Do this only after you've built the landscape
 * and, thus, know how many timeslices are required
 * per individual.
 */
void init_lifespan(t_genome* g, int t){
	g->gene_expr = (double *)malloc(g->ngenes*t*sizeof(double));
	/* gene_expr[gene->id * ntimeslices + timeslice] = value ranging from 0.0 to 1.0 */
	g->expr_timeslices = t;
}

/*
 * Builds gene expression array IN PLACE, rather than new.
 *
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
