#include "common.h"

/* This main method is a sandbox for testing the objects
 * and methods of simgenome-c.
 */
int main( int argc, const char* argv[] )
{
	srand(time(0));

	// Prints each argument on the command line.
	for( int i = 0; i < argc; i++ )
	{	printf( "(scratch 8) arg %d: %s\n", i, argv[i] );
	}

	settings *ss = make_settings();
	//char urspath[] = "/Users/victor/Applications/simgenome-c-beta/examples/five.urs";
	char urspath[] = "/Users/victor/Applications/simgenome-c-beta/examples/urs.time092013.txt";
	ss->urspath = urspath;
	char psampath[] = "/Users/victor/Applications/simgenome-c-beta/examples/five.psam";
	ss->psampath = psampath;

	//psam *p = make_psam(6,4);
	//printf("(line 12) %d %d\n", p->nsites, p->nstates);
	//rand_init_psam(p);
	//print_psam(p);
	//shuffle_psam(p);
	//print_psam(p);

	ss->niid = 100;

	t_gene** mygenes;
	int ngenes;
	mygenes = read_genes_from_file(ss, ngenes);

	printf("(scratch 30) I found %d genes\n", ngenes);
//	for(int ii=0; ii<ngenes; ii++){
//		printf("\n Gene %d:\n", ii);
//		print_urs( mygenes[ii]->urs, mygenes[ii]->urslen);
//		print_psam( mygenes[ii]->dbd );
//	}

	t_genome *gn = make_genome(ngenes, mygenes, ss);
	printf("\n\n(scratch 39) I built a genome from the genes.\n", ngenes);
//	for(int ii=0; ii < gn->ngenes; ii++){
//		printf("\n Gene %d:\n", ii);
//		print_urs( gn->genes[ii]->urs, gn->genes[ii]->urslen);
//		print_psam( gn->genes[ii]->dbd );
//	}

	printf("\n. I'm building a population....\n");
	t_pop *pop = make_population(1, gn, ss);
	//printf("(scratch 48) %d\n", pop->ngenomes);
	//print_population(pop, ss);

	printf("\n. I'm building a landscape....\n");
	t_landscape *l = (t_landscape *)malloc(sizeof(t_landscape));
	char rulepath[] = "/Users/victor/Applications/simgenome-c-beta/examples/five.rules.txt";
	ss->rulepath = rulepath;
	t_ruleset** rs = read_rulesets_from_file(ss, l->nrulesets, l->ntime);
	l->rulesets = rs;
	printf("\n. I found %d rulesets.\n", l->nrulesets);
	for (int ii=0; ii < l->nrulesets; ii++) {
		printf("\n. ruleset %d:\n", ii);
		for (int jj=0; jj < l->rulesets[ii]->ninputs; jj++){
			printf("\n. input %d start %d stop %d gid %d expr %f\n", jj, l->rulesets[ii]->inputs[jj]->start, l->rulesets[ii]->inputs[jj]->stop, l->rulesets[ii]->inputs[jj]->gid, l->rulesets[ii]->inputs[jj]->expr_level);
		}
		for (int jj=0; jj < l->rulesets[ii]->nrules; jj++){
			printf("\n. rule %d timepoint %d repid %d expr %f type %d weight %f\n", jj, l->rulesets[ii]->rules[jj]->timepoint, l->rulesets[ii]->rules[jj]->repid, l->rulesets[ii]->rules[jj]->expr_level, l->rulesets[ii]->rules[jj]->rule_type, l->rulesets[ii]->rules[jj]->weight);
		}
	}

	printf("\n. I'm building a G.A.\n");
	t_ga* ga = make_ga();
	ga->pop = pop;
	ga->l = l;

	printf("\n. I'm running the G.A.\n");
	runsim(ga, ss);

	//double f = get_fitness(pop->genomes[0], l, ss);


	printf("\n. I'm freeing the population....\n");
	free_pop( pop );

	/* This next block should raise an exception
	 * because the population was garbage collected.
	 */
//	for(int ii=0; ii < pop->genomes[0]->ngenes; ii++){
//		printf("\n Gene %d:\n", ii);
//		print_urs( pop->genomes[0]->genes[ii]->urs, pop->genomes[0]->genes[ii]->urslen);
//		print_psam( pop->genomes[0]->genes[ii]->dbd );
//	}

	//printf("(scratch 23) %f %d %f\n", ss->pwmlenmu, ss->pwmlenmumax, ss->ddgmu);

//	ss->verbosity = 101;
//	mutate_psam(p, ss);
//	printf("After mutate_psam:\n");
//	print_psam(p);
//
//	int seqlen = 6;
//	int *seq;
//	seq = (int *)malloc(seqlen*sizeof(int));
//	for (int ii=0; ii<seqlen; ii++) {
//		seq[ii] = (int)((float)rand()/(float)RAND_MAX);
//	}
//	double aff = get_affinity(p, seq, 6);
//	printf("(scratch 30) aff = %f\n", aff);
//	free(seq);
//
//	t_gene *g;
//	g = make_gene(6, 1000);
//	rand_init_psam( g->dbd );
//	g->urs = get_random_seq(1000);
//
//	printf("Random Gene:\n");
//	print_psam( g->dbd );
//	print_urs( g->urs, g->urslen );
//	aff = get_affinity(g->dbd, seq, 6);
//	printf("(scratch 44) aff = %f\n", aff);

//	printf("\n. Random Genome:\n");
//	t_genome *gn;
//	gn = make_genome_random(NGENES_DEFAULT);
//	for (int ii=0; ii<NGENES_DEFAULT; ii++){
//		printf("\n. Gene %d\n", ii);
//		if (gn->genes[ii]->has_dbd) {
//			print_psam( gn->genes[ii]->dbd );
//			if (gn->genes[ii]->is_repressor){
//				printf(" -> is a repressor\n");
//			}
//		}
//		//print_urs( g->urs, g->urslen );
//	}

//	t_pop *pop;
//	pop = make_population(4, -1, ss);
//	char outpath[] = "./psams.txt";
//	write_psams(pop->genomes[0], ss, outpath);
//
//	t_ptable *pt;
//	pt = make_ptable(4, 3, 1000);
//	print_ptable(pt);

}
