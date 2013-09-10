#include "common.h"

int main( int argc, const char* argv[] )
{
	srand(time(0));

	// Prints each argument on the command line.
	for( int i = 0; i < argc; i++ )
	{
		printf( "(scratch 8) arg %d: %s\n", i, argv[i] );
	}

	//psam *p = make_psam(6,4);
	//printf("(line 12) %d %d\n", p->nsites, p->nstates);
	//rand_init_psam(p);
	//print_psam(p);
	//shuffle_psam(p);
	//print_psam(p);

	settings *ss = make_settings();
	char urspath[] = "/Users/victor/Applications/simgenome-c-beta/examples/five.urs";
	ss->urspath = urspath;
	char psampath[] = "/Users/victor/Applications/simgenome-c-beta/examples/five.psam";
	ss->psampath = psampath;

	t_gene** mygenes;
	int ngenes;
	mygenes = read_genes_from_file(ss, ngenes);

	printf("(scratch 30) I found %d genes\n", ngenes);
	for(int ii=0; ii<ngenes; ii++){
		printf("\n Gene %d:\n", ii);
		print_urs( mygenes[ii]->urs, mygenes[ii]->urslen);
		print_psam( mygenes[ii]->dbd );
	}

	return 1;

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

	t_pop *pop;
	pop = make_population(4, -1, ss);
	char outpath[] = "./psams.txt";
	write_psams(pop->genomes[0], ss, outpath);

	t_ptable *pt;
	pt = make_ptable(4, 3, 1000);
	//print_ptable(pt);

}
