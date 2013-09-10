#include "common.h"

int main( int argc, const char* argv[] )
{
	srand(time(0));

	// Prints each argument on the command line.
	for( int i = 0; i < argc; i++ )
	{
		printf( "(scratch 8) arg %d: %s\n", i, argv[i] );
	}

	psam *p = make_psam(6,4);
	//printf("(line 12) %d %d\n", p->nsites, p->nstates);
	rand_init(p);
	print_psam(p);
	shuffle_psam(p);
	print_psam(p);

	settings *ss = make_settings();
	printf("(scratch 23) %f %d %f\n", ss->pwmlenmu, ss->pwmlenmumax, ss->ddgmu);

	ss->verbosity = 101;
	mutate_psam(p, ss);
	printf("After mutate_psam:\n");
	print_psam(p);

	int seqlen = 6;
	int *seq;
	seq = (int *)malloc(seqlen*sizeof(int));

	for (int ii=0; ii<seqlen; ii++) {
		seq[ii] = (int)((float)rand()/(float)RAND_MAX);
	}
	double aff = get_affinity(p, seq, 6);
	printf("(scratch 30) aff = %f\n", aff);
}
