#include "common.h"

int main( int argc, const char* argv[] )
{
	srand(time(0));

	// Prints each argument on the command line.
	for( int i = 0; i < argc; i++ )
	{
		printf( "(scratch 8) arg %d: %s\n", i, argv[i] );
	}

	//printf("(scratch 11)\n");

	psam *p = make_psam(10,4);
	//printf("(line 12) %d %d\n", p->nsites, p->nstates);
	rand_init(p);
	print_psam(p);

}
