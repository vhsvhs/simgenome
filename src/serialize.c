#include "common.h"

/* Writes the population and all its genomes and genes to a file. */
void serialize_population(t_pop* pop, settings* ss){
	char* gc;
	gc = (char*)malloc(10*sizeof(char));
	sprintf(gc, "%d", ss->gen_counter);

	char* p = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat(
			strcat(
			strcat(
			strcat(p, ss->outdir),
			"/POPS/pop.gen"),
			gc),
			".save.txt"
			);

	FILE *fo; /* File for psam specs */
	fo = fopen(p,"w");
	if (fo == NULL) {
	  fprintf(stderr, "Error: can't open output file %s!\n",
			  p);
	  exit(1);
	}

	if (ss->verbosity > 2){
		printf("\n. The population was saved to %s\n", p);
	}

	fprintf(fo, "Hello World.\n");

	fclose(fo);
	free(gc);
	free(p);
}

t_pop* deserialize_population(char* fpath, settings* ss){
	t_pop* retp;
	return retp;
}
