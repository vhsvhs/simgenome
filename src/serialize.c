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

	fprintf(fo, "N genomes: %d\n", pop->ngenomes);

	for (int ii = 0; ii < pop->ngenomes; ii++){
		for (int jj = 0; jj < pop->genomes[ii]->ngenes; jj++){
			fprintf(fo, "ID %d gene %d %b %d %s %d\n",
					ii, jj, pop->genomes[ii]->genes[jj],
					(pop->genomes[ii]->genes[jj]->has_dbd)?"true":"false",
					pop->genomes[ii]->genes[jj]->reg_mode,
					pop->genomes[ii]->genes[jj]->urs,
					pop->genomes[ii]->genes[jj]->urslen);
		}
	}


	fclose(fo);
	free(gc);
	free(p);
}

t_pop* deserialize_population(char* fpath, settings* ss){
	t_pop* retp;
	return retp;
}
