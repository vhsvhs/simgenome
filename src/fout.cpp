#include "common.h"

void write_psams(t_genome *gn, settings *ss, char* outpath){
	FILE *fp;
	fp = fopen(outpath,"w");
	if (fp == NULL) {
	  fprintf(stderr, "Error: can't open output file %s!\n",
			  outpath);
	  exit(1);
	}

	for (int gg=0; gg < gn->ngenes; gg++){
		if (gn->genes[gg]->has_dbd) {
			psam *this_dbd = gn->genes[gg]->dbd;

			fprintf(fp, "Gene %d\n", gg);
			for (int ii=0; ii < this_dbd->nsites; ii++){ // sites
				for (int jj=0; jj < this_dbd->nstates; jj++){ //states
					char state = int2nt( jj );
					double value = this_dbd->data[ii*this_dbd->nstates + jj];
					fprintf(fp, "%c (%f)\t", state, value );
				}
				fprintf(fp, "\n");
			}
		}
	}
	fclose(fp);
}
