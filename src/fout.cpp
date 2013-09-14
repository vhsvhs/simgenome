#include "common.h"

void build_output_folders(settings* ss){
	char *tmp;
	tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat( strcat(tmp, ss->outdir), "/FITNESS/");
	//printf("\n fout 7%s\n", tmp);
	if (!Filexists(tmp)){
		mkdir(tmp, 0700);
	}
	free(tmp);

	tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat( strcat(tmp, ss->outdir), "/LOGS/");
	if (!Filexists(tmp)){
		mkdir(tmp, 0700);
	}
	free(tmp);

	tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat( strcat(tmp, ss->outdir), "/EXPR/");
	if (!Filexists(tmp)){
		mkdir(tmp, 0700);
	}
	free(tmp);
}

/* Write all the PSAMs in genome gn to the file located at outpath. */
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


void log_fitness(double* f, int len, settings* ss){
	//printf("fout 33\n");
	char* g;
	g = (char*)malloc(10*sizeof(char));
	sprintf(g, "%d", ss->gen_counter);

	//printf("fout 37\n");

	char* p = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat(
		strcat(
			strcat(
					strcat(p, ss->outdir),
			"/FITNESS/fitness.gen"),
		g),
	".txt");

	//printf("fout 46\n");

	FILE *fp;
	fp = fopen(p, "w");
	if (fp == NULL) {
	  fprintf(stderr, "Error: can't open output file %s!\n",
			  p);
	  exit(1);
	}


	for (int ii=0; ii<len; ii++){
		//printf("%d %f\n", ii, f[ii]);
		fprintf(fp, "ID %d, %f\n", ii, f[ii]);
	}
	fclose(fp);
	free(g);
	free(p);
}

void log_expr(t_genome *g, int t, t_ptable* ptable, settings* ss){
	char* gc;
	gc = (char*)malloc(10*sizeof(char));
	sprintf(gc, "%d", ss->gen_counter);

	char* gs;
	gs = (char*)malloc(4*sizeof(char));
	sprintf(gs, "%d", g->id);

	char* ts;
	ts = (char*)malloc(100*sizeof(char));
	sprintf(ts, "%d", t);

	char* p = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat(
		strcat(
			strcat(
			strcat(
			strcat(
			strcat(p, ss->outdir),
			"/EXPR/expr.gen"),
			gc),
		".id"),
		gs),
	".txt");

	FILE *fp;
	if (t == 0){
		fp = fopen(p, "w");
	}
	else{
		fp = fopen(p, "a");
	}
	if (fp == NULL) {
	  fprintf(stderr, "Error: can't open output file %s!\n",
			  p);
	  exit(1);
	}



	for (int jj = 0; jj < g->ngenes; jj++){
		char* gs;
		gs = (char*)malloc(100*sizeof(char));
		sprintf(gs, "%d", jj);

		char *header = (char *)malloc(MAXLEN*sizeof(char));
		strcat(
			strcat(
				strcat(
				strcat(
					strcat( header, ". t "),
					ts),
				" gene "),
				gs),
			"\n");
		//printf("fout.cpp genome %d 150: %s\n", g->id, header);
		fprintf(fp, "%s", header);

		// Lines look like this:
		// site 1 :        0 a 0.632     1 r 0.368
		// where a is for activator
		// and r is for repressor
		for (int site = 0; site <  g->genes[jj]->urslen; site++){
			char* sss = (char*)malloc(10*sizeof(char));
			sprintf(sss, "%d", site);
			fprintf(fp, "site %s :", sss);

			double sump = ptable->cpr[site];
			for (int tf = 0; tf < g->ntfs; tf++){
				char* tfs = (char*)malloc(10*sizeof(char));
				sprintf(tfs, "%d", tf);

				char* rms = (char*)malloc(10*sizeof(char));
				if (g->genes[jj]->has_dbd && g->genes[jj]->reg_mode == 0){
					sprintf(rms, "r");
				}
				else if (g->genes[jj]->has_dbd && g->genes[jj]->reg_mode == 1){
					sprintf(rms, "a");
				}

				double sumt = ptable->cpt[site*ptable->M + tf];

				//if (sumt < 100 || sumt > 100){
				//	printf("fout:177 %d %d %d %f\n",site, tf, ptable->M, sumt);
				//}

				double this_val = 0.0;
				if (sump > 0.0 && sumt > 0.0){
					this_val = sumt/sump;
				}

				fprintf(fp, "\t%s", tfs);
				fprintf(fp, " %s", rms);
				fprintf(fp, " %f", this_val);

				free(tfs);
				free(rms);
			}
			fprintf(fp, "\n");

		}

		/* Reset the header */
		for (int qq = 0; qq < MAXLEN; qq++){
			header[qq] = NULL;
		}

		free(header);
		free(gs);
	}
	free(ts);
	free(gc);
	free(gs);
	free(p);

	fclose(fp);


}

void log_cofactor(t_genome *g, settings* ss){

}


void log_dbds(t_genome *g, settings* ss){

}

