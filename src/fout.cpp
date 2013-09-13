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

void log_expr(t_genome *g, int gid, int t, settings* ss){
	char* gc;
	gc = (char*)malloc(10*sizeof(char));
	sprintf(gc, "%d", ss->gen_counter);

	char* ts;
	ts = (char*)malloc(4*sizeof(char));
	sprintf(ts, "%d", t);

	char* gs;
	gs = (char*)malloc(4*sizeof(char));
	sprintf(gs, "%d", gid);

	char* us;
	us = (char*)malloc(g->genes[gid]->urslen*sizeof(char));
	for (int ii = 0; ii < g->genes[gid]->urslen; ii++){
		us[ii] = int2nt( g->genes[gid]->urs[ii] );
	}

	char* p = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat(
		strcat(
			strcat(
			strcat(
			strcat(
			strcat(
			strcat(
					strcat(p, ss->outdir),
					"/EXPR/expr.gen"),
			gc),
		".t"),
		ts),
		".id"),
		gs),
	".txt");

	FILE *fp;
	fp = fopen(p, "w");
	if (fp == NULL) {
	  fprintf(stderr, "Error: can't open output file %s!\n",
			  p);
	  exit(1);
	}

	char *header = (char *)malloc(MAXLEN*sizeof(char));
	strcat(
	strcat(
	strcat(
	strcat(
	strcat(
	strcat(
	strcat( header, ". t "),
	ts),
	" gene "),
	gs),
	" URS: "),
	us),
	"\n");
	fprintf(fp, "%s", header);

	fclose(fp);
	free(gc);
	free(ts);
	free(gs);
	free(us);
	free(p);
}

void log_cofactor(t_genome *g, settings* ss){

}


void log_dbds(t_genome *g, settings* ss){

}

