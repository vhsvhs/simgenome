#include "common.h"

void build_output_folders(settings* ss){
	char *tmp;
	tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat( strcat(tmp, ss->outdir), "/FITNESS/");
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
	strcat( strcat(tmp, ss->outdir), "/OCCUPANCY/");
	if (!Filexists(tmp)){
		mkdir(tmp, 0700);
	}
	free(tmp);

	tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat( strcat(tmp, ss->outdir), "/POPS/");
	if (!Filexists(tmp)){
		mkdir(tmp, 0700);
	}
	free(tmp);

	tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat( strcat(tmp, ss->outdir), "/MATING/");
	if (!Filexists(tmp)){
		mkdir(tmp, 0700);
	}
	free(tmp);

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

	FILE *fp;
	fp = fopen(p, "w");
	if (fp == NULL) {
	  fprintf(stderr, "Error: can't open output file %s!\n",
			  p);
	  exit(1);
	}

	for (int ii=0; ii<len; ii++){
		//printf("%d %f\n", ii, f[ii]);
		fprintf(fp, "%d %f\n", ii, f[ii]);
	}
	fclose(fp);
	if (ss->verbosity > 10){
		printf("\n. Fitness values were written to %s\n", p);
	}
	free(g);
	free(p);
}

void log_occupancy(t_genome *g, int gid, int t, int rid, t_ptable* ptable, settings* ss){
	char* gc;
	gc = (char*)malloc(10*sizeof(char));
	sprintf(gc, "%d", ss->gen_counter);

	char* gs;
	gs = (char*)malloc(4*sizeof(char));
	sprintf(gs, "%d", g->id);

	char *ge;
	ge = (char*)malloc(10*sizeof(char));
	sprintf(ge, "%d", gid);

	char* ts;
	ts = (char*)malloc(100*sizeof(char));
	sprintf(ts, "%d", t);

	char* ris;
	ris = (char*)malloc(100*sizeof(char));
	sprintf(ris, "%d", rid);

	char* p = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat(
		strcat(
		strcat(
		strcat(
			strcat(
			strcat(
			strcat(
			strcat(p, ss->outdir),
			"/OCCUPANCY/occ.gen"),
			gc),
		".id"),
		gs),
		".gene"),
		ge),
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


	char *header = (char *)malloc(MAXLEN*sizeof(char));
	strcat(
		strcat(
			strcat(
			strcat(
				strcat( header, ". time "),
				ts),
				" rid "),
				ris),
		"\n");
	//printf("fout.cpp genome %d 150: %s\n", g->id, header);
	fprintf(fp, "%s", header);

	// Lines look like this:
	// site 1 :        0 a 0.632     1 r 0.368
	// where a is for activator
	// and r is for repressor
	for (int site = 0; site <  g->genes[gid]->urslen; site++){
		char* sss = (char*)malloc(10*sizeof(char));
		sprintf(sss, "%d", site);
		fprintf(fp, "site %s :", sss);

		double sump = ptable->cpr[site];
		for (int tf = 0; tf < g->ntfs; tf++){
			char* tfs = (char*)malloc(10*sizeof(char));
			sprintf(tfs, "%d", tf);

			char* rms = (char*)malloc(10*sizeof(char));
			if (g->genes[gid]->has_dbd && g->genes[gid]->reg_mode == 0){
				sprintf(rms, "r");
			}
			else if (g->genes[gid]->has_dbd && g->genes[gid]->reg_mode == 1){
				sprintf(rms, "a");
			}

			double sumt = ptable->cpt[site*ptable->M + tf];


			double this_val = 0.0;
			if (sump > 0.0 && sumt > 0.0){
				this_val = sumt/sump;
			}
			else if (sump < 0.0 || sumt < 0.0){
				this_val = 0.0;
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
	free(ts);
	free(gc);
	free(ge);
	free(ris);
	free(p);

	fclose(fp);
}


void log_cofactor(t_genome *g, settings* ss, FILE* fo){
	char* gc;
	gc = (char*)malloc(10*sizeof(char));
	sprintf(gc, "%d", ss->gen_counter);

	char* gs;
	gs = (char*)malloc(4*sizeof(char));
	sprintf(gs, "%d", g->id);

	char* p = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat(
		strcat(
			strcat(
			strcat(
			strcat(
			strcat(p, ss->outdir),
			"/COOP/coop.gen"),
			gc),
		".id"),
		gs),
	".txt");

	FILE *fp;
	if (fo == NULL) {
		fp = fopen(p, "w");
		if (fp == NULL) {
		  fprintf(stderr, "Error: can't open output file %s!\n",
				  p);
		  exit(1);
		}
	}
	else{
		fp = fo;
	}

	if (fo != NULL) {
		fprintf(fp, "\nCOOP Genome %d :\n", g->id);
	}

	for (int ii = 0; ii < g->ntfs; ii++){
		for (int jj = 0; jj < g->ntfs; jj++){
			fprintf(fp, "%d %d %f\n", ii, jj, g->genes[ii]->tfcoop[jj]);
		}
	}

	if (fo == NULL) {
		fclose(fp);
	}
	free(gc);
	free(gs);
}

/* Writes the URSs for genome g.
 * If fo is not NULL, then the data will be written to fo.
 */
void log_urs(t_genome *g, settings* ss, FILE* fo){
	char* gc;
	gc = (char*)malloc(10*sizeof(char));
	sprintf(gc, "%d", ss->gen_counter);
	char* gs;
	gs = (char*)malloc(4*sizeof(char));
	sprintf(gs, "%d", g->id);

	char* p = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat(
		strcat(
			strcat(
			strcat(
			strcat(
			strcat(p, ss->outdir),
			"/URS/urs.gen"),
			gc),
		".id"),
		gs),
	".txt");


	FILE* fp;
	if (fo == NULL) {
		fp = fopen(p, "w");
		if (fp == NULL) {
		  fprintf(stderr, "Error: can't open output file %s!\n",
				  p);
		  exit(1);
		}
	}
	else{
		fp = fo;
	}

	if (fo != NULL) {
		fprintf(fp, "\nURS Genome %d :\n", g->id);
	}
	for (int ii = 0; ii < g->ngenes; ii++){
		fprintf(fp, ">%d %s\n", ii, g->genes[ii]->name);

		for (int jj = 0; jj < g->genes[ii]->urslen; jj++){
			fprintf(fp, "%c", int2nt(g->genes[ii]->urs[jj]) );
		}
		fprintf(fp, "\n");
	}

	if (fo == NULL){
		fclose(fp);
	}
	free(p);
	free(gc);
	free(gs);
}

void log_dbds(t_genome *g, settings* ss, FILE* fo){
	char* gc;
	gc = (char*)malloc(10*sizeof(char));
	sprintf(gc, "%d", ss->gen_counter);

	char* gs;
	gs = (char*)malloc(4*sizeof(char));
	sprintf(gs, "%d", g->id);

	char* p = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat(
		strcat(
			strcat(
			strcat(
			strcat(
			strcat(p, ss->outdir),
			"/DBD/dbd.gen"),
			gc),
		".id"),
		gs),
	".txt");

	FILE *fp;
	if (fo == NULL) {
		fp = fopen(p, "w");
		if (fp == NULL) {
		  fprintf(stderr, "Error: can't open output file %s!\n",
				  p);
		  exit(1);
		}
	}
	else{
		fp = fo;
	}

	if (fo != NULL) {
		fprintf(fp, "\nPSAM Genome %d :\n", g->id);
	}
	for (int ii = 0; ii < g->ngenes; ii++){
		if (g->genes[ii]->has_dbd == true){
			fprintf(fp, "Gene %d %d\n", ii, g->genes[ii]->reg_mode);

			for (int jj = 0; jj < g->genes[ii]->dbd->nsites; jj++){
				for (int kk = 0; kk < g->genes[ii]->dbd->nstates; kk++){
					fprintf(fp, "%f ", g->genes[ii]->dbd->data[jj*g->genes[ii]->dbd->nstates + kk]);
				}
				fprintf(fp, "\n");
			}
		}
	}


	if (fo == NULL){
		fclose(fp);
	}
	free(gs);
	free(gc);
}



