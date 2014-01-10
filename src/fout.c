#include "common.h"

void build_output_folders(settings* ss){
	if (ss->verbosity > 3){
		printf("\n. I'm building output folders in %s\n", ss->outdir);
	}

	if (!Filexists(ss->outdir)){
		mkdir(ss->outdir, 0700);
	}

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
	strcat( strcat(tmp, ss->outdir), "/MUTATIONS/");
	if (!Filexists(tmp)){
		mkdir(tmp, 0700);
	}
	free(tmp);

	tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat(strcat(tmp, ss->outdir), "/LOGS/expression.txt");
	ss->file_expr_log = fopen(tmp, "w");
	if (ss->file_expr_log == NULL) {
	  fprintf(stderr, "Error: can't open output file %s!\n", tmp);
	  fprintf(stderr, "fout 121: %s\n", ss->outdir);
	  exit(1);
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

	if (ss->enable_timelog){
		tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
		strcat(strcat(tmp, ss->outdir), "/LOGS/time.txt");
		ss->file_time_log = fopen(tmp, "w");
		if (ss->file_time_log == NULL) {
		  fprintf(stderr, "Error: can't open output file %s!\n", tmp);
		  fprintf(stderr, "fout 72: %s\n", ss->outdir);
		  exit(1);
		}
		free(tmp);
	}

}


/* Writes the file FITNESS/fitness.genX.txt for generation X,
 * and updates the file LOGS/generations.txt for all generations so far.
 */
void log_fitness(double* f, double* er, int len, settings* ss){

	double minf = min(f, len);
	double maxf = max(f, len);
	double meanf = mean(f, len);
	double errf = sderr(f, len);

	/*
	 * Part 1: write FITNESS/fitness.genX.txt */
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
	char mark = ' ';
	for (int ii=0; ii<len; ii++){
		double this_f = f[ii];
		if (this_f == maxf){ mark = '*'; }
		else{ mark = ' '; }
		fprintf(fp, "%d %f %f %c\n", ii, f[ii], er[ii], mark);
	}
	fclose(fp);
	if (ss->verbosity > 2){
		printf("\n. Fitness values were written to %s\n", p);
	}
	free(g);
	free(p);

	/*
	 * Part 2: update LOGS/generations.txt */
	p = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat( strcat(p, ss->outdir), "/LOGS/generations.txt");
	FILE *fo;

	/* Open a new generation log file for the first generation. */
	if (ss->gen_counter == ss->start_gen){
		fo = fopen(p, "w");
	}
	/* Otherwise, append to the existing generation log. */
	else{
		fo = fopen(p, "a");
	}
	if (fo == NULL) {
	  fprintf(stderr, "Error: can't open output file %s!\n", p);
	  free(p);
	  exit(1);
	}

	// gen, max, min, mean
	fprintf(fo, "gen: %d\tmax= %f\tmin= %f\tmean= %f\terr= %f\n",
			ss->gen_counter, maxf, minf, meanf, errf);

	fclose(fo);
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
		fprintf(fp, "site %d (%c):", site, int2nt(g->genes[gid]->urs[site]) );

		double sump = ptable->cpr[site];
		for (int tf = 0; tf < g->ntfs; tf++){
			char* tfs = (char*)malloc(10*sizeof(char));
			sprintf(tfs, "%d", tf);

			char* rms = (char*)malloc(10*sizeof(char));
			if (g->genes[tf]->has_dbd && g->genes[tf]->reg_mode == 0){
				sprintf(rms, "r");
			}
			else if (g->genes[tf]->has_dbd && g->genes[tf]->reg_mode == 1){
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
		header[qq] = 0;
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

void log_timegen(settings* ss, int gen){
	if (ss->file_time_log != NULL) {
		if (gen == 0){
			fprintf(ss->file_time_log,
					"gen.\tt_tot\tt_f\tt_makept\tt_fillpt\tt_samplept\n");
		}
		int delta = ss->t_stopgen - ss->t_startgen;
		float secs = ((float)delta)/CLOCKS_PER_SEC;
		float sumf = ((float)ss->t_sumf)/CLOCKS_PER_SEC;
		float makept = ((float)ss->t_summakept)/CLOCKS_PER_SEC;
		float fillpt = ((float)ss->t_sumfillpt)/CLOCKS_PER_SEC;
		float samplept = ((float)ss->t_sumsamplept)/CLOCKS_PER_SEC;
		fprintf(ss->file_time_log,
				"%d\t%f\t%f\t%f\t%f\t%f\n", gen, secs, sumf, makept, fillpt, samplept);
	}
}

