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
	strcat( strcat(tmp, ss->outdir), "/OCCUPANCY/");
	if (!Filexists(tmp)){
		mkdir(tmp, 0700);
	}
	free(tmp);

	tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat( strcat(tmp, ss->outdir), "/COOP/");
	if (!Filexists(tmp)){
		mkdir(tmp, 0700);
	}
	free(tmp);

	tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat( strcat(tmp, ss->outdir), "/DBD/");
	if (!Filexists(tmp)){
		mkdir(tmp, 0700);
	}

	tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat( strcat(tmp, ss->outdir), "/POPS/");
	if (!Filexists(tmp)){
		mkdir(tmp, 0700);
	}

	tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat( strcat(tmp, ss->outdir), "/MATING/");
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

void log_occupancy(t_genome *g, int gid, int t, t_ptable* ptable, settings* ss){
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
				strcat( header, ". t "),
				ts),
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
	free(p);

	fclose(fp);
}

void log_cofactor(t_genome *g, settings* ss){
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
	fp = fopen(p, "w");
	if (fp == NULL) {
	  fprintf(stderr, "Error: can't open output file %s!\n",
			  p);
	  exit(1);
	}

	//printf("\n fout 256 g->ntfs = %d %s", g->ntfs, p);
	fprintf(fp, "a->\tb\t@d\tmultiplier\n");
	for (int ii = 0; ii < g->ntfs; ii++){
		for (int jj = 0; jj < g->ntfs; jj++){
			for (int dd = 0; dd < ss->maxgd; dd++){
				//printf("\n. fout 260 ntfs %d ii %d jj %d dd %d", g->ntfs, ii, jj, dd);
				fprintf(fp, "%d\t%d\t%d\t%f\n", ii, jj, dd, g->genes[ii]->gamma[jj*ss->maxgd + dd] );
			}
		}
	}

	free(gc);
	free(gs);
	free(p);

	fclose(fp);
}


void log_dbds(t_genome *g, settings* ss){
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
	fp = fopen(p, "w");
	if (fp == NULL) {
	  fprintf(stderr, "Error: can't open output file %s!\n",
			  p);
	  exit(1);
	}


	for (int ii = 0; ii < g->ngenes; ii++){
		if (g->genes[ii]->has_dbd == true){
			psam *x = g->genes[ii]->dbd;
			fprintf(fp, "Gene %d\n", ii);
			for (int jj = 0; jj < x->nsites; jj++){
				fprintf(fp, "site %d A: %f C: %f G: %f T: %f\n",
						jj,
						x->data[ii*x->nstates],
						x->data[ii*x->nstates+1],
						x->data[ii*x->nstates+2],
						x->data[ii*x->nstates+3]
						);
			}
		}
	}

	fclose(fp);
	free(gs);
	free(gc);
	free(p);
}
