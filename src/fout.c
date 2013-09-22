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

	tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat( strcat(tmp, ss->outdir), "/URS/");
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
		fprintf(fp, "\nGenome %d COOPs:\n", g->id);
	}
	fprintf(fp, "a->\tb\t@d\tmultiplier\n");
	for (int ii = 0; ii < g->ntfs; ii++){
		for (int jj = 0; jj < g->ntfs; jj++){
			for (int dd = 0; dd < ss->maxgd; dd++){
				//printf("\n. fout 260 ntfs %d ii %d jj %d dd %d", g->ntfs, ii, jj, dd);
				fprintf(fp, "%d\t%d\t%d\t%f\n", ii, jj, dd, g->genes[ii]->gamma[jj*ss->maxgd + dd] );
			}
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
		fprintf(fp, "\nGenome %d URSs:\n", g->id);
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
		fprintf(fp, "\nGenome %d PSAMs\n", g->id);
	}
	for (int ii = 0; ii < g->ngenes; ii++){
		if (g->genes[ii]->has_dbd == true){
			psam *x = g->genes[ii]->dbd;
			fprintf(fp, "Gene %d %d\n", ii, g->genes[ii]->reg_mode);
			for (int jj = 0; jj < x->nsites; jj++){
				fprintf(fp, "%f %f %f %f\n",
						x->data[ii*x->nstates],
						x->data[ii*x->nstates+1],
						x->data[ii*x->nstates+2],
						x->data[ii*x->nstates+3]
						);
			}
		}
	}

	if (fo == NULL){
		fclose(fp);
	}
	free(gs);
	free(gc);
}

/* Writes the population and all its genomes and genes to a file. */
void log_population(t_pop* pop, settings* ss){
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
		fprintf(fo, "Genome: %d\n", ii);
		for (int jj = 0; jj < pop->genomes[ii]->ngenes; jj++){
			fprintf(fo, "ID %d gene %d %s %d %d\n", //%s %d\",
					ii, jj,
					(pop->genomes[ii]->genes[jj]->has_dbd)?"regulator":"reporter",
					(pop->genomes[ii]->genes[jj]->reg_mode),
					pop->genomes[ii]->genes[jj]->urslen);
		}
	}

	for (int ii = 0; ii < pop->ngenomes; ii++){
		log_urs(pop->genomes[ii], ss, fo);
		log_dbds(pop->genomes[ii], ss, fo);
		log_cofactor(pop->genomes[ii], ss, fo);
	}

	fclose(fo);
	free(gc);
	free(p);
}

