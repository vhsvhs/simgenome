#include "common.h"

/* Makes memory for one gene
 *
 * NOTE: The returned gene will not yet have memory allocated
 * for the gamma array (i.e. the co-factor affinity matrix).
 * The gamma array can be build manualy, or by using the
 * function build_coop followed by init_coop.
 * */
t_gene* make_gene(int psamlen, int urslen) {
	t_gene *g;
	g = (t_gene *)malloc(1*sizeof(t_gene));
	g->name = (char *)malloc(GENE_NAME_MAX*sizeof(char));
	g->urslen = urslen;
	g->urs = (int *)malloc(urslen*sizeof(int));
	if (psamlen > 0){
		g->has_dbd = true;
		g->dbd = make_psam(psamlen, N_STATES);
	}
	else{
		g->has_dbd = false;
	}
	g->gammalen = 0;
	g->tfcooplen = 0;
	g->reg_mode = 0;
	return g;
}

/* Copies a gene 'from' into an entirely new gene, which is returned */
t_gene* copy_gene(t_gene* from) {
	t_gene* to = (t_gene*)malloc(sizeof(t_gene));

	to->id = from->id;
	// This strcpy is safe because gene names are always GENE_NAME_MAX characters long.

	/* name */
	to->name = (char *)malloc(GENE_NAME_MAX*sizeof(char));
	strcpy( to->name, from->name );

	/* urs */
	to->urs = (int *)malloc(from->urslen*sizeof(int));
	for(int ii=0; ii<from->urslen; ii++){
		to->urs[ii] = from->urs[ii];
	}
	to->urslen = from->urslen;

	/* psam */
	to->has_dbd = from->has_dbd;
	if (to->has_dbd == true){
		to->dbd = make_psam(from->dbd->nsites, from->dbd->nstates);
		copy_psam( to->dbd, from->dbd );
	}
	to->reg_mode = from->reg_mode;


	if (to->has_dbd == true){
		if (from->tfcoop > 0) {
			/* We could, alternatively, use the method "build_coop" here,
			 * but that would require that we know how many TFs and the max
			 * cofactor distance... both are variables not stores in the t_gene
			 * struct.  Rather, we'll just allocate arrays that are the same
			 * size as the 'from' gene.
			 */
			to->gamma = (double*)malloc(from->gammalen*sizeof(double));
			to->tfcoop = (double*)malloc(from->tfcooplen*sizeof(double));
		}
		to->gammalen = from->gammalen;
		to->tfcooplen = from->tfcooplen;

		for(int ii=0; ii < from->gammalen; ii++){
			to->gamma[ii] = from->gamma[ii];
		}

		for(int ii=0; ii < from->tfcooplen; ii++){
			to->tfcoop[ii] = from->tfcoop[ii];
		}
	}
	else{
		to->gammalen = 0;
		to->tfcooplen = 0;
	}


	return to;

}

void build_coop(t_gene* g, int ntfs, int maxd){
	g->tfcooplen = ntfs;
	g->tfcoop = (double*)malloc(g->tfcooplen*sizeof(double));
	g->gammalen = ntfs*maxd;
	g->gamma = (double*)malloc(g->gammalen*sizeof(double));
}

void init_coop(t_gene* g){
	for (int ii = 0; ii < g->tfcooplen; ii++){
		g->tfcoop[ii] = 0.0;
	}
	for (int ii = 0; ii < g->gammalen; ii++){
		g->gamma[ii] = 0.0;
	}
}

/* f is the relative affinity for cofactor f.
 * d is the distance to the cofactor.
 */
double coopfunc(double f, int d){
	double val = 1 + f * exp( (-1)*(d^2)/V_RATE_OF_COOP_DECAY );
	//printf("gene 105 coopfunc, d= %d ret= %f\n", d, val);
	return val;
}

/* Recalculates the gamma array, given values in the tfcoop array.
 * ntfs = n possible cofactors
 * maxd = the maximum distance over which cofactor interacions can occur.
 *
 * This method assumes that you've already allocated memory
 * for g->tfcoop and g->gamma.
 */
void calc_gamma(t_gene* g, int ntfs, int maxd){
	for (int ii = 0; ii < ntfs; ii++){
		for (int jj = 0; jj < maxd; jj++){
			g->gamma[ii*maxd + jj] = coopfunc( g->tfcoop[ii], jj );
			//printf("\n calc_gamma tf %d %d %f\n", ii, jj, g->gamma[ii*maxd + jj]);
		}
	}
}

/* Garbage collection for a gene */
void free_gene(t_gene* g){
	free(g->name);
	free(g->urs);
	if (g->has_dbd){
		free_psam(g->dbd);
		free(g->dbd);
		free(g->gamma);
		free(g->tfcoop);
	}
}



void print_urs(int* urs, int urslen) {
	char state;
	state = '-';
	printf("URS: ");
	for (int ii=0; ii<urslen; ii++){
		if (urs[ii] == 0)   {state='A';}
		else if (urs[ii]==1){state='C';}
		else if (urs[ii]==2){state='G';}
		else if (urs[ii]==3){state='T';}
		printf("%c", state);
	}
	printf("\n");
}

int* get_random_seq(int len){
	int *seq;
	seq = (int *)malloc(len*sizeof(int));
	for (int ii=0; ii<len; ii++) {
		double rand_f = (float)rand()/(float)RAND_MAX;
		seq[ii] = (int)( rand_f * N_STATES );
	}
	return seq;
}

/* Returns an array of pointers to genes.
 * Upon completion, ngenes is also set to the correct value. */
t_gene** read_genes_from_file(settings *ss, int &ngenes) {

	FILE *fp; /* File for psam specs */
	fp = fopen(ss->psampath,"r");
	if (fp == NULL) {
	  fprintf(stderr, "Error: can't open PSAM file %s!\n",
			  ss->psampath);
	  exit(1);
	}

	FILE *fu; /* File for urs specs */
	fu = fopen(ss->urspath,"r");
	//printf("\n gene 90: urspath=%s\n", ss->urspath);
	if (fu == NULL) {
	  fprintf(stderr, "Error: can't open URS file %s!\n",
			  ss->urspath);
	  exit(1);
	}

	/* How many genes?
	 * ...count the FASTA-formatted URS definitions */
	ngenes = 0;
	char line[MAXLEN];
	while (  fgets(line, MAXLEN, fu)  ){
		if (line[0] == '>'){
			ngenes += 1;
		}
		const char* tokens[MAX_TOKENS] = {};
		tokens[0] = strtok(line, " ");
		if (tokens[0]){ // zero if line is blank
			for (int ii=0; ii<MAX_TOKENS; ii++){
				tokens[ii] = strtok(0, " "); //subsequent tokens
				if (!tokens[ii]){ // No more tokens
					break;
				}
			}
		}
	}

	/* Get URSes */
	rewind(fu); /* Set the fu file pointer back to the start */
	int *urslengths = (int *)malloc(ngenes*sizeof(int));
	char **urses = (char **)malloc(ngenes*sizeof(char*));
	int this_gene = -1;
	while (  fgets(line, MAXLEN, fu) ){
		if (line[0] == '>'){
			this_gene += 1;
		}
		else {
			int ii = 0;
			int count = 0;
			while (line[ii] == 'A'
					|| line[ii] == 'C'
					|| line[ii] == 'T'
					|| line[ii] == 'G')
			{
				count++;
				ii++;
			}
			if (this_gene < ngenes){
				//printf("gene 133, this_gene %d count= %d\n", this_gene, count);
				urslengths[this_gene] = count;
				urses[this_gene] = (char *)malloc(urslengths[this_gene]*sizeof(char));
				for(int ii=0; ii< urslengths[this_gene]; ii++){
					urses[this_gene][ii] = line[ii];
				}
			}
		}
	}
	fclose(fu);

	if (ss->verbosity > 20){
		printf("\n. I found %d gene definitions.\n", ngenes);
	}

	/* Get PSAMS */
	double **psams = (double **)malloc(ngenes*sizeof(double*)); // key = gene ID, value = 1-d array of doubles
	bool *has_dbd = (bool *)malloc(ngenes*sizeof(bool));
	int *psamlengths = (int *)malloc(ngenes*sizeof(int)); // key = gene ID, value = length of psam
	bool *reg_modes = (bool *)malloc(ngenes*sizeof(bool));
	for(int ii=0; ii<ngenes; ii++){
		has_dbd[ii] = false;
		psamlengths[ii] = 0;
		reg_modes[ii] = 0;
	}

	this_gene = -1;
	int count = 0;
	while ( fgets(line, MAXLEN, fp)  ){
		if (line[0] == '#'){
			continue;
		}

		//printf("(settings 90) line=%s\n", line);
		char * token = strtok(line, " ");
		//printf("(settings 92) line=%s\n", line);
		if (token){ // zero if line is blank
			if (token[0] == '\n'){
				continue;
			}
			else if (token[0] == '#') {
				continue;
			}
			else if (token[0] == 'G'
					&& token[1] == 'e'
					&& token[2] == 'n'
					&& token[3] == 'e') {
				token = strtok(NULL, " ");
				this_gene = atoi( token );
				if (this_gene < ngenes){
					//printf("\n. gene 272 - adding dbd to gene %d\n", this_gene);
					has_dbd[this_gene] = true;
					token = strtok(NULL, " ");
					reg_modes[ this_gene ]= (bool)atoi( token );
					count = 0;
					psams[this_gene] = (double *)malloc(MAX_PSAM_LEN*sizeof(double)); //reset this_psam
					psamlengths[this_gene] = 0;
				}
			}
			else if (this_gene < ngenes){ // a line with 4 floating-point numbers
				//printf("\n. gene 282 reading PSAM for gene %d\n", this_gene, token);

				// State 0
				float value = atof( token );
				//printf("(settings 112a) value=%f this_gene=%d count=%d\n", value, this_gene, count);
				psams[this_gene][count] = value;
				count += 1;

				// States 1,2,3
				for (int ii=1; ii<4; ii++){ // ii = state
					value = atof( strtok(NULL, " ") ); //subsequent tokens
					//printf("(settings 112b) value=%f this_gene=%d count=%d\n", value, this_gene, count);
					psams[this_gene][count] = value;
					count++;
				}
				psamlengths[this_gene]++;
			}
		}
	}
	fclose(fp);

	/* genes is an array of t_gene*
	 * We'll fill this array with t_gene objects. */
	t_gene** genes = (t_gene**)malloc(ngenes*sizeof(t_gene));

	/* Build the genes, using the URSes and PSAMs we previously found. */
	int ntfs = 0;
	for (int ii=0; ii<ngenes; ii++){ // ii = gene
		genes[ii] = make_gene(psamlengths[ii], urslengths[ii]);
		genes[ii]->id = ii;
		for (int jj=0; jj<urslengths[ii]; jj++){ // jj = site
			genes[ii]->urs[jj] = nt2int( urses[ii][jj] );
		}
		if (has_dbd[ii] == true){
			ntfs++;
			for (int jj=0; jj<N_STATES*psamlengths[ii]; jj++){
				genes[ii]->dbd->data[jj] = psams[ii][jj];
				genes[ii]->has_dbd = true;
				genes[ii]->reg_mode = reg_modes[ii];
			}
		}
		if (ss->verbosity > 10){
			printf("\n+ Gene %d, URS (%d sites)", ii, genes[ii]->urslen);
			if (genes[ii]->has_dbd == true){
				printf(", PSAM (%d sites)", genes[ii]->dbd->nsites);
			}
			printf("\n");
		}
	}

	/* Setup cofactor variables */
	for (int ii=0; ii<ntfs; ii++){
		build_coop(genes[ii], ntfs, ss->maxgd);
		init_coop(genes[ii]);
		calc_gamma(genes[ii], ntfs, ss->maxgd);
	}

	/* Cleanup Memory */
	for (int ii=0; ii < ngenes; ii++){
		free(urses[ii]);
	}
	//for (int ii = 0; ii < ntfs; ii++){
	//	free(psams[ii]);
	//}
	free(urses);
	//free(psams);
	free(urslengths);
	free(psamlengths);
	free(reg_modes);
	free(has_dbd);

	return genes;
}
