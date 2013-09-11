#include "common.h"

/* Makes memory for one gene */
t_gene* make_gene(int psamlen, int urslen) {
	t_gene *g;
	g = (t_gene *)malloc(1*sizeof(t_gene));
	g->name = (char *)malloc(GENE_NAME_MAX*sizeof(char));
	g->urslen = urslen;
	g->urs = (int *)malloc(urslen*sizeof(int));
	g->dbd = make_psam(psamlen, N_STATES);
	return g;
}

/* Garbage collection for a gene */
void free_gene(t_gene* g){
	free(g->name);
	free(g->urs);
	free_psam(g->dbd);
	free(g->dbd);
}

/* Copies a gene, assuming that memory has already been allocated. */
void copy_gene(t_gene* to, t_gene* from) {
	to->id = from->id;
	// This strcpy is safe because gene names are always GENE_NAME_MAX characters long.
	strcpy( to->name, from->name );

	// Ensure that the URS lenghts match
	if (to->urslen != from->urslen){
		free(to->urs);
		to->urs = (int *)malloc(from->urslen*sizeof(int));
	}

	for(int ii=0; ii<from->urslen; ii++){
		to->urs[ii] = from->urs[ii];
	}

	to->urslen = from->urslen;
	to->has_dbd = from->has_dbd;
	copy_psam( to->dbd, from->dbd );
	to->reg_mode = from->reg_mode;
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
	if (fu == NULL) {
	  fprintf(stderr, "Error: can't open URS file %s!\n",
			  ss->urspath);
	  exit(1);
	}

	/* How many genes? */
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

	t_gene** genes = (t_gene**)malloc(ngenes*sizeof(t_gene));

	/* Get URSes */
	int *urslengths = (int *)malloc(ngenes*sizeof(int));
	char **urses = (char **)malloc(ngenes*sizeof(char*));
	rewind(fu); /* Set the fu file pointer back to the start */
	int this_gene = -1;
	while (  fgets(line, MAXLEN, fu)  ){
		if (line[0] == '>'){
			this_gene += 1;
		}
		else {
			int ii = 0;
			int count = 0;
			while(line[ii] != '\n'){
				count++;
				ii++;
			}
			urslengths[this_gene] = count;
			urses[this_gene] = (char *)malloc(urslengths[this_gene]*sizeof(char));
			for(int ii=0; ii< urslengths[this_gene]; ii++){
				urses[this_gene][ii] = line[ii];
				//printf("(settings 77) urs %d = %s\n", this_gene, urses[this_gene]);
			}
		}
	}

	/* Get PSAMS */
	bool *has_dbd = (bool *)malloc(ngenes*sizeof(bool));
	for(int ii=0; ii<ngenes; ii++){
		has_dbd[ii] = false;
	}
	int *psamlengths = (int *)malloc(ngenes*sizeof(int)); // key = gene ID, value = length of psam
	double **psams = (double **)malloc(ngenes*sizeof(double*)); // key = gene ID, value = 1-d array of doubles
	bool *reg_modes = (bool *)malloc(ngenes*sizeof(bool));
	this_gene = -1;
	int count = 0;
	while ( fgets(line, MAXLEN, fp)  ){
		//printf("(settings 88) line=%s\n", line);

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
			else if (token[0] == 'p'
					&& token[1] == 's'
					&& token[2] == 'a'
					&& token[3] == 'm') {
				token = strtok(NULL, " ");
				this_gene = atoi( token );
				has_dbd[this_gene] = true;
				token = strtok(NULL, " ");
				reg_modes[ this_gene ]= (bool)atoi( token );
				count = 0;
				psams[this_gene] = (double *)malloc(MAX_PSAM_LEN*sizeof(double)); //reset this_psam
				psamlengths[this_gene] = 0;
			}
			else{ // a line with 4 floating-point numbers
				//printf("(settings 103) line=%s token=%s\n", line, token);
				float value = atof( token );
				//printf("(settings 104) this_gene = %d, value=%f count=%d\n",this_gene, value, count);
				psams[this_gene][count] = value;
				//printf("(settings 109) value=%f\n", value);
				count++;
				for (int ii=1; ii<4; ii++){ // ii = state
					value = atof( strtok(NULL, " ") ); //subsequent tokens
					//printf("(settings 112) value=%f\n", value);
					psams[this_gene][count] = value;
					count++;
				}
				psamlengths[this_gene]++;
			}
		}
	}

	/* Build the genes, using the URSes and PSAMs we previously found. */
	for (int ii=0; ii<ngenes; ii++){ // ii = gene
		genes[ii] = make_gene(psamlengths[ii], urslengths[ii]);
		for (int jj=0; jj<urslengths[ii]; jj++){ // jj = site
			genes[ii]->urs[jj] = nt2int( urses[ii][jj] );
		}
		if (has_dbd[ii]){
			for (int jj=0; jj<N_STATES*psamlengths[ii]; jj++){
				genes[ii]->dbd->data[jj] = psams[ii][jj];
			}
			genes[ii]->reg_mode = reg_modes[ii];
		}
	}

	return genes;
}
