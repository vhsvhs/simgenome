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
