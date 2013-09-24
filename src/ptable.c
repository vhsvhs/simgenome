#include "common.h"

t_ptable* make_ptable(int M, int D, int L) {
	t_ptable *p;
	p = (t_ptable *)malloc(1*sizeof(t_ptable));
	p->M = M;
	p->D = D;
	p->L = L;
	p->dim1 = M*(M+1)*D;
	p->dim2 = (M+1)*D;
	p->dim3 = D;
	p->cpa = (double *)malloc(L*p->dim1*sizeof(double));
	for (int ii=0; ii<L*p->dim1; ii++){
		p->cpa[ii] = 0.0;
	}
	p->cpt = (double *)malloc(M*L*sizeof(double));
	for (int ii=0; ii<M*L; ii++){
		p->cpt[ii] = 0.0;
	}
	p->cpr = (double *)malloc(L*sizeof(double));
	for (int ii=0; ii<L; ii++){
		p->cpr[ii] = 0.0;
	}
	return p;
}

void free_ptable( t_ptable* p){
	free(p->cpa);
	free(p->cpt);
	free(p->cpr);
}

void print_ptable(t_ptable *p){
	printf("\n");
	printf("CPT\nsite\tTR\tvalue\n");
	for (int ii=0; ii<p->L; ii++){
		printf("Site %d ", ii);
		for (int jj=0; jj<p->M; jj++){
			printf(" tf %d cpt=%f\n", jj, p->cpt[ii*p->M + jj]);
		}
		printf("\n");
	}
	printf("\n");
	printf("CPR\nsite\tvalue\n");
	for (int ii=0; ii<p->L; ii++){
		printf("Site %d cpr=%f\n", ii, p->cpr[ii]);
	}
}

/* Samples from the ptable's CDF at site s.
 * The resultant i, j, and d values
 * will be located in reti, retj, and retd, respectively.
 */
void ptable_sample(t_ptable* pt, int s, int& reti, int& retj, int& retd){
	/* This following block is essentially a CDF sampler */
	double totp = pt->cpr[s]; // total binding energy at this site
	//double randp = (float)rand()/(float)RAND_MAX * totp;
	double randp = drand() * totp;
	int x1 = s*pt->dim1;
	double sump = 0.0;
	reti = 0;
	retj = 0;
	retd = 0;
	for (int qq = 0; qq < pt->M; qq++){ // qq = TF on the left
		int x2 = qq * pt->dim2;
		for (int rr = 0; rr < pt->M + 1; rr++){ // rr = TF (or empty) on the right
			int x3 = rr * pt->dim3;
			for(int ss = 0; ss < pt->D; ss++){ // ss = distance between TFs
				sump += pt->cpa[x1 + x2 + x3 + ss];
				if (sump > randp){
					reti = qq;
					retj = rr;
					retd = ss;
					return;
				}
			}
		}
	}
}
