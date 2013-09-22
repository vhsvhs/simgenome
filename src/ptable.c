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
