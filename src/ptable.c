#include "common.h"

t_afftable* make_afftable(int ntfs, int nsites){
	t_afftable *t;
	t = (t_afftable *)malloc(1*sizeof(t_afftable));
	t->tf_site_aff = (double *)malloc(ntfs*nsites*sizeof(double));
	for (int ii=0; ii < ntfs*nsites; ii++){
		t->tf_site_aff[ii] = 0.0;
	}
	return t;
}

void free_afftable( t_afftable *t){
	free(t->tf_site_aff);
	free(t);
}

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

#ifdef CDF_CPA
	p->cpacdf = (double *)malloc(L*p->dim1*sizeof(double));
	for (int ii=0; ii<L*p->dim1; ii++){
		p->cpacdf[ii] = 0.0;
	}
#endif


	return p;
}

void free_ptable( t_ptable* p){
	free(p->cpa);
	free(p->cpt);
	free(p->cpr);
#ifdef CDF_CPA
	free(p->cpacdf);
#endif
	free(p);
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
	/*
	 * To-do: approximately 97% of the computational runtime is spent
	 * in this method, according to profiling analysis.
	 * How can I make this function faster?
	 */

	double totp = pt->cpr[s]; // total binding energy at this site
	double randp = drand() * totp;

#ifndef CDF_CPA
	double sump = 0.0;
#endif
	reti = 0;
	retj = 0;
	retd = 0;

	for (int ii=s*pt->dim1; ii < s*pt->dim1 + pt->dim1; ii++) {
#ifdef CDF_CPA
		if (pt->cpacdf[ii] > randp){
#else
		sump += pt->cpa[ii];
		if (sump > randp){
#endif
			reti = ii%pt->dim1 / pt->dim2;
			retj = ii%pt->dim2 / pt->dim3;
			retd = ii%pt->dim3;
			return;
		}

	}

	//Depricated code, replaced by above on February 15th, 2014
	/*
	int x1 = s * pt->dim1;
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
	*/
}

