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

	p->trancpa = NULL;

	return p;
}

void free_ptable( t_ptable* p){
	free(p->cpa);
	free(p->cpt);
	free(p->cpr);
	if (p->trancpa != NULL){
		printf("\n. ptable32\n");
		free(p->trancpa);
		printf("\n. ptable34\n");
	}
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
	double totp = pt->cpr[s]; // total binding energy at this site
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

/*
 * Inverse Transform Sampling on the CPA table
 * December 2013: an alternate to repeated CDF sampling.
 *
 */
void build_tran(t_ptable* pt){
	printf("\n. 93: building tran_cpa\n");
	pt->trancpa = (int *)malloc(pt->L*TRANPT_SIZE*3*sizeof(int));
	for (int ii=0; ii<pt->L; ii++){
		for (int jj=0; jj<TRANPT_SIZE*3;jj++){
			pt->trancpa[ii*TRANPT_SIZE*3 + jj] = 0;
		}
	}
	printf("\n. 100: done building tran_cpa\n");
}

void tran_ptable(t_ptable* pt)
{
	printf("\n. 105: entered tran_ptable\n");
	for (int s = 0; s < pt->L; s++){ // for each site
		double totp = pt->cpr[s]; // total P at this site
		int x1 = s*pt->dim1;
		double sump = 0.0;
		int lastindex = 0;
		int index = 0;
		for (int qq = 0; qq < pt->M; qq++){ // qq = TF on the left
			int x2 = qq * pt->dim2;
			for (int rr = 0; rr < pt->M + 1; rr++){ // rr = TF (or empty) on the right
				int x3 = rr * pt->dim3;
				for(int ss = 0; ss < pt->D; ss++){
					sump += pt->cpa[x1 + x2 + x3 + ss];
					index = int((sump / totp) * TRANPT_SIZE) * 3;
					// fill-in any p values we jumped over...
					for (int jj=lastindex+3; jj<index; jj += 3){
						printf("121a: jj = %d\n", jj);
						pt->trancpa[s*TRANPT_SIZE*3 + jj] = pt->trancpa[s*lastindex];
						pt->trancpa[s*TRANPT_SIZE*3 + jj + 1] = pt->trancpa[s*lastindex + 1];
						pt->trancpa[s*TRANPT_SIZE*3 + jj + 2] = pt->trancpa[s*lastindex + 2];
					}
					printf("121b: index = %d\n", index);
					pt->trancpa[s*TRANPT_SIZE*3 + index] = qq;
					pt->trancpa[s*TRANPT_SIZE*3 + index + 1] = rr;
					pt->trancpa[s*TRANPT_SIZE*3 + index + 2] = ss;
					lastindex = index;
				}
			}
		}
		for (int jj=lastindex+3; jj<TRANPT_SIZE*3; jj += 3){
			printf("121c: jj = %d\n", jj);
			pt->trancpa[s*TRANPT_SIZE*3 + jj] = pt->trancpa[s*lastindex];
			pt->trancpa[s*TRANPT_SIZE*3 + jj + 1] = pt->trancpa[s*lastindex + 1];
			pt->trancpa[s*TRANPT_SIZE*3 + jj + 2] = pt->trancpa[s*lastindex + 2];
		}
	}
	printf("\n. 138: done with tran_ptable\n");
	printf("");
}

void tran_sample(t_ptable* pt, int s, int& reti, int& retj, int& retd){
	printf("\n. 142: entered tran_sample\n");
	int randindex = int(drand() * 3*TRANPT_SIZE);
	reti = pt->trancpa[s*3*TRANPT_SIZE + randindex];
	retj = pt->trancpa[s*3*TRANPT_SIZE + randindex + 1];
	retd = pt->trancpa[s*3*TRANPT_SIZE + randindex + 2];
	printf("\n. 147: done tran_sample %d %d %d %d\n",s, reti, retj, retd);
}
