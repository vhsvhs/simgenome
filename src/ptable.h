#define CDF_CPA 1

typedef struct __ProbTable {
	double* cpa;
	double* cpt;
	double* cpr;

	int dim1;
	int dim2;
	int dim3;

	int M; /* n transcription factors */
	int L; /* length of upstream regulatory sequence */
	int D; /* maximum interaction distance between two cofactors */

#ifdef CDF_CPA
	double *cpacdf;
#endif

}t_ptable;

t_ptable* make_ptable(int M, int D, int L);
void free_ptable( t_ptable* p);
void print_ptable(t_ptable *p);
void ptable_sample(t_ptable* pt, int s, int& reti, int& retj, int& retd);

typedef struct __AffTable{
	int length;
	double* tf_site_aff; /* Records each TFs affinity to each site in an URS. */
	// tf_site_aff[ tf_i * nsites + site_j ]

}t_afftable;

t_afftable* make_afftable(int ntfs, int nsites);
void free_afftable( t_afftable *t);


