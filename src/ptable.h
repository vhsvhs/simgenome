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

	int* trancpa; /* Inverse transform of CPR */

}t_ptable;

t_ptable* make_ptable(int M, int D, int L);
void free_ptable( t_ptable* p);
void print_ptable(t_ptable *p);
void ptable_sample(t_ptable* pt, int s, int& reti, int& retj, int& retd);
void build_tran(t_ptable* pt);
void tran_ptable(t_ptable* pt);
void tran_sample(t_ptable* pt, int s, int& reti, int& retj, int& retd);
