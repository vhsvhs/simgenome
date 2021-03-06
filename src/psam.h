typedef struct __PSAM {
	int nstates; // How many characters in the alphabet?
	int nsites; // How many sites in this PSAM?
	double *data;
}psam;

psam* make_psam(int nsites, int nstates);
void free_psam(psam* p);
void copy_psam(psam *to, psam *from );
void shuffle_psam(psam *p);
void rand_init_psam(psam *p);
void print_psam(psam *p);
double get_affinity(psam *p, int *seq, int seqlen, int startsite);
