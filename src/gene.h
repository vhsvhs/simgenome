typedef struct __Gene {
	int id;
	char* name;
	int* urs;
	int urslen;
	bool has_dbd;
	psam* dbd;
	int reg_mode; /* 0 = repressor, 1 = activator */

	// to-do: add code to init these arrays
	double* tfcoop; /* tfcoop[gene] = my relative affinity for gene */
	int tfcooplen;

	/* gamma[gene*maxgd + d]  = my cofactor affinity
	 * for gene at distance d, given the maximum gd distance
	 * in the settings.
	 */
	double* gamma;
	int gammalen;
}t_gene;

t_gene* make_gene(int psamlen, int urslen);
t_gene* copy_gene(t_gene* from);
void build_coop(t_gene* g, int ntfs, int maxd);
void init_coop(t_gene* g);
double coopfunc(double f, int d);
void calc_gamma(t_gene* g, int ntfs, int maxd);
void free_gene(t_gene* g);
int* get_random_seq(int len);
void print_urs(int* urs, int urslen);
t_gene** read_genes_from_file(settings *ss, int &ngenes);
