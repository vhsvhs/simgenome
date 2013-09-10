typedef struct __Gene {
	int id;
	char* name;
	int* urs;
	int urslen;
	bool has_dbd;
	psam* dbd;
	bool is_repressor;

	// to-do: add code to init these arrays
	double* tfcoop; /* tfcoop[gene] = my relative affinity for gene */
	int tfcooplen;
	double* gamma; /* gamma[gene][d] = my cofactor affinity for gene at distance d */
	int gammalen;
}t_gene;

t_gene* make_gene(int psamlen, int urslen);
int* get_random_seq(int len);
void print_urs(int* urs, int urslen);
