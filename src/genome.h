typedef struct __Genome {
	int id;
	t_gene** genes; /* An array of pointers to genes */
	int ngenes; /* N total genes */
	int ntfs; /* N transcription factors */
	double* gene_expr; /* gene_expr[gene->id * ntimeslices + timeslice] = value ranging from 0.0 to 1.0 */
	int expr_timeslices; //to-do: build gene_expr and expr_timeslices
	bool is_elite;
	int* r; /* Used in calc_prob_tables to temporarily store the lengths of PSAMs */
}t_genome;

t_genome* make_genome(int ngenes, t_gene** ingenes, settings *ss);
t_genome* make_genome_default(int ngenes, settings *ss);
t_genome* make_genome_random(int ngenes);
void init_lifespan(t_genome* g, int t);

void copy_genome(t_genome* to, t_genome* from);
void free_genome(t_genome* gn);

