typedef struct __Genome {
	int id;
	t_gene** genes; /* An array of pointers to genes */
	int ngenes;
	double* gene_expr; /* gene_expr[gene->id][timeslice] = value ranging from 0.0 to 1.0 */
	bool is_elite;

}t_genome;

t_genome* make_genome(int ngenes);
t_genome* make_genome_default(int ngenes, settings *ss);
t_genome* make_genome_random(int ngenes);
