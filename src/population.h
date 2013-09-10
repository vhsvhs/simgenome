typedef struct __Population {
	t_genome** genomes;
	int ngenomes;
}t_pop;

t_pop* make_population(int popsize, int ngenes, settings *ss);
