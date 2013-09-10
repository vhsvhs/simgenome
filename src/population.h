typedef struct __Population {
	t_genome** genomes;
	int ngenomes;
	double min_fitness; /* The lowest fitness value from the last generation */
	double max_fitness; /* The highest fitness value from the last generation */
}t_pop;

t_pop* make_population(int popsize, int ngenes, settings *ss);

void reproduce(t_pop* pop, settings* ss);
