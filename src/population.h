typedef struct __Population {
	t_genome** genomes;
	int ngenomes;
	double min_fitness; /* The lowest fitness value from the last generation */
	double max_fitness; /* The highest fitness value from the last generation */
}t_pop;

t_pop* make_population_random(int popsize, int ngenes, settings *ss);
t_pop* make_population(int popsize, t_genome* gn, settings *ss);
void free_pop(t_pop* pop);
void reproduce(t_pop* pop, settings* ss);
void print_population(t_pop* pop, settings* ss);
