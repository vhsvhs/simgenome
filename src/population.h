typedef struct __Population {
	t_genome** genomes;
	int ngenomes;
}t_pop;

t_pop* make_population_random(int popsize, int ngenes, settings *ss);
t_pop* make_population(t_genome* gn, settings *ss);
void free_pop(t_pop* pop);
t_pop* reproduce(t_pop* pop, settings* ss, double* f);
t_genome* mate(t_genome* par1, t_genome* par2);
void print_population(t_pop* pop, settings* ss);

