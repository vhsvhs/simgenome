typedef struct __GeneticAlgorithm {
	t_pop* pop;
	t_landscape* l;
}t_ga;

t_ga* make_ga();
void runsim(t_ga* ga, settings* ss);
void get_fitness_stats(double *f, int len, double &max, double &min, double &mean, double &median, double &stddev);

