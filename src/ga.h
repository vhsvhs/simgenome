typedef struct __GeneticAlgorithm {
	t_pop* pop;
	t_landscape* l;
}t_ga;

t_ga* make_ga();
void runsim(t_ga* ga, settings* ss);
void free_ga(t_ga* ga);
bool is_world_consistent(t_ga* ga);
