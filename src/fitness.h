double get_fitness(t_genome* g, t_landscape* l, settings* ss);
double get_expr_modifier(t_genome *g, int gid, int t, settings *ss);
void fill_prob_table(t_genome *g, int gid, t_ptable *ret, int t, settings *ss);
double prob_expr(t_genome *g, int gid, t_ptable *pt, int t, settings *ss);
