double get_fitness(t_genome* g, t_landscape* l, settings* ss, double& this_er);
double get_expr_modifier(t_genome *g, int gid, int t, int rid, settings *ss);
void fill_prob_table(t_genome *g, int gid, t_ptable *ret, int t, settings *ss);
double prob_expr(t_genome *g, int gid, t_ptable *pt, int t, settings *ss, double* tf_k);
