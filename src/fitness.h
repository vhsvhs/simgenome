#define MEMO_AFFINITY 1

double get_fitness(t_genome* g, t_landscape* l, settings* ss, double& this_er);
double get_expr_modifier(t_genome *g, int gid, int t, int rid, settings *ss);

#ifdef MEMO_AFFINITY
void fill_prob_table(t_genome *g, int gid, t_ptable *ret, t_afftable *afft, int t, settings *ss);
double prob_expr(t_genome *g, int gid, t_ptable *pt, t_afftable *afft, int t, settings *ss, double* tf_k);
#else
void fill_prob_table(t_genome *g, int gid, t_ptable *ret, int t, settings *ss);
double prob_expr(t_genome *g, int gid, t_ptable *pt, int t, settings *ss, double* tf_k);

#endif

#ifdef PTHREADS
double niid_sample(struct pthread_args* args);
#endif
