void build_output_folders(settings* ss);
void write_psams(t_genome *gn, settings *ss, char* outpath);
void log_fitness(double* f, double* er, int len, settings* ss);
void log_occupancy(t_genome *g, int gid, int t, int rid, t_ptable* ptable, settings* ss);
void log_k(t_genome *g, int gid, int t, int rid, double* tf_k, settings* ss);
void log_cofactor(t_genome *g, settings* ss, FILE* fo);
void log_dbds(t_genome *g, settings* ss, FILE* fo);
void log_urs(t_genome *g, settings* ss, FILE* fo);
void log_timegen(settings* ss, int gen);
