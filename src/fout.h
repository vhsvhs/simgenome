void build_output_folders(settings* ss);
void write_psams(t_genome *gn, settings *ss, char* outpath);
void log_fitness(double* f, int len, settings* ss);
void log_occupancy(t_genome *g, int t, t_ptable* ptable, settings* ss);
void log_cofactor(t_genome *g, settings* ss);
void log_dbds(t_genome *g, settings* ss);
