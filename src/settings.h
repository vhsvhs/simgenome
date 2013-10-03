#include <getopt.h>

typedef struct __Settings {
	int verbosity;
	/* Mutation rates */
	bool do_mutation;
	double urs_mu_rate;
	double psam_mu_rate;
	double psamlenmu;
	int psamlenmumax;
	double ddgmu;

	char* outdir; /* The directory into which output is written */

	char* psampath;
	char* urspath;
	char* cooppath;
	char* rulepath;
	char* poppath;
	bool load_save_pop;

	FILE* file_expr_log;

	bool inherit_expression;

	int gen_counter; /* a counter that gets updated during the run */
	int max_gens; /* the maximum number of generations */

	int maxgd; /* The maximum distance, in sites, over which two co-factors can interact. */
	int niid;
	int popsize;

	int maxtime;
	double elite_proportion;


	/* Random-Init Stuff: */
	bool build_random_population;
	int ngenes;
	int urslen;
	int nreg;

	double pe_scalar;
	double fitness_scalar;

	double growth_rate;
	double decay_rate;

	/* run_clean: if true = remove all files from the output directory
	 * before starting the genetic algorithm.  This prevents
	 * "stale" files from previous runs form polluting
	 * the output of this run.
	 */
	bool run_clean;


}settings;

settings* make_settings();
void free_settings(settings* ss);
void read_cli(int argc, char **argv, settings* ss);
void read_path_from_cli(char* target, bool build);
void print_splash();
void print_settings(settings* ss);
