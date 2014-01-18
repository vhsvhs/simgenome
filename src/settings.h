#include <getopt.h>
#include <time.h>

typedef struct __Settings {
	int verbosity;

	/* Mutation-related parameters */
	bool do_mutation;
	double urs_mu_rate; // mean mutation rate for URSs
	double psam_mu_rate; // mean mutation rate for PSAMs
	double mu_stdev; // mutation rate standard deviation

	double urslenmu; // how often are URS indels?
	double urslensd;
	double urslensizemu; // mean size of URS indel
	double urslensizesd;

	double psamlenmu; // how often are PSAM indels?
	double psamlensd;
	double psamlensizemu; // mean size of PSAM indel
	double psamlensizesd;
	int psamlenmumax;

	double ddgmu;

	/* Paths, directories, and output-related stuff */
	char* outdir; /* The directory into which output is written */

	char* psampath;
	char* urspath;
	char* cooppath;
	char* rulepath;
	char* poppath;
	bool load_save_pop;

	FILE* file_expr_log;

	/* Enable paternal inheritence? */
	bool inherit_expression;

	int start_gen; /* The starting generation */
	int gen_counter; /* a counter that gets updated during the run */
	int max_gens; /* the maximum number of generations */
	bool no_sex;

	int maxgd; /* The maximum distance, in sites, over which two co-factors can interact. */
	int niid;
	int popsize;
	int randseed; /* The seed for srand() */

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

	/* Tuning-related parameters */
	bool enable_timelog;
	char* timelogpath;
	FILE* file_time_log;
	clock_t t_startmain;
	clock_t t_stopmain;
	clock_t t_startga;
	clock_t t_stopga;
	clock_t t_startgen;
	clock_t t_stopgen;
	clock_t t_startf; // time to calculate fitness for this generation
	int t_sumf;
	clock_t t_startmakept; // time to make P tables
	int t_summakept;
	clock_t t_startfillpt; // time to fille P tables
	int t_sumfillpt;
	clock_t t_startsamplept; // time to sample form the P tables
	int t_sumsamplept;


	/* Debug-related / backdoor features */
	bool use_tran_sampling;

}settings;

settings* make_settings();
void free_settings(settings* ss);
void read_cli(int argc, char **argv, settings* ss);
void read_path_from_cli(char* target);
void print_splash();
void print_settings(settings* ss);
