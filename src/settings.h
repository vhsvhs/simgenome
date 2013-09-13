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
	char* rulepath;

	bool inherit_expression;

	int gen_counter; /* a counter that gets updated during the run */
	int max_gens; /* the maximum number of generations */

	int maxgd; /* The maximum distance, in sites, over which two co-factors can interact. */
	int niid;



	double pe_scalar;

	double growth_rate;
	double decay_rate;

}settings;

settings* make_settings();
void read_cli(int argc, char **argv, settings* ss);
void read_path_from_cli(char* target, bool build);
void print_splash();
void print_settings(settings *ss);
