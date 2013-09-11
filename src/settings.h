typedef struct __Settings {
	int verbosity;
	/* Mutation rates */
	double pwmlenmu;
	int pwmlenmumax;
	double ddgmu;

	char* outdir; /* The directory into which output is written */

	char* psampath;
	char* urspath;
	char* rulepath;

	bool inherit_expression;

	int gen_counter;

	int maxgd; /* The maximum distance, in sites, over which two co-factors can interact. */
	int niid;

	double pe_scalar;

	double growth_rate;
	double decay_rate;

}settings;

settings* make_settings();

