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

}settings;

settings* make_settings();

