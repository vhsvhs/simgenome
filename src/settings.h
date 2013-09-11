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

}settings;

settings* make_settings();

