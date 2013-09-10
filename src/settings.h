typedef struct __Settings {
	int verbosity;
	/* Mutation rates */
	double pwmlenmu;
	int pwmlenmumax;
	double ddgmu;

	char* outdir;

}settings;

settings* make_settings();
