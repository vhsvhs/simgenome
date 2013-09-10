typedef struct __Settings {
	int verbosity;
	/* Mutation rates */
	double pwmlenmu;
	int pwmlenmumax;
	double ddgmu;
}settings;

settings* make_settings();
