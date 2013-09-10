typedef struct __Settings {
	int verbosity;
	/* Mutation rates */
	double pwmlenmu;
	int pwmlenmumax;
	double ddgmu;

	char* outdir; /* The directory into which output is written */

	char* psampath;
	char* urspath;

}settings;

settings* make_settings();
t_gene** read_genes_from_file(settings *ss, int &ngenes);
