#include "common.h"

settings* make_settings(){
	settings *ss;
	ss = (settings *)malloc(1*sizeof(settings));

	ss->verbosity = DEF_VERBOSITY; //(int)DEF_VERBOSITY;
	ss->log_occupancy = false;
	ss->log_k = false;
	ss->log_sample_stride = LOG_SAMPLE_STRIDE;

	ss->do_mutation = 1;
	ss->urs_mu_rate = URSMU; // mean subs per seq site
	ss->psam_mu_rate = PSAMMU;  // mean subs per psam site
	ss->mu_stdev = MU_STDEV;

	ss->psamlenmu = PSAMLENMU;
	ss->psamlensd = PSAMLENSD;
	ss->psamlensizemu = PSAMLENSIZEMU;
	ss->psamlensizesd = PSAMLENSIZESD;
	ss->psamlenmumax = PSAMLENMUMAX;

	ss->ddgmu = DDGMU;
	ss->ddgmusize = DDGMUSIZE;
	ss->ddgmusd = DDGMUSD;

	ss->urslenmu = URSLENMU;
	ss->urslensd = URSLENSD;
	ss->urslensizemu = URSLENSIZEMU;
	ss->urslensizesd = URSLENSIZESD;

	ss->outdir = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	ss->psampath = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	ss->urspath = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	ss->cooppath = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	ss->rulepath = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	ss->poppath = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	ss->timelogpath = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	ss->load_save_pop = false;

	ss->popsize = POPSIZE;

	ss->inherit_expression = false;
	ss->start_gen = 0;
	ss->gen_counter = 0;
	ss->max_gens = MAX_GENS;

	ss->maxgd = MAX_GD;
	ss->niid = NIID;
	ss->randseed = time(0);
	ss->maxtime = MAX_TIME;
	ss->elite_proportion = ELITE_PROPORTION;
	ss->no_sex = false;

	/* Variables for random-init: */
	ss->build_random_population = false;
	ss->ngenes = NGENES_DEFAULT;
	ss->nreg = NGENES_DEFAULT/2;
	ss->urslen = 1000;

	ss->pe_scalar = PE_SCALAR;
	ss->fitness_scalar = FITNESS_SCALAR;

	ss->growth_rate = GROWTH_RATE;
	ss->decay_rate = DECAY_RATE;

	ss->run_clean = false;

	ss->enable_timelog = true;
//	ss->use_tran_sampling = false; // Not finished implemeting this feature: use transitive CDF sampling,

	ss->mu_reporters_only = false;
	ss->mu_regulators_only = false;
#ifdef PTHREADS
	ss->n_pthreads = 4;
#endif

	return ss;
}

void free_settings(settings* ss){
	free(ss->outdir);
	free(ss->psampath);
	free(ss->urspath);
	free(ss->cooppath);
	free(ss->rulepath);
	free(ss->poppath);
	free(ss->timelogpath);
	fclose(ss->file_expr_log);
	if (ss->enable_timelog){
		fclose(ss->file_time_log);
	}
}

void read_cli(int argc, char **argv, settings* ss){
	/* This struct is for the getopts library */
	struct option longopts[] =
	{
			{"verbosity", 	required_argument, 	NULL, 	0},
			{"outdir",		required_argument,	NULL,	1},
			{"psampath",	required_argument,	NULL,	2},
			{"urspath",		required_argument,	NULL,	3},
			{"cooppath",	required_argument,	NULL,	4},
			{"rulepath",	required_argument,	NULL,	5},
			{"poppath", 	required_argument,  NULL,	6}, // read one genome, make the population a copy of this genome.

			{"nomu", 		no_argument, 	    NULL, 	100}, // disable mutations
			{"psamlenmu", 	required_argument, 	NULL, 	101},
			{"psamlenmumax", required_argument, NULL, 	102},
			{"coopmu", 		required_argument, 	NULL, 	103},
			{"urs_mu",		required_argument,	NULL,	104},
			{"psam_mu",		required_argument,	NULL,	105},
			{"growth_scalar",required_argument,	NULL,	106}, // at each timeslice, expression delta = growth_rate * pe_scalar
			{"decay_scalar", required_argument,	NULL,	107}, // at each timeslice, expression delta = decay_rate * pe_scalar;
			{"mu_stdev",	 required_argument, NULL,	108},
			{"urslensd",	required_argument,  NULL, 	109},
			{"psamlensd",	required_argument,	NULL,	110},


			{"niid",		required_argument, 	NULL,	200},
			{"maxgen",		required_argument, 	NULL,	201},
			{"startgen",	required_argument,	NULL,	202},
			{"popsize",		required_argument,	NULL,	203},
			{"maxgd",		required_argument,	NULL,	204}, // maximum co-factor distance
			{"elite_prop",	required_argument,	NULL,	205},
			{"no_sex",		no_argument,		NULL,	206}, // disables parental crossover during reproduction.

			{"pe_scalar", 	required_argument,	NULL,	250}, // in fitness.c, pe = (1 / (1+exp(-1*ss->pe_scalar*(sum_act-sum_rep) ) ) );
			{"f_scalar", 	required_argument,	NULL,	251}, // this_fit = exp( ss->fitness_scalar * error);

			{"run_clean",	no_argument,		NULL, 	300}, // erase previous output files, defaults to keep old output.
			{"randseed",	required_argument,	NULL,	301},

			/* For Random-Init of Population: */
			{"randompop",	no_argument,		NULL,	400},
			{"ngenes",		required_argument,	NULL,	401}, // for random initialization only
			{"urslen", 		required_argument,	NULL,	402},
			{"nreg",		required_argument,	NULL,	403},

			{"time",		no_argument,		NULL,	500},
			{"log_occupancy", no_argument,		NULL,	501},
			{"log_k",		no_argument,		NULL,	502},
			{"log_sample_stride", required_argument,NULL,503},
//			{"tran_cdf",	no_argument,		NULL,	600},

			{"mu_reporters_only", no_argument,	NULL,	1001},
			{"mu_regulators_only", no_argument,	NULL,	1002},

#ifdef PTHREADS
			{"n_pthreads",	no_argument,		NULL,	601},
#endif
			{0,0,0,0}
	};

	/* Here I use the getopt library to parse the command-line
	 * arguments.  See getopt documentation for details on the
	 * method getopt_long.
	 */
	int c;
	while((c = getopt_long(argc, argv,
			 //"qi:load_save_popd:m:b:n:t:f:zk:v:c:a:u:ho:s:x:g:l:ep",
			"",
			longopts, NULL)) != -1)
		{
		switch(c)
		{
			case 0:{
				ss->verbosity = atoi(optarg);
				break;
			}
			case 1:{ // output directory
				read_path_from_cli(ss->outdir);
				break;
			} // end case 1
			case 2:{
				read_path_from_cli(ss->psampath);
				break;
			}
			case 3:{
				read_path_from_cli(ss->urspath);
				break;
			}
			case 4:{
				read_path_from_cli(ss->cooppath);
				break;
			}
			case 5:{
				read_path_from_cli(ss->rulepath);
				break;
			}
			case 6:{
				read_path_from_cli(ss->poppath);
				ss->load_save_pop = true;
				break;
			}


			case 100:{
				ss->do_mutation = false;
				break;
			}
			case 101:{
				ss->psamlenmu = atof(optarg);
				break;
			}
			case 102:{
				ss->psamlenmumax = atoi(optarg);
				break;
			}
			case 103:{
				ss->ddgmu = atof(optarg);
				break;
			}
			case 104:{
				ss->urs_mu_rate = atof(optarg);
				break;
			}
			case 105:{
				ss->psam_mu_rate = atof(optarg);
				break;
			}
			case 106:{
				ss->growth_rate = atof(optarg);
				break;
			}
			case 107:{
				ss->decay_rate = atof(optarg);
				break;
			}
			case 108:{
				ss->mu_stdev = atof(optarg);
				break;
			}
			case 109:{
				ss->urslensd = atof(optarg);
				break;
			}
			case 110:{
				ss->psamlensd = atof(optarg);
				break;
			}

			case 200:{
				ss->niid = atoi(optarg);
				break;
			}
			case 201:{
				ss->max_gens = atoi(optarg);
				break;
			}
			case 202:{
				ss->start_gen = atoi(optarg);
				break;
			}
			case 203:{
				ss->popsize = atoi(optarg);
				break;
			}
			case 204:{
				ss->maxgd = atoi(optarg);
				if (ss->maxgd < 1){
					ss->maxgd = 1;
				}
				break;
			}
			case 205:{
				ss->elite_proportion = atof(optarg);
				break;
			}
			case 206:{
				ss->no_sex = true;
				break;
			}

			case 250:{
				ss->pe_scalar = atof(optarg);
				break;
			}
			case 251:{
				ss->fitness_scalar = atof(optarg);
				break;
			}

			case 300:{
				ss->run_clean = true;
				break;
			}
			case 301:{
				ss->randseed = atof(optarg);
				break;
			}


			case 400:{
				ss->build_random_population = true;
				break;
			}
			case 401:{
				ss->ngenes = atoi(optarg);
				break;
			}
			case 402:{
				ss->urslen = atoi(optarg);
				break;
			}
			case 403:{
				ss->nreg = atoi(optarg);
				break;
			}

			case 500:{
				ss->enable_timelog = true;
				break;
			}
			case 501:{
				ss->log_occupancy = true;
				break;
			}
			case 502:{
				ss->log_k = true;
				break;
			}
			case 503:{
				ss->log_sample_stride = atoi(optarg);
				if (ss->log_sample_stride < 1){
					printf("\n. You cannot use a log sample stride with a value less than 1");
				}
				break;
			}

//			case 600:{
//				ss->use_tran_sampling = true;
//				break;
//			}
#ifdef PTHREADS
			case 601:
			{
				ss->n_pthreads = atof(optarg);
				if (ss->n_pthreads < 1){
					printf("\n. You cannot have n_pthreads < 1.\n. Restoring n_pthreads to 1.");
					ss->n_pthreads = 1;
				}
				break;
			}
#endif

			case 1001:
			{
				ss->mu_reporters_only = true;
				break;
			}
			case 1002:
			{
				ss->mu_regulators_only = true;
				break;
			}

			case 1000:
			{
				printf("\n. I found --hello\n");
				break;
			}
			default:
			{
				//usage();
				break;
			}
		} // end of switch

	} // end of while c = getopt


	/* Now we do some post-parsing business */
	if (ss->run_clean == true){
		char *qq = (char*)malloc(FILEPATH_LEN_MAX*sizeof(char));
		if (ss->verbosity > 60){
			printf("\n. I'm wiping the output folder clean.\n");
		}
		strcat( strcat( strcat(qq, "rm -rf "), ss->outdir), "/*");
		system( qq );
	}

	build_output_folders(ss);
}

/* Reads a file path or directory path from the OPTARG variable
 * from the getopts library, and then writes that path to
 * target.
 */
void read_path_from_cli(char* target) {
	char *tmp;
	tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	if(strlen(optarg) > FILEPATH_LEN_MAX) {
			char choix;
			strcpy (tmp, "\n. The file name'");
			strcat (tmp, optarg);
			strcat (tmp, "' is too long.\n");
			printf("%s",tmp);
			printf("\n. Type any key to exit.\n");
			if(!scanf("%c",&choix)) exit(0);
			exit(0);
	  }

// error checking for existence of path object
// should occur outside this function.
//	else if (! Filexists (optarg) ) {
//			char choix;
//			strcpy (tmp, "\n. Sorry, the file or directory '");
//			strcat (tmp, optarg);
//			strcat (tmp, "' doesn't exist.\n");
//			printf("%s",tmp);
//			printf("\n. Type any key to exit.\n");
//			if(!scanf("%c",&choix)) exit(0);
//			exit(0);
//	  }

	strcpy(target, optarg);
	free(tmp);
}

void print_splash(){
	printf("\n");
	printf("==========================================\n");
	printf("SIMREG:\n");
	printf("    simulated directed evolution\n");
	printf("    of transcription regulatory circuits.\n");
	printf("    Version %s\n", __VERSION);
	printf("\n");
	printf("Written by Victor Hanson-Smith\n");
	printf("    University of California, San Francisco\n");
}

void print_settings(settings *ss){
	/* Write settings to STDOUT */
	if (ss->verbosity > 0){
		printf("\n");
		printf("==========================================\n");
		printf("Current Settings:\n");
		printf(". verbosity: %d\n", ss->verbosity);
		printf(". output directory: \t%s\n", ss->outdir);
		printf(". fitness rules: \t%s\n", ss->rulepath);
		if (ss->load_save_pop == true){
			printf(". starting with the saved population: \t%s\n", ss->poppath);
		}
		else if(ss->build_random_population == true){
			printf(". population: \trandom\n");
		}
		else{
			printf(". PSAM definitions: \t%s\n", ss->psampath);
			printf(". URS definitions: \t%s\n", ss->urspath);
		}
		printf("\n");
		printf(". start at generation: \t%d\n", ss->gen_counter);
		printf(". stop at generation: \t%d\n", ss->max_gens);
		printf(". population size: \t%d\n", ss->popsize);
		printf("\n");
		if (ss->do_mutation) {
			if (ss->mu_reporters_only){
				printf(". mutations: enabled only for reporter genes\n");
			}
			else if (ss->mu_regulators_only){
				printf(". mutations: enabled only for regulator genes\n");
			}
			else{
				printf(". mutations: enabled\n");
			}
			printf(". URS mu rate: \t\t%f\n", ss->urs_mu_rate);
			printf(". PSAM mu rate: \t%f\n", ss->psam_mu_rate);
			printf(". URS indel rate: \t%f\n", ss->urslenmu);
			printf(". PSAM indel rate: \t%f\n", ss->psamlenmu);
			printf(". cofactor mu rate: \t%f\n", ss->ddgmu);
			printf(". 'elite' proportion: \t%f\n", ss->elite_proportion);
		}
		else {
			printf(". mutations: disabled\n");
		}
		if (ss->no_sex){
			printf(". meiosis: disabled\n");
		}
		printf("\n");
		printf(". n I.I.D. samples: \t%d\n", ss->niid);
		printf(". random seed: \t\t%d\n", ss->randseed);
		printf(". pe scalar: \t\t%f\n", ss->pe_scalar);
		printf(". max. on rate: \t%f\n", ss->growth_rate);
		printf(". max. off rate: \t%f\n", ss->decay_rate);
		printf(". max. cofactor dist.: \t%d sites\n", ss->maxgd);
		printf("==========================================\n");
	}

	/* Write settings to log file */
	if (ss->verbosity > 0){
		char* p = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
		strcat(
				strcat(p, ss->outdir),
		"/LOGS/settings.txt");

		FILE *fp;
		fp = fopen(p, "w");
		if (fp == NULL) {
		  fprintf(stderr, "Error: can't open output file %s!\n",
				  p);
		  exit(1);
		}
		fprintf(fp,"\n");
		fprintf(fp,"==========================================\n");
		fprintf(fp,"Current Settings:\n");
		fprintf(fp,". verbosity: %d\n", ss->verbosity);
		fprintf(fp,". output directory: %s\n", ss->outdir);
		fprintf(fp,". fitness rules: %s\n", ss->rulepath);
		if (ss->load_save_pop == true){
			fprintf(fp,". starting with the saved population: %s\n", ss->poppath);
		}
		else if(ss->build_random_population == true){
			fprintf(fp,". population: random\n");
		}
		else{
			fprintf(fp,". PSAMs: %s\n", ss->psampath);
			fprintf(fp,". URSs: %s\n", ss->urspath);
		}
		fprintf(fp,"\n");
		fprintf(fp,". start at generation: %d\n", ss->gen_counter);
		fprintf(fp,". stop at generation: %d\n", ss->max_gens);
		fprintf(fp,". population size: %d\n", ss->popsize);
		fprintf(fp,"\n");
		if (ss->do_mutation) {
			if (ss->mu_reporters_only){
				fprintf(fp, ". mutations: enabled only for reporter genes\n");
			}
			else if (ss->mu_regulators_only){
				fprintf(fp, ". mutations: enabled only for regulator genes\n");
			}
			else{
				fprintf(fp, ". mutations: enabled\n");
			}
			fprintf(fp, ". URS mu rate: %f\n", ss->urs_mu_rate);
			fprintf(fp, ". PSAM mu rate: %f\n", ss->psam_mu_rate);
			fprintf(fp, ". URS indel rate: %f\n", ss->urslenmu);
			fprintf(fp, ". PSAM indel rate: %f\n", ss->psamlenmu);
			fprintf(fp, ". cofactor mu rate: %f\n", ss->ddgmu);
			fprintf(fp, ". 'elite' proportion: %f\n", ss->elite_proportion);
		}
		else {
			fprintf(fp,". mutations: disabled\n");
		}
		fprintf(fp,"\n");
		fprintf(fp,". n I.I.D. samples: %d\n", ss->niid);
		fprintf(fp, ". random seed: %d\n", ss->randseed);
		fprintf(fp,". pe scalar: %f\n", ss->pe_scalar);
		fprintf(fp,". max. on rate: %f\n", ss->growth_rate);
		fprintf(fp,". max. off rate: %f\n", ss->decay_rate);
		fprintf(fp,". max. cofactor distance: %d sites\n", ss->maxgd);
		fprintf(fp,"==========================================\n");

		fclose(fp);
		free(p);
	}

}
