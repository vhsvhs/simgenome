#include "common.h"

settings* make_settings(){
	settings *ss;
	ss = (settings *)malloc(1*sizeof(settings));

	ss->verbosity = DEF_VERBOSITY; //(int)DEF_VERBOSITY;

	ss->psamlenmu = PWMLENMU;
	ss->psamlenmumax = PWMLENMUMAX;
	ss->ddgmu = DDGMU;

	ss->outdir = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	ss->psampath = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	ss->urspath = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	ss->cooppath = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	ss->rulepath = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	ss->poppath = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	ss->load_save_pop = false;

	ss->popsize = POPSIZE;

	ss->inherit_expression = false;
	ss->gen_counter = 0;
	ss->max_gens = MAX_GENS;

	ss->maxgd = MAX_GD;
	ss->niid = NIID;
	ss->maxtime = MAX_TIME;

	ss->do_mutation = 1;
	ss->urs_mu_rate = URSMU; // subs per seq site
	ss->psam_mu_rate = PSAMMU;  // subs per psam site

	ss->pe_scalar = PE_SCALAR;

	ss->growth_rate = GROWTH_RATE;
	ss->decay_rate = DECAY_RATE;

	ss->run_clean = false;

	return ss;
}

void free_settings(settings* ss){
	free(ss->outdir);
	free(ss->psampath);
	free(ss->urspath);
	free(ss->cooppath);
	free(ss->rulepath);
	free(ss->poppath);
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
			{"poppath", required_argument,NULL,	6}, // read one genome, make the population a copy of this genome.


			{"nomu", 		required_argument, 	NULL, 	100}, // disable mutations
			{"psamlenmu", 	required_argument, 	NULL, 	101},
			{"psamlenmumax", required_argument, NULL, 	102},
			{"ddgmu", 		required_argument, 	NULL, 	103},
			{"urs_mu",		required_argument,	NULL,	104},
			{"psam_mu",		required_argument,	NULL,	105},

			{"niid",		required_argument, 	NULL,	200},
			{"maxgen",		required_argument, 	NULL,	201},
			{"startgen",	required_argument,	NULL,	202},
			{"popsize",		required_argument,	NULL,	203},
			{"maxgd",		required_argument,	NULL,	204}, // maximum co-factor distance

			{"run_clean",	no_argument,		NULL, 	300}, // erase previous output files

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
				read_path_from_cli(ss->outdir, true);
				build_output_folders(ss);
				break;
			} // end case 1
			case 2:{
				read_path_from_cli(ss->psampath, false);
				break;
			}
			case 3:{
				read_path_from_cli(ss->urspath, false);
				break;
			}
			case 4:{
				read_path_from_cli(ss->rulepath, false);
				break;
			}
			case 5:{
				read_path_from_cli(ss->rulepath, false);
				break;
			}
			case 6:{
				read_path_from_cli(ss->poppath, false);
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

			case 200:{
				ss->niid = atoi(optarg);
				break;
			}
			case 201:{
				ss->max_gens = atoi(optarg);
				break;
			}
			case 202:{
				ss->gen_counter = atoi(optarg);
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


			case 300:{
				ss->run_clean = true;
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
		strcat( strcat( strcat(qq, "rm -rf "), ss->outdir), "/*");
		system( qq );
		build_output_folders(ss);
	}
}

/* Reads a file path or directory path from the OPTARG variable
 * from the getopts library, and then writes that path to
 * target.
 */
void read_path_from_cli(char* target, bool build) {
	char *tmp;
	tmp = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	if(strlen(optarg) > FILEPATH_LEN_MAX)
	  {
			char choix;
			strcpy (tmp, "\n. The file name'");
			strcat (tmp, optarg);
			strcat (tmp, "' is too long.\n");
			printf("%s",tmp);
			printf("\n. Type any key to exit.\n");
			if(!scanf("%c",&choix)) exit(0);
			exit(0);
	  }
	else if (! Filexists (optarg) && build == false)
	  {
			char choix;
			strcpy (tmp, "\n. Sorry, the file or directory '");
			strcat (tmp, optarg);
			strcat (tmp, "' doesn't exist.\n");
			printf("%s",tmp);
			printf("\n. Type any key to exit.\n");
			if(!scanf("%c",&choix)) exit(0);
			exit(0);
	  }
	else if( ! Filexists (optarg) && build == true){
		/* Make the file, if it doesn't exist */
		printf("\n. I'm making %s\n", optarg);
		mkdir(optarg, 0700);
	}
	else
	  {
			strcpy(target, optarg);
	  }
	free(tmp);
}

void print_splash(){
	printf("\n");
	printf("==========================================\n");
	printf("SIMREG:\n");
	printf("    simulated directed evolution\n");
	printf("    of transcription regulatory circuits.\n");
	printf("\n");
	printf("Written by Victor Hanson-Smith\n");
	printf("    University of California, San Francisco\n");
}

void print_settings(settings *ss){
	if (ss->verbosity > 0){
		printf("\n");
		printf("==========================================\n");
		printf("Current Settings:\n");
		printf(". verbosity: %d\n", ss->verbosity);
		printf(". output directory: %s\n", ss->outdir);
		printf(". PSAMs: %s\n", ss->psampath);
		printf(". URSs: %s\n", ss->urspath);
		printf(". fitness rules: %s\n", ss->rulepath);
		printf("\n");
		printf(". starting generation: %d\n", ss->gen_counter);
		printf(". generation limit: %d\n", ss->max_gens);
		printf(". population size: %d\n", ss->popsize);
		printf("\n");
		if (ss->do_mutation) { printf(". mutations: enabled\n"); }
		else { printf(". mutations: disabled\n"); }
		printf(". URS mu rate: %f\n", ss->urs_mu_rate);
		printf(". PSAM mu rate: %f\n", ss->psam_mu_rate);
		printf(". cofactor mu rate: %f\n", ss->ddgmu);
		printf("\n");
		printf(". n I.I.D. samples: %d\n", ss->niid);
		printf(". pe scalar: %f\n", ss->pe_scalar);
		printf(". max. on rate: %f\n", ss->growth_rate);
		printf(". max. off rate: %f\n", ss->decay_rate);
		printf(". max. cofactor distance: %d sites\n", ss->maxgd);
		printf("==========================================\n");
	}

}
