#include "common.h"

settings* make_settings(){
	settings *ss;
	ss = (settings *)malloc(1*sizeof(settings));

	ss->verbosity = DEF_VERBOSITY; //(int)DEF_VERBOSITY;
	ss->psamlenmu = PWMLENMU;
	ss->psamlenmumax = PWMLENMUMAX;
	ss->ddgmu = DDGMU;

	ss->outdir = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	ss->inherit_expression = false;
	ss->gen_counter = 0;
	ss->max_gens = MAX_GENS;

	ss->maxgd = MAX_GD;
	ss->niid = NIID;

	ss->do_mutation = true;
	ss->urs_mu_rate = URSMU; // subs per seq site
	ss->psam_mu_rate = PSAMMU;  // subs per psam site

	ss->pe_scalar = PE_SCALAR;

	ss->growth_rate = GROWTH_RATE;
	ss->decay_rate = DECAY_RATE;

	return ss;
}

void read_cli(int argc, char **argv, settings* ss){
	/* This struct is for the getopts library */
	struct option longopts[] =
	{
			{"verbosity", 	required_argument, 	NULL, 	0},
			{"outdir",		required_argument,	NULL,	1},
			{"psampath",	required_argument,	NULL,	2},
			{"urspath",		required_argument,	NULL,	3},
			{"rulepath",	required_argument,	NULL,	4},


			{"do_mu", 		no_argument, 		NULL, 	100},
			{"pwmlenmu", 	required_argument, 	NULL, 	101},
			{"pwnlenmumax", required_argument, 	NULL, 	102},
			{"ddgmu", 		required_argument, 	NULL, 	103},
			{"urs_mu",		required_argument,	NULL,	104},
			{"psam_mu",		required_argument,	NULL,	105},



			{0,0,0,0}
	};

	/* Here I use the getopt library to parse the command-line
	 * arguments.  See getopt documentation for details on the
	 * method getopt_long.
	 */
	int c;
	while((c = getopt_long(argc, argv,
			 //"qi:d:m:b:n:t:f:zk:v:c:a:u:ho:s:x:g:l:ep",
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
				read_path_from_cli(ss->rulepath);
				break;
			}


			case 100:{
				ss->do_mutation = true;
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
}

/* Reads a file path or directory path from the OPTARG variable
 * from the getopts library, and then writes that path to
 * target.
 */
void read_path_from_cli(char* target) {
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
	else if (! Filexists (optarg))
	  {
			char choix;
			strcpy (tmp, "\n. Sorry, the directory '");
			strcat (tmp, optarg);
			strcat (tmp, "' doesn't exist.\n");
			printf("%s",tmp);
			printf("\n. Type any key to exit.\n");
			if(!scanf("%c",&choix)) exit(0);
			exit(0);
	  }
	else
	  {
			strcpy(target, optarg);
			printf("\n. I found a valid path: %s\n", target);
	  }
	free(tmp);
}



void print_settings(settings *ss){
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
	if (ss->do_mutation) { printf(". mutations: enabled\n"); }
	else { printf(". mutations: disabled\n"); }
	printf(". URS mu rate: %f\n", ss->urs_mu_rate);
	printf(". PSAM mu rate: %f\n", ss->psam_mu_rate);
	printf(". cofactor mu rate: %f\n", ss->ddgmu);
	// to-do: indels
	printf("\n");
	printf(". n I.I.D. samples: %d\n", ss->niid);
	printf(". pe scalar: %f\n", ss->pe_scalar);
	printf(". max. on rate: %f\n", ss->growth_rate);
	printf(". max. off rate: %f\n", ss->decay_rate);
	printf("==========================================\n");

}
