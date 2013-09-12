#include "common.h"

settings* make_settings(){
	settings *ss;
	ss = (settings *)malloc(1*sizeof(settings));

	ss->verbosity = DEF_VERBOSITY; //(int)DEF_VERBOSITY;
	ss->pwmlenmu = PWMLENMU;
	ss->pwmlenmumax = PWMLENMUMAX;
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
