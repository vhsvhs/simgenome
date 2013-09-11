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

	ss->maxgd = MAX_GD;
	ss->niid = NIID;

	ss->pe_scalar = PE_SCALAR;

	ss->growth_rate = GROWTH_RATE;
	ss->decay_rate = DECAY_RATE;

	return ss;
}
