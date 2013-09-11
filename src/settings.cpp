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

	return ss;
}
