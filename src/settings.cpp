#include "common.h"

settings* make_settings(){
	settings *ss;
	ss = (settings *)malloc(1*sizeof(settings));

	ss->verbosity = DEF_VERBOSITY; //(int)DEF_VERBOSITY;
	ss->pwmlenmu = PWMLENMU;
	ss->pwmlenmumax = PWMLENMUMAX;
	ss->ddgmu = DDGMU;
	return ss;
}
