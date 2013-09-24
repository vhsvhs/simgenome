#include "common.h"

/* Builds a rule structure */
t_rule* make_rule(int timepoint, int repid, double expr, int ruletype, double weight){
	t_rule* r;
	r = (t_rule*)malloc(1*sizeof(t_rule));
	r->timepoint = timepoint;
	r->repid = repid;
	r->expr_level = expr;
	r->rule_type = ruletype;
	r->weight = weight;
	return r;
}

/* Builds an gene input structure */
t_input* make_input(int start, int stop, int gid, double expr){
	t_input* i;
	i = (t_input*)malloc(1*sizeof(t_input));
	i->start = start;
	i->stop = stop;
	i->gid = gid;
	i->expr_level = expr;
	return i;
}

t_ruleset* make_ruleset(int id, int nrules, int ninputs){
	t_ruleset* rs;
	rs = (t_ruleset*)malloc(1*sizeof(t_ruleset));
	rs->id = id;
	rs->nrules = nrules;
	rs->ninputs = ninputs;
	rs->rules = (t_rule**)malloc(nrules*sizeof(t_rule) );
	rs->inputs = (t_input**)malloc(ninputs*sizeof(t_input) );
	return rs;
}

/* Builds and initializes an array of t_ruleset objects,
 * using the rules and inputs in the file at ss->rulepath.
 *
 * Upon completion, ret_n will hold the number of t_ruleset object
 * that were built into the return value.
 */

t_ruleset** read_rulesets_from_file(settings* ss, int &ret_n, int &ntime){
	FILE *fr; /* File for psam specs */
	fr = fopen(ss->rulepath,"r");
	if (fr == NULL) {
	  fprintf(stderr, "Error: can't open fitness rule file %s!\n",
			  ss->rulepath);
	  exit(1);
	}

	//printf("\n. landscape 36\n");

	/* Pass 1: count the number of rulesets */
	int maxrs = -1;
	int nrs = 0; // n rulesets
	char line[MAXLEN];
	while ( fgets(line, MAXLEN, fr) ){
		//printf("landscape 43: line:%s\n", line);
		const char* tokens[MAX_TOKENS] = {};
		tokens[0] = strtok(line, " ");
		if (tokens[0][0] == '#'){
			continue;
		}
		else if (tokens[0][0] == '\n') {
			continue;
		}
		else if (tokens[0]){ // Does this token contain content?
			//printf("landscape 47: tokens[0]:%s\n", tokens[0]);
			int this_rs = atoi( strtok(NULL, " ") );
			if (this_rs > maxrs){
				maxrs = this_rs;
				nrs++;
			}
		}
	}

	//printf("\n. landscape 54 maxrs= %d nrs= %d\n", maxrs, nrs);
	//exit(1);
	/* Pass 2: count the number of rules and inputs for each ruleset */
	int *nrules = (int *)malloc(nrs*sizeof(int));
	int *ninputs = (int *)malloc(nrs*sizeof(int));
	for (int ii=0; ii<nrs; ii++){
		nrules[ii] = 0;
		ninputs[ii] = 0;
	}
	//printf("\n. landscape 68 nrules[nrs-1]= %d ninputs[nrs-1]= %d\n", nrules[nrs-1], ninputs[nrs-1]);
	rewind(fr); // return fr to the start of the file.
	while ( fgets(line, MAXLEN, fr) ){
		const char* tokens[MAX_TOKENS] = {};
		tokens[0] = strtok(line, " ");
		if (tokens[0][0] == '#'){
			continue;
		}
		if (tokens[0][0] == '\n') {
			continue;
		}
		else if (tokens[0]){
			int this_rs = atoi( strtok(NULL, " " ) );
			int geneid = atoi( strtok(NULL, " ") );
			if (tokens[0][0] == 'R' &&
					tokens[0][1] == 'U' &&
					tokens[0][2] == 'L' &&
					tokens[0][3] == 'E') {
				nrules[this_rs] += 1;
			}
			else if (tokens[0][0] == 'I' &&
					tokens[0][1] == 'N' &&
					tokens[0][2] == 'P' &&
					tokens[0][3] == 'U' &&
					tokens[0][4] == 'T') {
				ninputs[this_rs] += 1;
			}
		}
	}

	//printf("\n. landscape 96 nrules[nrs-1]= %d ninputs[nrs-1]= %d\n", nrules[nrs-1], ninputs[nrs-1]);

	/* Final pass: read and save the rules and inputs */

	// An array of pointers to t_ruleset objects:
	t_ruleset** rulesets;
	rulesets = (t_ruleset**)malloc(nrs*sizeof(t_ruleset));
	// key = ruleset ID, value = the number of inputs for this ruleset
	int* ruleset_countinput;
	ruleset_countinput = (int *)malloc(nrs*sizeof(int));
	// key = ruleset ID, value = the number of rules for this ruleset
	int* ruleset_countrule;
	ruleset_countrule = (int *)malloc(nrs*sizeof(int));
	/* For each rule set... */
	for(int qq = 0; qq < nrs; qq++) {
		t_ruleset *trs;
		trs = make_ruleset(qq, nrules[qq], ninputs[qq]);
		rulesets[qq] = trs;

		// initializing these to zero, for use later in this function
		ruleset_countinput[qq] = 0;
		ruleset_countrule[qq] = 0;
	}
	//printf("\n. landscape 126\n");
	rewind(fr);

	/* . . . and then finally we build and parse the t_ruleset object:
	 */
	while ( fgets(line, MAXLEN, fr) ){
		const char* tokens[MAX_TOKENS] = {};
		tokens[0] = strtok(line, " "); //tokens 0
		if (tokens[0][0] == '#'){
			continue;
		}
		if (tokens[0][0] == '\n') {
			continue;
		}
		else if (tokens[0]){

			// ruleset ID:
			int this_rs = atoi( strtok(NULL, " " ) ); // tok 1

			// gene ID:
			int geneid = atoi( strtok(NULL, " ") ); // tok 2

			// RULE
			if (tokens[0][0] == 'R' &&
					tokens[0][1] == 'U' &&
					tokens[0][2] == 'L' &&
					tokens[0][3] == 'E') {
				int timepoint = atoi( strtok(NULL, " ") ); // tokens 3
				if (timepoint > ntime){
					ntime = timepoint;
				}
				double expr = atof( strtok(NULL, " ") ); // tokens 4
				int ruletype = atoi( strtok(NULL, " ") ); //tokens 5
				double weight = atof( strtok(NULL, " ") ); // tokens 6
				rulesets[this_rs]->rules[ruleset_countrule[this_rs]] = make_rule(timepoint, geneid, expr, ruletype, weight);
				if (ss->verbosity > 3){
					if (ruletype == 0){
						printf("   + Fitness Optimum: time %d, gene %d, expr > %f [weight = %f]\n",
							timepoint, geneid, expr, weight);
					}
					else
					{
						printf("   + Fitness Optimum: time %d, gene %d, expr < %f [weight = %f]\n",
							timepoint, geneid, expr, weight);
					}
				}
				ruleset_countrule[this_rs]++;
			}
			// INPUT
			else if (tokens[0][0] == 'I' &&
					tokens[0][1] == 'N' &&
					tokens[0][2] == 'P' &&
					tokens[0][3] == 'U' &&
					tokens[0][4] == 'T') {
				int timestart = atoi(strtok(NULL, " ") ); // tokens 1
				int timestop = atoi(strtok(NULL, " ") ); // tokens 2
				if (timestop > ntime){
					ntime = timestop;
				}
				double expr = atof(strtok(NULL, " ") ); // tokens 3
				rulesets[this_rs]->inputs[ ruleset_countinput[this_rs] ] = make_input(timestart, timestop, geneid, expr);
				ruleset_countinput[this_rs]++;
			}
		}
	}
	ntime += 1;// to offset for 0.
	ss->maxtime = ntime;
	ret_n = nrs;
	return rulesets;
}


void free_landscape(t_landscape* l){
	for (int ii = 0; ii < l->nrulesets; ii++){
		for (int jj = 0; jj < l->rulesets[ii]->ninputs; jj++){
			free( l->rulesets[ii]->inputs[jj] );
		}
		for (int jj = 0; jj < l->rulesets[ii]->nrules; jj++){
			free( l->rulesets[ii]->rules[jj] );
		}
		free( l->rulesets[ii]->inputs );
		free( l->rulesets[ii]->rules );
		free( l->rulesets[ii] );
	}
	free( l->rulesets );
	free(l);

}
