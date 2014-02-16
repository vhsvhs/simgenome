#include "common.h"

t_ga* make_ga(){
	t_ga* ga;
	ga = (t_ga*)malloc(sizeof(t_ga));
	return ga;
}

/* Assumes that ga has already been built and
 * that ga->l and ga->pop exist.
 */
void runsim(t_ga* ga, settings* ss){

	if (ss->enable_timelog){
		ss->t_startga = clock();
	}

	/* For each generation */
	int start_gen = ss->start_gen;
	for (int ii = start_gen;
			ii < (start_gen + ss->max_gens);
			ii++){
		ss->gen_counter = ii;

		if (ss->enable_timelog){
			ss->t_startgen = clock();
			ss->t_summakept = 0;
			ss->t_sumfillpt = 0;
			ss->t_sumprob_expr = 0;

		}

		if (ss->verbosity > 2){
			printf("\n\n=================\n");
			printf(" Generation %d\n", ii);
			printf("=================\n");
		}

		/*
		 * Pre-fitness Setup
		 */
		if (ii == 0 || ii == ss->start_gen){ // Only do this for generation zero
			for(int ii = 0; ii < ga->pop->ngenomes; ii++){
				// allocate memory AND init gene expression
				build_lifespan(ga->pop->genomes[ii], ga->l->ntime);
			}
		}
		else{
			for (int gid = 0; gid < ga->pop->ngenomes; gid++){
				// just re-init gene expression (don't re-allocate memory)
				reset_lifespan( ga->pop->genomes[gid]);
			}
		}

		/*
		 * FITNESS
		 */

		/* f is an array of fitness values, one for each individual */
		double* f = (double *)malloc(ga->pop->ngenomes * sizeof(double));
		double* er = (double *)malloc(ga->pop->ngenomes * sizeof(double));
		for (int qq = 0; qq < ga->pop->ngenomes; qq++){
			f[qq] = 0.0;
			er[qq] = 0.0;
		}
		double maxf = 0;
		double minf = 1.0;
		double meanf = 0;
		double medianf = 0;
		double stdf = 0;


		if (ss->enable_timelog){
			ss->t_startf = clock();
		}

		/* Consider each individual */

		for (int gid=0; gid < ga->pop->ngenomes; gid++){
			/*
			 * GET FITNESS of an individual
			 */
			double this_er = 0.0;
			f[gid] = get_fitness(ga->pop->genomes[gid], ga->l, ss, this_er);
			//
			// to-do: get fitness for each regulatory problem, in addition to fitness overall regulatory problems.
			//
			er[gid] = this_er;
		}

		if (ss->enable_timelog){
			ss->t_sumf += clock() - ss->t_startf;
		}


		/* Using the new fitness values, mark elite individuals */
		mark_elite(ga->pop, f, ss);

		printf("==================================\n");
		printf(". Generation %d Fitness:\n", ii);
		for (int gid = 0; gid < ga->pop->ngenomes; gid++){
			printf("\t ID %d, f= %f", gid, f[gid]);
			if (ga->pop->genomes[gid]->is_elite){
				printf(" (elite)");
			}
			printf("\n");
		}
		printf("==================================\n");

		/* Save the fitness stats, and serialize the population */
		log_fitness(f, er, ga->pop->ngenomes, ss);

		/* Save the population to disk. */
		serialize_population(ga->pop, ss);

		/* Reproduce and mutate, but not on the final generation. */
		if (ii < start_gen + ss->max_gens - 1) {
			ga->pop = reproduce( ga->pop, ss, f);
			mutate(ga->pop, ss);
		}

		if (ss->enable_timelog){
			ss->t_stopgen = clock();
			log_timegen(ss, ii);
		}
	}

	if (ss->enable_timelog){
		ss->t_stopga = clock();
	}
}

void free_ga(t_ga* ga){
	free_pop(ga->pop);
	free_landscape(ga->l);
}


bool is_world_consistent(t_ga* ga){
	for (int ii = 0; ii < ga->l->nrulesets; ii++){

		/* For each input */
		for (int jj = 0; jj < ga->l->rulesets[ii]->ninputs; jj++){
			int this_gid = ga->l->rulesets[ii]->inputs[jj]->gid;
			/* Does gene jj exist? */
			if (this_gid >= ga->pop->genomes[0]->ngenes){
				printf("\n. Hmmm, somethng is wrong with your input definitions.");
				printf("\n. Gene %d does not exist.\n", this_gid);
				return false;
			}
		}
		/* for each rule */
		for (int jj = 0; jj < ga->l->rulesets[ii]->nrules; jj++){
			int this_gid = ga->l->rulesets[ii]->rules[jj]->repid;
			if (this_gid >= ga->pop->genomes[0]->ngenes){
				printf("\n. Hmmm, somethng is wrong with your rules.");
				printf("\n. Gene %d does not exist.\n", this_gid);
				return false;
			}
		}
	}
	return true;
}

