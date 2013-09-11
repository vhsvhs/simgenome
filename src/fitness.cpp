#include "common.h"

/* Returns the fitness of the genome g on the
 * landscape described by l.
 */
double get_fitness(t_genome* g, t_landscape* l, settings* ss){



	/* Build the r vector */
	int* r = (int *)malloc(g->ngenes*sizeof(int));
	for (int ii=0; ii<g->ngenes; ii++){
		if (g->genes[ii]->has_dbd){
			r[ii] = g->genes[ii]->dbd->nsites;
		}
		else {
			r[ii] = -1;
		}
	}

	/* Consider each rule set */
	for (int rid=0; rid<l->nrulesets; rid++) {

		/* Initialize gene expression levels */
		if (!ss->inherit_expression){
			for (int gid=0; gid < g->ngenes; gid++){
				for (int t=0; t<l->ntime; t++){
					g->gene_expr[gid*g->ngenes + t] = MINIMUM_EXPRESSION_LEVEL;
				}
			}
		}

		//printf("fitness 33, l->ntime= %d\n", l->ntime);

		/* For each timeslice */
		for (int t=0; t< l->ntime; t++){

			/* Manually set expression levels
			 * for genes defined in the rule set
			 * where its applicable for timeslice t.
			 */
			for (int ii=0; ii< l->rulesets[rid]->ninputs; ii++){
				t_input* x = l->rulesets[rid]->inputs[ii];
				if (t >= x->start && t < x->start){
					g->gene_expr[ x->gid*g->ngenes + t] = x->expr_level;
				}
			}

			/* Knock-out any genes that have been specified
			 * in the KO array.
			 */

			//TBD


			/* Print a special line for the 0th timeslice. */
			if (t == 0 && ss->verbosity > 1){
				for (int gid=0; gid< g->ngenes; gid++){
					char mark = ' ';
					if (g->genes[gid]->has_dbd && g->genes[gid]->reg_mode == 0){
						mark = 'r';
					}
					else if (g->genes[gid]->has_dbd && g->genes[gid]->reg_mode == 1) {
						mark = 'a';
					}
					printf("r: %d\tgenr: %d\ttime: 0\tID: %d\tgene: %d\t%c\texpr: %f\n",
							rid, ss->gen_counter, g->id, gid, mark, g->gene_expr[gid*g->ngenes]);
				}
			}


			/* For each gene, update is expression level based on
			 * the delta-G of binding
			 */
			for (int gid=0; gid< g->ngenes; gid++){
				/* To-do: if gid is in the knock-out, then actually knock it out here.*/

				// to-do:
				//double pe = get_expr_modifier(g, gid, ss);

				// continue here

			}
		}
	}

	return 0.0;
}
