#include "common.h"

/* Returns the fitness of the genome g on the
 * landscape described by l.
 */
double get_fitness(t_genome* g, t_landscape* l, settings* ss, double& this_er){

	/* Build the r vector */
	g->r = (int *)malloc(g->ngenes*sizeof(int));
	for (int ii=0; ii<g->ngenes; ii++){
		if (g->genes[ii]->has_dbd){
			g->r[ii] = g->genes[ii]->dbd->nsites;
		}
		else {
			g->r[ii] = -1;
		}
	}

	double my_fit = 0.0;

	/* Consider every rule set */
	for (int rid=0; rid < l->nrulesets; rid++) {

		/* Initialize gene expression levels */
		if (!ss->inherit_expression){
			for (int gid=0; gid < g->ngenes; gid++){
				for (int t=0; t<l->ntime; t++){
					g->gene_expr[gid*g->expr_timeslices + t] = MINIMUM_EXPRESSION_LEVEL;
				}
			}
		}

		/* For each timeslice */
		for (int t=0; t < l->ntime; t++){

			//printf("\n. TIME %d\n", t);

			/* Manually set expression levels for genes defined in the rule set
			 * where its applicable for timeslice t.
			 */
			for (int ii=0; ii< l->rulesets[rid]->ninputs; ii++){
				t_input* x = l->rulesets[rid]->inputs[ii];
				//printf("\n 46 %d %d %d\n", x->gid, x->start, x->stop);
				if (t >= x->start && t <= x->stop){
					g->gene_expr[ x->gid * g->expr_timeslices + t] = x->expr_level;
				}
			}

			/* To-do: Knock-out any genes that have been specified
			 * in the KO array.
			 */

			//TBD


			/* If t is not the last timeslice,
			 * then update expression.
			 */

			if (t < l->ntime-1)
			{
				/* For each gene, update its expression level based on
				 * the delta-G of binding
				 */
				for (int gid=0; gid< g->ngenes; gid++){
					/* To-do: if gid is in the knock-out, then actually knock it out here.*/
					double pe = get_expr_modifier(g, gid, t, rid, ss);

					//printf("\n. fitness.c 70 for gene %d, pe = %f",gid,  pe);

					if (pe > 0.0){
						pe = ss->growth_rate * pe;
					}
					else if (pe < 0.0){
						pe = ss->decay_rate * pe;

					}

					//printf("\n. fitness.c 70a gene %d, pe = %f\n",gid,  pe);

					double new_expr_level = g->gene_expr[gid * g->expr_timeslices + t] + pe;
					if (new_expr_level < MINIMUM_EXPRESSION_LEVEL){
						new_expr_level = MINIMUM_EXPRESSION_LEVEL;
					}
					else if (new_expr_level > MAXIMUM_EXPRESSION_LEVEL)
					{
						new_expr_level = MAXIMUM_EXPRESSION_LEVEL;
					}

					g->gene_expr[gid * g->expr_timeslices + t+1] = new_expr_level;

					//printf("\n. fitness.c 70b for gene %d, pe = %f, new_expr_level= %f\n",gid,  pe, new_expr_level);
				}
			}

			if (ss->verbosity > 2){
				for (int gid=0; gid< g->ngenes; gid++){
					char mark = ' ';
					if (g->genes[gid]->has_dbd && g->genes[gid]->reg_mode == 0){
						mark = 'r';
					}
					else if (g->genes[gid]->has_dbd && g->genes[gid]->reg_mode == 1) {
						mark = 'a';
					}

					double delta;

					/* Note: yes I realize this is horrible code style to repeat the same string construction (with some variation)
					 * four times.  Consolidating this code block is on my to-do list.
					 */
					if (t > 0){
						delta = g->gene_expr[gid * g->expr_timeslices + t] - g->gene_expr[gid * g->expr_timeslices + t - 1];
						printf("r: %d\tgenr: %d\ttime: %d\tID: %d\tgene: %d\t%c\texpr: %f\tdelta: %f\n",
							rid, ss->gen_counter, t, g->id, gid, mark, g->gene_expr[gid*g->expr_timeslices + t], delta);
						fprintf(ss->file_expr_log, "r: %d\tgenr: %d\ttime: %d\tID: %d\tgene: %d\t%c\texpr: %f\tdelta: %f\n",
							rid, ss->gen_counter, t, g->id, gid, mark, g->gene_expr[gid*g->expr_timeslices + t], delta);
					}
					else{
						printf("r: %d\tgenr: %d\ttime: %d\tID: %d\tgene: %d\t%c\texpr: %f\n",
							rid, ss->gen_counter, t, g->id, gid, mark, g->gene_expr[gid*g->expr_timeslices + t]);
						fprintf(ss->file_expr_log, "r: %d\tgenr: %d\ttime: %d\tID: %d\tgene: %d\t%c\texpr: %f\n",
							rid, ss->gen_counter, t, g->id, gid, mark, g->gene_expr[gid*g->expr_timeslices + t]);
					}

				}
			}

		} // end for t

		/* Now score this genome's gene expression vs. this regulatory problem */
		double sum_of_wt = 0.0; // sum of rule weights
		for (int rr = 0; rr < l->rulesets[rid]->nrules; rr++){
			sum_of_wt += l->rulesets[rid]->rules[rr]->weight;
		}

		/* For each rule in this regulatory problem */
		for (int rr = 0; rr < l->rulesets[rid]->nrules; rr++){
			t_rule* rul = l->rulesets[rid]->rules[rr];
			double obs_expr = g->gene_expr[ rul->repid
			                                * g->expr_timeslices
			                                + rul->timepoint];
			double rid_error = 0.0; // we'll incrementally add to this error counter

			/* Rule type 0: expression must be greater than the rule. */
			if (rul->rule_type == 0){
				if (obs_expr < rul->expr_level){
					//error += fabs( obs_expr - rul->expr_level) / (MAXIMUM_EXPRESSION_LEVEL-MINIMUM_EXPRESSION_LEVEL);
					// new error function: January 2014:
					rid_error += rul->expr_level / obs_expr - 1.0;
				}
			}

			/* Rule type 1: expression must be less than the rule. */
			if (rul->rule_type == 1){
				if (obs_expr > rul->expr_level){
					//error += fabs( obs_expr - rul->expr_level) / (MAXIMUM_EXPRESSION_LEVEL-MINIMUM_EXPRESSION_LEVEL);
					// new error function: January 2014:
					rid_error += obs_expr / rul->expr_level - 1.0;
				}
			}
			//printf("\n. debug fitness 162: rid %d sumriderr %f ruletype %d rulexpr %f obsexpr %f\n",
			//		rr, rid_error, rul->rule_type, rul->expr_level, obs_expr);
			this_er += rid_error * (rul->weight / sum_of_wt);
		} // end for each rule
	} // end for ruleset
	my_fit = exp(ss->fitness_scalar * this_er);

	//my_fit /= l->nrulesets;

	free(g->r);

	return my_fit;
}

/* Returns a floating-point value, the expression level of gene gid,
 * given the expression levels in g->gene_expr,
 * and the relative DNA-binding affinities of all
 * transcription factors in g.
 */
double get_expr_modifier(t_genome *g, int gid, int t, int rid, settings *ss){

	/*
	 * Make the P table
	 */
	if (ss->enable_timelog){
		ss->t_startmakept = clock();
	}
	t_ptable *ptable = make_ptable( g->ntfs, ss->maxgd, g->genes[gid]->urslen);

	if (ss->enable_timelog){
		ss->t_summakept += clock() - ss->t_startmakept;
		ss->t_startfillpt = clock();
	}

#ifdef MEMO_AFFINITY
	t_afftable *afft = make_afftable(g->ntfs, g->genes[gid]->urslen);
#endif
	/*
	 * Fill the P table with probability values
	 */
#ifdef MEMO_AFFINITY
	fill_prob_table(g, gid, ptable, afft, t, ss);
#else
	fill_prob_table(g, gid, ptable, t, ss);
#endif
	if (ss->enable_timelog){
		ss->t_sumfillpt += clock() - ss->t_startfillpt;
		ss->t_startprob_expr = clock();
	}

	/*
	 * tf_k records the binding energy of each TF for gene 'gid'
	 */
	double* tf_k = (double *)malloc(g->ntfs*sizeof(double));
	for (int ii=0; ii<g->ntfs; ii++){
		tf_k[ii] = 0.0;
	}

	/*
	 * Sample from the P table's CDF
	 */
#ifdef MEMO_AFFINITY
	double pe = prob_expr(g, gid, ptable, afft, t, ss, tf_k);
#else
	double pe = prob_expr(g, gid, ptable, t, ss, tf_k);
#endif
	pe = pe - 0.5;

	if (ss->enable_timelog){
		ss->t_sumprob_expr += clock() - ss->t_startprob_expr;
	}

	/*
	 * Log occupancy and/or binding energies, if logging is enabled.
	 */
	if (ss->log_occupancy){
		if (ss->gen_counter%ss->log_sample_stride == 0){
			log_occupancy(g, gid, t, rid, ptable, ss);
		}
	}
	if (ss->log_k){
		if (ss->gen_counter%ss->log_sample_stride == 0){
			log_k(g, gid, t, rid, tf_k, ss);
		}
	}

	/*
	 * Or, log only occupancies after each Nth generation.
	 * This option can consume much less disk space.
	 */
//	else if (ss->verbosity > 3){
//		if (ss->gen_counter%CONFIG_SAMPLE_STRIDE == 0){
//			log_occupancy(g, gid, t, rid, ptable, ss);
//		}
//	}

#ifdef MEMO_AFFINITY
	free_afftable(afft);
#endif
	free(tf_k);
	free_ptable(ptable);
	return pe;
}

/* Uses a dynamic algorithm to fill values into all cells
 * of the ptable.  The algorithm works in place, and upon
 * completion, the new values will be written into ret.
 */
#ifdef MEMO_AFFINITY
void fill_prob_table(t_genome *g, int gid, t_ptable *ret, t_afftable *afft, int t, settings *ss){
#else
void fill_prob_table(t_genome *g, int gid, t_ptable *ret, int t, settings *ss){
#endif
	for (int xx = 0; xx < ret->L; xx++){ // for each site in the upstream regulatory sequence

		//double debug1 = g->genes[1]->gamma[1*ss->maxgd];
		//printf("\n. debug 244: %f\n", debug1);


		double sum_cpr = 0.0;
		for (int ii = 0; ii < g->ntfs; ii++){ // for each transcription factor

			double aff = get_affinity(g->genes[ii]->dbd, g->genes[gid]->urs, g->genes[gid]->urslen, xx);
#ifdef MEMO_AFFINITY
			afft->tf_site_aff[ii*g->genes[gid]->urslen + xx] = aff;
#endif
			double sum_cpt = 0.0; // the sum of cpt values over this jj.
			for (int jj = 0; jj < g->ntfs + 1; jj++){ // for co-factors, +1 is the empty slot, i.e. no cofactor.
				for (int dd = 0; dd < ret->D; dd++){ // for possible distance between TF and co-factor

					double value = 0.0;
					if (ret->L-xx < g->r[ii]){
						/* Case 1: TF i cannot fit here. */
						value = 0.0;
						ret->cpa[xx*ret->dim1 + ii*ret->dim2 + jj*ret->dim3 + dd] = value;
						if (ss->verbosity > 100){
							printf("site %d case 1: tfs %d %d distance %d p= %f\n",
									xx, ii, jj, dd, value);
						}
					}
					if (jj < g->ntfs) // jj is real, not the empty slot.
					{
						if (dd < MIN_TF_SEPARATION){
							/* Case 4: the distance between ii and jj is too small. */
							value = 0.0;
							ret->cpa[xx*ret->dim1 + ii*ret->dim2 + jj*ret->dim3 + dd] = value;
							if (ss->verbosity > 100){
								printf("site %d case 4: tfs %d %d distance %d p= %f\n",
										xx, ii, jj, dd, value);
							}
						}
						else if (ret->L - g->r[ii] - dd - g->r[jj] < 0){
							/* Case 3: ii and jj cannot both fit.
							 *
							 */
							value = 0.0;
							ret->cpa[xx*ret->dim1 + ii*ret->dim2 + jj*ret->dim3 + dd] = value;
							if (ss->verbosity > 100){
								printf("site %d case 3: tfs %d %d distance %d p= %f\n",
										xx, ii, jj, dd, value);
							}
						}
						else{
							/* Case 5: ii and jj can both fit here. */
							value = g->gene_expr[ii*g->expr_timeslices + t]
							                            * aff
							                            * g->genes[ii]->gamma[jj*ss->maxgd + dd];
							ret->cpa[xx*ret->dim1 + ii*ret->dim2 + jj*ret->dim3 + dd] = value;
							sum_cpt += value;
							sum_cpr += value;
							if (ss->verbosity > 100){
								printf("site %d case 5: tfs %d %d distance %d p= %f\n",
										xx, ii, jj, dd, value);
							}
						}
					}
					else if (jj == g->ntfs){
						/* Case 6: no cofactor */
						value = g->gene_expr[ii*g->expr_timeslices + t] * aff;
						ret->cpa[xx*ret->dim1 + ii*ret->dim2 + jj*ret->dim3 + dd] = value;
						sum_cpt += value;
						sum_cpr += value;
						if (ss->verbosity > 100){
							printf("site %d case 6: tfs %d %d distance %d p= %f\n",
									xx, ii, jj, dd, value);
						}
					}
				}
			} // end for jj

			ret->cpt[xx * ret->M + ii] = sum_cpt;
			//printf("fitness 249 %d %d %f\n", xx, ii, ret->cpt[xx * ret->M + ii]);

		} // end for ii
		ret->cpr[xx] = sum_cpr;
	}

#ifdef CDF_CPA
	// Precompute the C.D.F. for the cpa distribution.
	for (int ii=0; ii < ret->L*ret->dim1; ii++){
		if (ii == 0){ ret->cpacdf[ii] = ret->cpa[ii]; }
		else if (ii%ret->dim1 == 0){ ret->cpacdf[ii] = ret->cpa[ii];}
		else {
			ret->cpacdf[ii] = ret->cpa[ii] + ret->cpacdf[ii-1];
		}
	}
#endif

}

#ifdef PTHREADS

struct pthread_args {
	int pthread_id;
	double* return_pes;
	t_genome *g;
	int gid;
	t_ptable *pt;
	t_afftable* afft;
	int t;
	settings *ss;
	double* tf_k;
};


/* Returns a partial value to add to pe_sum */
double niid_sample(struct pthread_args* args){
	t_genome* g = args->g;
	int gid = args->gid;
	t_ptable* pt = args->pt;
	t_afftable* afft = args->afft;
	int t = args->t;
	settings* ss = args->ss;
	double* tf_k = args->tf_k;

	int n_samples = ss->niid / ss->n_pthreads;

	args->return_pes = 0.0;
	for (int ii=0; ii<n_samples; ii++){
		double sum_act = 0.0; // sum of activation energy
		double sum_rep = 0.0; // sum of repression energy
		int s = 0; // site counter
		while (s < g->genes[gid]->urslen){
			if (pt->cpr[s] == 0.0){
				s += 1;
			}
			else{

				int reti, retj, retd;

				/* Sample from the CDF. . . */
				ptable_sample(pt, s, reti, retj, retd);
				// So, gene reti will bind at site, with retj as a cofactor retd distance away.

				// affinity of reti for the sequence starting at site
	#ifdef MEMO_AFFINITY
				double aff = afft->tf_site_aff[reti*g->genes[gid]->urslen + s];
	#else
				double aff = get_affinity(g->genes[reti]->dbd,
						g->genes[gid]->urs,
						g->genes[gid]->urslen,
						s);
	#endif




				if (g->genes[reti]->reg_mode == 0){
					sum_rep += aff;
				}
				else if (g->genes[reti]->reg_mode == 1){
					sum_act += aff;
				}
				/* Advance the site counter */
				s += g->r[reti];
				s += retd;
				if (retj < g->ntfs){
					s += g->r[retj];
				}

				// February 2014:
				tf_k[reti] += aff;
			}

		} // end while s < urslen

		double pe = (1 / (1+exp(-1*ss->pe_scalar*(sum_act-sum_rep) ) ) );

		//printf("\n. fitness 315b - sum_act= %f, sum_rep= %f pe= %f\n", sum_act, sum_rep, pe);

		return_pe += pe;
	}
}
#endif




/* Returns the expression of gene gid, given the probtable ret */
#ifdef MEMO_AFFINITY
double prob_expr(t_genome *g, int gid, t_ptable *pt, t_afftable* afft, int t, settings *ss, double* tf_k){
#else
double prob_expr(t_genome *g, int gid, t_ptable *pt, int t, settings *ss, double* tf_k){
#endif
	double pe_sum = 0.0;

#ifdef PTHREADS
	pthread_attr_t attr;
	pthread_t* threads = (pthread_t *)malloc(ss->n_pthreads*sizeof(pthread_t));
	double* return_pes = (double *)malloc(ss->n_pthreads*sizeof(double));
	int rc;
	int ii;
	void *status;
    /* Initialize and set thread detached attribute */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	/* Arguments for the pthread function */
	struct pthread_args args;
	args.return_pes = return_pes[ii];
	args.g = g;
	args.gid = gid;
	args.pt = pt;
	args.afft = afft;
	args.t = t;
	args.ss = ss;
	args.tf_k = tf_k;
    for(ii=0; ii<ss->n_pthreads; ii++){
		 rc = pthread_create(&threads[ii], NULL, &niid_sample, (void *)&args));
		 if (rc){
			 printf("ERROR; return code from pthread_create() is %d\n", rc);
			 exit(-1);
		 }
	}

	pthread_attr_destroy(&attr);
	for(ii=0; ii<ss->n_pthreads; ii++) {
		rc = pthread_join(threads[tt], &status);
		if (rc) {
			printf("ERROR; return code from pthread_join() is %d\n", rc);
			exit(-1);
		}
		printf("Main: completed join with thread %ld having a status of %ld\n",t,(long)status);"
		pe_sum += return_pes[ii];
	}
	 pthread_exit(NULL);
	 free(threads);
	 free(return_pes);
#else

	int urslen = g->genes[gid]->urslen;
	for (int ii = 0; ii < ss->niid; ii++){ // for each IID sample
		double sum_act = 0.0; // sum of activation energy
		double sum_rep = 0.0; // sum of repression energy
		int s = 0; // site counter
		while (s < urslen){
			if (pt->cpr[s] == 0.0){
				s += 1;
			}
			else{

				int reti, retj, retd;

				/* Sample from the CDF. . . */
				ptable_sample(pt, s, reti, retj, retd);
				// So, gene reti will bind at site, with retj as a cofactor retd distance away.

				// affinity of reti for the sequence starting at site
#ifdef MEMO_AFFINITY
				double aff = afft->tf_site_aff[reti*g->genes[gid]->urslen + s];
#else
				double aff = get_affinity(g->genes[reti]->dbd,
						g->genes[gid]->urs,
						g->genes[gid]->urslen,
						s);
#endif




				if (g->genes[reti]->reg_mode == 0){
					sum_rep += aff;
				}
				else if (g->genes[reti]->reg_mode == 1){
					sum_act += aff;
				}
				/* Advance the site counter */
				s += g->r[reti];
				s += retd;
				if (retj < g->ntfs){
					s += g->r[retj];
				}

				// February 2014:
				tf_k[reti] += aff;
			}

		} // end while s < urslen

		double pe = (1 / (1+exp(-1*ss->pe_scalar*(sum_act-sum_rep) ) ) );

		//printf("\n. fitness 315b - sum_act= %f, sum_rep= %f pe= %f\n", sum_act, sum_rep, pe);

		pe_sum += pe;
#endif
	}

	for (int ii=0; ii<g->ntfs; ii++){
		tf_k[ii] /= ss->niid;
	}

	return pe_sum / (double)ss->niid;
}



