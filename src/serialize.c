#include "common.h"

/* Writes the population and all its genomes and genes to a file. */
void serialize_population(t_pop* pop, settings* ss){
	char* gc;
	gc = (char*)malloc(10*sizeof(char));
	sprintf(gc, "%d", ss->gen_counter);

	char* p = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat(
			strcat(
			strcat(
			strcat(p, ss->outdir),
			"/POPS/pop.gen"),
			gc),
			".save.txt"
			);

	FILE *fo; /* File for psam specs */
	fo = fopen(p,"w");
	if (fo == NULL) {
	  fprintf(stderr, "Error: can't open output file %s!\n",
			  p);
	  exit(1);
	}

	if (ss->verbosity > 2){
		printf("\n. The population was saved to %s\n", p);
	}

	fprintf(fo, "Written with version %s\n", __VERSION);
	fprintf(fo, "===============================================\n");

	/* For each genome, write a summary about the number and size of genes.
	 * This information is useful for de-serializing the population because
	 * we'll need to malloc a bunch of memory, requiring us to know the size
	 * of the population and the number of genes, etc.
	 */
	fprintf(fo, "N Genomes: %d\n", pop->ngenomes);
	for (int ii = 0; ii < pop->ngenomes; ii++){
		fprintf(fo, "Genome: %d -- %d total genes, %d regulator genes\n",
				ii,
				pop->genomes[ii]->ngenes,
				pop->genomes[ii]->ntfs);

		for (int jj = 0; jj < pop->genomes[ii]->ngenes; jj++){
			if (pop->genomes[ii]->genes[jj]->has_dbd){ /* REGULATOR GENE */
				int psamlen = 0;
				if (pop->genomes[ii]->genes[jj]->has_dbd){
					psamlen = pop->genomes[ii]->genes[jj]->dbd->nsites;
				}
				fprintf(fo, "\tGene: %d role: %s urslen: %d reg_mode: %d psamlen: %d\n",
						jj,
						(pop->genomes[ii]->genes[jj]->has_dbd)?"regulator":"reporter",
						pop->genomes[ii]->genes[jj]->urslen,
						pop->genomes[ii]->genes[jj]->reg_mode,
						psamlen);
			}
			else { /* REPORTER GENE */
				int psamlen = 0;
				if (pop->genomes[ii]->genes[jj]->has_dbd){
					psamlen = pop->genomes[ii]->genes[jj]->dbd->nsites;
				}
				fprintf(fo, "\tGene: %d role: %s urslen: %d \n",
						jj,
						(pop->genomes[ii]->genes[jj]->has_dbd)?"regulator":"reporter",
						pop->genomes[ii]->genes[jj]->urslen);
			}
		}
	}

	/* Write a detailed description of the URSs, DBDs, and co-factor interactions
	 * within each genome.
	 */
	for (int ii = 0; ii < pop->ngenomes; ii++){
		log_urs(pop->genomes[ii], ss, fo);
		log_dbds(pop->genomes[ii], ss, fo);
		log_cofactor(pop->genomes[ii], ss, fo);
	}

	fclose(fo);
	free(gc);
	free(p);
}


t_pop* deserialize_population(settings* ss){

	char* p = (char *)malloc(FILEPATH_LEN_MAX*sizeof(char));
	strcat(p, ss->poppath);

	FILE *fi;
	fi = fopen(p,"r");
	if (fi == NULL) {
		printf("\n. I can't open the file %s.\n", p);
		fprintf(stderr, "Error: can't open output file %s!\n", p);
		exit(1);
	}

	printf("\n. Reading the saved population in %s\n", p);

	int popsize = 0;
	t_pop* pop; // the new population

	char match1 [] = "N Genomes: ";
	char match2 [] = "Genome: ";
	char match3 [] = "Gene: ";
	char match4 [] = "URS ";
	char match5 [] = "PSAM Genome ";
	char match6 [] = "COOP ";
	char match7 [] = "Gene ";

	int this_genome = 0; // temporary counters
	int this_gene = 0;
	int this_site = 0;
	int currmode = -1; // 0 = URS, 1 = PSAM, 2 = COOP

	char line[MAXLEN];
	while (  fgets(line, MAXLEN, fi)  ){

		/* Skip empty lines */
		if (strlen(line) < 2){
			continue;
		}

		/* N Genomes: */
		int x = match(line, match1);
		if (x > -1){
			//printf("seri 103 %s\n", line);
			const char* tokens[MAX_TOKENS] = {};
			tokens[0] = strtok(line, " ");
			tokens[1] = strtok(0, " ");
			tokens[2] = strtok(0, " ");
			popsize = atoi(tokens[2]);
			if (ss->verbosity > 2){
				printf("\n. There are %d unique genomes.", popsize);
			}
			pop = make_population_basic( popsize );
			continue;
		}

		/* Genome
		 *
		 */
		x = match(line, match2);
		if (x > -1){
			//printf("seri 117 %s\n", line);
			const char* tokens[MAX_TOKENS] = {};
			tokens[0] = strtok(line, " ");
			tokens[1] = strtok(0, " ");
			int this_id = atoi(tokens[1]);
			tokens[2] = strtok(0, " ");
			tokens[3] = strtok(0, " ");
			int this_ngenes = atoi(tokens[3]);
			tokens[4] = strtok(0, " ");
			tokens[5] = strtok(0, " ");
			tokens[6] = strtok(0, " ");
			int this_ntfs = atoi(tokens[6]);
			tokens[7] = strtok(0, " ");
			tokens[8] = strtok(0, " ");
			if (ss->verbosity > 3){
				printf("\n. Building genome %d", this_id);
			}
			pop->genomes[this_id] = make_genome(this_ngenes, NULL, ss);
			pop->genomes[this_id]->id = this_id;
			pop->genomes[this_id]->ntfs = this_ntfs;
			this_genome = this_id;
			continue;
		}

		/* Found a Gene:
		 * id = 1
		 * reg_mode = 7
		 * urslen = 5
		 * psamlen = 9
		 * */
		x = match(line, match3);
		if (x > -1){
			//printf("seri 149 %s\n", line);
			const char* tokens[MAX_TOKENS] = {};
			tokens[0] = strtok(line, " ");
			tokens[1] = strtok(0, " ");
			int gid = atoi(tokens[1]);
			tokens[2] = strtok(0, " ");
			tokens[3] = strtok(0, " ");
			tokens[4] = strtok(0, " ");
			tokens[5] = strtok(0, " ");
			int this_urslen = atoi(tokens[5]);
			int this_reg_mode;
			int this_psamlen = 0;
			if (gid < pop->genomes[this_genome]->ntfs){
				tokens[6] = strtok(0, " ");
				tokens[7] = strtok(0, " ");
				this_reg_mode = atoi(tokens[7]);
				tokens[8] = strtok(0, " ");
				tokens[9] = strtok(0, " ");
				this_psamlen = atoi(tokens[9]);
			}
			//printf("\n. Building gene %d genome %d", gid, this_genome);
			pop->genomes[this_genome]->genes[gid] = make_gene(this_psamlen, this_urslen);
			pop->genomes[this_genome]->genes[gid]->id = gid;

			if (this_psamlen > 0){
				build_coop( pop->genomes[this_genome]->genes[gid],
								pop->genomes[this_genome]->ntfs,
								ss->maxgd);
				//init_coop(pop->genomes[this_genome]->genes[gid]);
				pop->genomes[this_genome]->genes[gid]->reg_mode = this_reg_mode;
			}

			this_gene = gid;
			continue;
		}

		x = match(line, match4); //URS
		if (x > -1){
			const char* tokens[MAX_TOKENS] = {};
			tokens[0] = strtok(line, " ");
			tokens[1] = strtok(0, " ");
			tokens[2] = strtok(0, " ");
			this_genome = atoi( tokens[2] );
			this_gene = 0;
			currmode = 0;
			continue;
		}

		x = match(line, match5); //PSAM
		if (x > -1){
			const char* tokens[MAX_TOKENS] = {};
			tokens[0] = strtok(line, " ");
			tokens[1] = strtok(0, " ");
			tokens[2] = strtok(0, " ");
			this_genome = atoi( tokens[2] );
			this_gene = 0;
			currmode = 1;
			continue;
		}

		x = match(line, match6); //COOP
		if (x > -1){
			const char* tokens[MAX_TOKENS] = {};
			tokens[0] = strtok(line, " ");
			tokens[1] = strtok(0, " ");
			tokens[2] = strtok(0, " ");
			this_genome = atoi( tokens[2] );
			currmode = 2;
			continue;
		}

		/* Process URS */
		if (line[0] != '\n' && currmode == 0){
			//printf("\n. seri 233 this_genome %d this_gene %d %s\n", this_genome, this_gene, line);
			const char* tokens[MAX_TOKENS] = {};
			// The FASTA taxa header:
			if (line[0] == '>'){
				tokens[0] = strtok(line, " ");
				tokens[0]++;
				this_gene = atoi( tokens[0] );
			}
			else {
				// For each site in the URS:
				int this_int;
				int count_len = 0;
				for(int ii=0; ii < MAXLEN; ii++){
					this_int = nt2int( line[ii] );
					if (this_int < 0 || this_int > 3){
						ii = MAXLEN;
					}
					else
					{
						pop->genomes[this_genome]->genes[this_gene]->urs[ii] = this_int;

						//printf("\n. seri 252 %c %d\n", line[ii], this_int);
						count_len += 1;
					}
					pop->genomes[this_genome]->genes[this_gene]->urslen = count_len;
				}

			}

		}
		/* Process PSAM */
		else if (line[0] != ' ' && currmode == 1){
			//printf("\n. seri 240, this_site %d PSAM %s\n", this_site, line);
			const char* tokens[MAX_TOKENS] = {};
			x = match(line, match7); // Gene X X
			if (x > -1){
				tokens[0] = strtok(line, " ");
				tokens[1] = strtok(0, " ");
				this_site = 0;
				int this_gene_id = atoi( tokens[1] );
				this_gene = this_gene_id;
				tokens[2] = strtok(0, " ");
				int this_regmod = atoi( tokens[2] );
				pop->genomes[this_genome]->genes[this_gene]->reg_mode = this_regmod;
				pop->genomes[this_genome]->genes[this_gene]->has_dbd = 1;
				continue;
			}
			else{
				//printf("\n.seri 279 line: %s\n", line);
				for (int ii=0; ii < N_STATES; ii++){ // ii = state
					if (ii == 0){ tokens[ii] = strtok(line, " "); }
					else{tokens[ii] = strtok(0, " "); }
					double value = atof( tokens[ii] ); //subsequent tokens
					int index = this_site*N_STATES + ii;
					pop->genomes[this_genome]->genes[this_gene]->dbd->data[index] = value;

				}
				//printf("\n. seri 289 this_site ++ %d\n", this_site );
				this_site += 1;
			}
		}

		/* Process COOP */
		else if (line[0] != ' ' && currmode == 2){
			// process COOP
			int from = atoi( strtok(line, " ") );
			int to = atoi( strtok(0, " "));
			fflush(stdout);
			double value = atof( strtok(0, " " ));
			pop->genomes[this_genome]->genes[from]->tfcoop[to] = value;
		}
	}
	fclose(fi);

	// Jan 2014:
	for (int gid=0; gid < pop->ngenomes; gid++){
		printf("\n. 328 gid: %d", gid);
		for (int geneid=0; geneid < pop->genomes[gid]->ntfs; geneid++){
			printf("\n. 330 geneid: %d", geneid);
			calc_gamma(pop->genomes[gid]->genes[geneid], pop->genomes[gid]->ntfs, ss->maxgd);
		}

	}

	return pop;
}
