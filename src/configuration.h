#define DEF_VERBOSITY 3

#define N_STATES 4 /* for nucleotide regulatory sequences */

/* char* constants */
#define GENE_NAME_MAX 100
#define FILEPATH_LEN_MAX 300
#define TMPFILE_MAX 1000

/* MAXLEN is the max line length for user-specified file.
 * This value should always be longer than URSLEN_DEFAULT
 */
#define MAXLEN 2000
#define MAX_TOKENS 20 // maximum expected number of tokens on an input file line
#define MAX_PSAM_LEN 100

#define MAX_DDG 5.0
#define MIN_DDG -5.0
#define MAX_GD 1
#define MIN_TF_SEPARATION 0

#define URSMU 0.01 // how often are URS substitutions?
#define PSAMMU 0.05 // how often are PSAM substitutions?
#define MU_STDEV 0.01

#define PSAMLENMU 0.05 // how often are PSAM indels?
#define PSAMLENSD 0.01
#define PSAMLENSIZEMU 2 // mean size of PSAM indel
#define PSAMLENSIZESD 1
#define PSAMLENMUMAX 2

#define URSLENMU 0.05 // how often are URS indels?
#define URSLENSD 0.01
#define URSLENSIZEMU 3 // mean size of URS indel
#define URSLENSIZESD 1

#define DDGMU 0.0
#define DDGMUSIZE 1.1
#define DDGMUSD 0.0

/* For random-init of the population: */
#define URSLEN_DEFAULT 1000
#define PSAMLEN_DEFAULT 6
#define NGENES_DEFAULT 10

#define MINIMUM_EXPRESSION_LEVEL 0.00001
#define MAXIMUM_EXPRESSION_LEVEL 1.0

#define POPSIZE 24
#define NIID 4000 // N I.I.D. samples

#define ELITE_PROPORTION 0.3

#define MAX_GENS 1000
#define MAX_TIME 1
#define LOG_SAMPLE_STRIDE 50 // every N generations, the configurations will be logged.  Else, they consume way too much disk space.

#define GROWTH_RATE 1.0 /* Controls how fast genes are turned on. */
#define DECAY_RATE 1.0 /* Controls how fast gene expression lowers when RNA Pol. is removed */

/*This rate gets used in the function coopfunc (in Landscape) to
control the rate of decay of the cooperative or competitive interactions
between TFs.  Big values facilitate cooperative binding across long distances,
whereas small values make the coopfunc dropoff quickly.*/
#define V_RATE_OF_COOP_DECAY 3

/* PE_SCALAR is inversely proportional to URS length.
Higher values make k_act higher, which makes it easer for genes to be expressed.
Lower values make k_act lower, which means that more bound sites are required
for genes to be activated.*/
#define PE_SCALAR 0.005

/* Larger values make the fitness function sharper around
 * optimal expression, smaller values spread the function's
 * hills further out over poor expression. */
#define FITNESS_SCALAR -0.0005

/* Output directory path defaults */
#define POPDIR "LOG_POPS"
#define EXPRDIR "LOG_EXPR"
#define FITNESS_DIR "LOG_FITNESS"
#define COOP_DIR "LOG_COOP"
#define DBD_DIR "LOG_DBD"
#define CONFIG_DIR "LOG_OCCUPANCY"


#define TRANPT_SIZE 10
