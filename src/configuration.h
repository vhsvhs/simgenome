#define DEF_VERBOSITY 2

#define N_STATES 4 /* for nucleotide regulatory sequences */

#define GENE_NAME_MAX 100
#define FILEPATH_LEN_MAX 300

/* MAXLEN is the max line length for user-specified file.
 * This value should always be longer than URSLEN_DEFAULT
 */
#define MAXLEN 1200
#define MAX_TOKENS 6
#define MAX_PSAM_LEN 100

#define MAX_DDG 2.0
#define MIN_DDG -2.0
#define MAX_GD 3
#define MIN_TF_SEPARATION 0

#define PWMLENMU 0.3
#define PWMLENMUMAX 2
#define DDGMU 0.5

#define URSLEN_DEFAULT 100
#define PSAMLEN_DEFAULT 6
#define NGENES_DEFAULT 10

#define MINIMUM_EXPRESSION_LEVEL 0.00001
#define MAXIMUM_EXPRESSION_LEVEL 1.0

#define NIID 4000 // N I.I.D. samples

#define GROWTH_RATE 1.0 /* Controls how fast genes are turned on. */
#define DECAY_RATE 1.0 /* Controls how fast gene expression lowers when RNA Pol. is removed */

/* PE_SCALAR is inversely proportion to URS length.
Higher values make k_act higher, which makes it easer for genes to be expressed.
Lower values make k_act lower, which means that more bound sites are required
for genes to be activated.*/
#define PE_SCALAR 0.01

/* Larger values make the fitness function sharper around
 * optimal expression, smaller values spread the function's
 * hills further out over poor expression. */
#define FITNESS_SCALAR -3.0

/* Output directory path defaults */
#define POPDIR "LOG_POPS"
#define EXPRDIR "LOG_EXPR"
#define FITNESS_DIR "LOG_FITNESS"
#define COOP_DIR "LOG_COOP"
#define DBD_DIR "LOG_DBD"
#define CONFIG_DIR "LOG_OCCUPANCY"
