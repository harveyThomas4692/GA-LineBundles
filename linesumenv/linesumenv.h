/******************************************************************************************/
/* linesumenv.h - environment for line bundle models */
/******************************************************************************************/
/* libraries */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/******************************************************************************************/
/* definitions for basic settings */
/******************************************************************************************/

#define MAXPS 7       /* maximal number of projective factors in ambient space */
#define MAXPOLS 7     /* maximal number of defining polynomials */
#define MAXRANK 25    /* maximal rank of line bundle sum */
#define MAXRANGE 30   /* maximal range of line bundle integers */
#define MAXCONE 20    /* maximal number of cone generatots */
#define MAXSYMM 120   /* maximal number of symmetry generators for permuting ambient space factors */
#define MAXWEIGHTS 10 /* maximal number of contributions to fitness function */
#define MAXORD 10     /* maximal number of possible symmetry orders */

/******************************************************************************************/
/* structures */
/******************************************************************************************/

struct CYdata                            /* structure to store CY data */
{
  int cicynum;                           /* number of CICY */
  int numps;                             /* number of projective ambient factors */
  int numpols;                           /* number of defining polynomials */
  int conf[MAXPS][MAXPOLS];              /* configuration matrix */
  int c2TX[MAXPS];                       /* second Chern class of tangent bundle */
  int numord;                            /* number of possible symmetry orders */
  int symmord[MAXORD];                   /* possible symmetry orders */
  int isec[MAXPS][MAXPS][MAXPS];         /* triple intersection numbers */
  int numkahlerconegen;                  /* number of Kahler cone generators */
  int kahlerconegen[MAXCONE][MAXPS];     /* Kahler cone generators */
  int numskahlerconegen;                 /* number of s-Kahler cone generators */
  int skahlerconegen[MAXCONE][MAXPS];    /* s-Kahler cone generators */
  int numsymm;                           /* number of symmmetry generators */
  int symm[MAXSYMM][MAXPS][MAXPS];       /* symmetry */
  int cohformula;                        /* whether or not a line bundle cohomology formula is available for this CY */
  int nonDiagSym;                        /* the number of non diagonal symmetry generators */
  int nonDiagSymGens[MAXSYMM][MAXPS][MAXPS]; /* the symmetry generators needed for checking equivariance, when symmetry isn't diagonal*/
};


struct envsettings                       /* structure to store environment settings */
{
  int slopezerometh;                     /* method to determine wheter slope 0 is satisfied */
  int rank;                              /* rank of line bundle sum */
  int symmorder;                         /* order of discrete symmetry */
  int wilsonline[2];                     /* Wilson line */
  int downstairs;                        /* whether or not to compute downstairs spectrum */
  int minentry;                          /* minimal entry of line bundle sum */
  int maxentry;                          /* maximal entry of line bundle sum */
  int numbits;                           /* number of bits used per entry */
  int numweights;                        /* number of weights for fitness function contributions */
  int weights[MAXWEIGHTS];               /* weights for fitness function contributions */
  float probdis[MAXRANGE];               /* probability distribution for line bundle integers */
  float termcond;                        /* state considered terminal if fitness >= this value */
};


struct linebundlesum                     /* structure for a line bundle sum */
{
  int rank;                              /* rank of line bundle sum */
  int entries[MAXPS][MAXRANK];           /* entries of line bundle sum */
};


struct lbmodel                           /* structure for a line bundle model */
{
  struct linebundlesum lbs;             /* line bundle sum */
  struct bitlist bl;                    /* bitlist for model */
  int ch2[MAXPS];                       /* ch2 of line bundle sum */
  int indlst[MAXRANK];                  /* indices of individual line bundles */
  int ind;                              /* total index of line bundle sum */
  int slope0;                           /* measure for slope: 0 if there is a solution, negative if not */
  int equiv;                            /* measure for equivariance of line bundle sum */        
  int nOX;                              /* number of trivial line bundles in line bundle sum */
  int nsplits;                          /* number of splits in line bundle sum */
  int cohL[MAXRANK][4];                 /* cohomology of line bundles in line bundle sum */
  int cohV[4];                          /* cohomology of line bundle sum */
  int cohLL[MAXRANK][4];                /* cohomology of line bundles in seconed wedge power of line bundle sum */
  int cohV2[4];                         /* cohomology of second wedge power of line bundle sum */
  int nHiggs;                           /* number of Higgs pairs */
  float valuelst[MAXWEIGHTS];           /* list of contributions to fitness */
};

/******************************************************************************************/
/* external variables */
/******************************************************************************************/

struct envsettings set;
struct CYdata CY;

/******************************************************************************************/
/* fixed data to compute line bundle cohomology */
/******************************************************************************************/

#include "cohdata.h"

/******************************************************************************************/
/* prototypes */
/******************************************************************************************/
/* auxiliary functions */

/* load in global settings and CY data from file */
void confenv(FILE *fp);

/* print global settings */
void printsettings();

/* print CY data */
void printCYdata();

/* print line bundle sum */
void printlbs(struct linebundlesum lbs);

/* write line bundle sum to file */
void fprintlbs(FILE *fp, struct linebundlesum lbs);

/* print a line bundle model */
void printlbm(struct lbmodel lbm);

/* write a line bundle model to a file */
void fprintlbm(FILE *fp, struct lbmodel lbm);

/******************************************************************************************/
/* functions for line bundle sums */

/* compute first Chern character of a line bundle sum */
void ch1lbs(int c1[MAXPS],struct linebundlesum lbs);

/* compute second Chern character of a line bundle sum */
void ch2lbs(int ch2[MAXPS], struct linebundlesum lbs);

/* compute indices of line bundles in line bundle sum */
void indlstlbs(int ind[MAXRANK], struct linebundlesum lbs);

/* compute total index line bundle sum */
int indlbs(struct linebundlesum lbs);

/* generate random line bundle sum */
struct linebundlesum randomlbs(int rank);

/* generate random line bundle sum with c1=0 */
struct linebundlesum randomlbs0(int rank);

/* test if matrix has positive and negative entries  */
int matposneg(int mat[MAXPS][MAXPS], int size);

/* test for slope zero of a rank 5 line bundle sum */
int slopezero5(struct linebundlesum lbs);

/* test if a'th line bundle in a line bundle sum is trivial */
int lbtrivial(struct linebundlesum lbs, int nP, int a);

/* test if a1'th and a2'th line bundles in a line bundle sum are the same */
int lbssame(struct linebundlesum lbs, int nP, int a1, int a2, int fac);

/* find multiplicities of line bundles in a line bundle sum */
void lbsmul(int mul[MAXRANK], struct linebundlesum lbs, int nP);

/* check if line bundle in position a of line bundle sum is greater than line bundle lb2 */
int lbgreater(struct linebundlesum lbs, int a, int lb2[MAXPS], int nP);

/*sort line bundles in a line bundle sum */
void sortlbs(struct linebundlesum *lbs);

/* check if two line bundle sums are equal modulo permutations of line bundles */
int lbsequal(struct linebundlesum lbs1, struct linebundlesum lbs2);

/* check if two line bundle sums are equivalent modulo permutations of line bundles and ambient symmetries */
int lbsequiv(struct linebundlesum lbs1, struct linebundlesum lbs2);

/* perform tensor product between lbs1 and conj*lbs2, where conj=1 or conj=-1 */
struct linebundlesum linetensor(struct linebundlesum lbs1, struct linebundlesum lbs2, int conj);

/* perform second wedge power of line bundle sum */
struct linebundlesum linewedge2(struct linebundlesum lbs);

/* compute line bundle cohomology */
int cohline(int coh[4], int lb[MAXPS]);

/* compute cohomology of line bundle sum */
int cohlbs(int coh[MAXRANK][4], struct linebundlesum lbs);

/******************************************************************************************/
/* line bundle cohomology formulae */

/* compute h0(L) for CICY 7862 */ 
int h07862(int k1, int k2, int k3, int k4);

/* compute h1(L) for CICY 7862 */
int h17862(int k1, int k2, int k3, int k4);

/* compute h0(L) for CICY 7447 */ 
int h07447(int k1, int k2, int k3, int k4, int k5);

/* compute h1(L) for CICY 7447 */ 
int h17447(int k1, int k2, int k3, int k4, int k5);

/******************************************************************************************/
/* main functions */

/* complete a line bundle sum to a model */
void completemodel(struct lbmodel *lbs, int meth);

/* complete a bitlist to a model */
void completebitlist(struct lbmodel *lbm, struct bitlist *bl);

/* check if two line bundle models are equivalent */
int lbmequiv(struct lbmodel *lbm1, struct lbmodel *lbm2);

/* compute line bundle models from list of bitlists and remove redundancies */
struct lbmodel * removeredlbm(struct bitlist *bl, int *len);

