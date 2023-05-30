/******************************************************************************************/
/* testing application of genetic code to line sum environment */
/******************************************************************************************/

/* change paths according to your set-up */
#include "/DIR/genetic.h"
#include "/DIR/linesumenv.h"
#include "/DIR/genetic.c"
#include "/DIR/linesumenv.c"

int main()
{

  struct linebundlesum lbs;
  struct lbmodel *lbm;
  struct bitlist *bl;
  struct population *evol;
  FILE *fp;
  int i, nterm;
  int maxRuns = 150; 
  int satur[maxRuns];

 
  /* read in options and CY data from files */
  fp=fopen("CYdata/CICY7862","r");
  confenv(fp);
  fclose(fp);

  /* struct bitlist * searchenv(int numevol, int numgen, int popsize, int *nterm, int *saturate, int meth, int numcuts,
			   int keepfitest, float mutrate, float alpha, int monitor) */

  bl = searchenv(maxRuns, 200/*500*/, 300, &nterm, satur, ROULETTE, STDNUMCUTS, KEEPFITEST,0.003,STDALPHA-1.5,MONITORON);
  lbm = removeredlbm(bl,&nterm);
  
  /* evolve a random initial population */
  /* evolvepop(randompop(300),300,ROULETTE,STDNUMCUTS,KEEPFITEST,0.003,STDALPHA-1.5,MONITORON); */
  /*evol=evolvepop(randompop(300),300,ROULETTE,STDNUMCUTS,1,0.003,STDALPHA-1.5,1); */
  
  /* extract terminal states from evolution */
  /* bl=termstatesred(evol,300,&nterm); */
  /* printf("no. of different terminal states: %i\n",nterm); */
  /* lbm=removeredlbm(bl,&nterm); */
  /* printf("no. of non-redundant terminal states: %i\n",nterm); */
   
  /* write terminal states to a file */
   fp=fopen("data/4071TermsHydra-Sym2-ResRange.m","a");
   fprintf(fp, "BREAK");
   fprintf(fp, "\n");
   for (i=0; i<nterm; i++) {
      fprintlbm(fp,lbm[i]);
      fprintf(fp, "\n");
   } 
   fclose(fp);

   /* fp=fopen("data/7862Saturtest.m","a"); */
   /* fprintf(fp, "BREAK"); */
   /* fprintf(fp, "\n"); */
   /* for (i=0; i<maxRuns; i++) { */
   /*     fprintf(fp, "%d",satur[i]); */
   /*     fprintf(fp,"\n"); */
   /* } */
   /* fclose(fp);  */


}  
