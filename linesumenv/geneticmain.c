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
 
  /* read in options and CY data from files */
  fp=fopen("CYdata/CICY7862","r");
  confenv(fp);
  fclose(fp);
  
  /* evolve a random initial population */
  evol=evolvepop(randompop(300),300,ROULETTE,STDNUMCUTS,KEEPFITEST,0.003,STDALPHA-1.5,MONITORON);

  /* remove redundancy */
  
  
  /* extract terminal states from evolution */
  bl=termstatesred(evol,300,&nterm);
  printf("no. of different terminal states: %i\n",nterm);
  lbm=removeredlbm(bl,&nterm);
  printf("no. of non-redundant terminal states: %i\n",nterm);
   
  /* write terminal states to a file */
  fp=fopen("data/test.m","a");
  for (i=0; i<nterm; i++) fprintlbm(fp,lbm[i]);
  fclose(fp);	   
}  
