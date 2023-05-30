/******************************************************************************************/
/*lbstools.c - computing properties of line bundle sums */
/******************************************************************************************/


/******************************************************************************************/
/* compute first Chern character of a line bundle sum */

void ch1lbs(int ch1[MAXPS], struct linebundlesum lbs)
{

  int i, a;
  extern struct CYdata CY;

  for (i=0; i<CY.numps; i++) {
    ch1[i]=0;
    for (a=0; a < lbs.rank; a++) ch1[i]=ch1[i]+lbs.entries[i][a];
  }
  
}


/******************************************************************************************/
/* compute second Chern character of a line bundle sum */

void ch2lbs(int ch2[MAXPS], struct linebundlesum lbs)
{

  int i, j, k, a;
  extern struct CYdata CY;

  for (i=0; i<CY.numps; i++) {
    ch2[i]=0;
    for (j=0; j<CY.numps; j++)
      for (k=0; k<CY.numps; k++)
	for (a=0; a<lbs.rank; a++) ch2[i]=ch2[i]+(CY.isec)[i][j][k]*lbs.entries[j][a]*lbs.entries[k][a];
    ch2[i]=ch2[i]/2;
  }

}


/******************************************************************************************/
/* compute indices of line bundles in line bundle sum */

void indlstlbs(int ind[MAXRANK], struct linebundlesum lbs)
{

  int i, j, k, a;
  extern struct CYdata CY;

  for (a=0; a<lbs.rank; a++) {
    ind[a]=0;
    for (i=0; i<CY.numps; i++) {
      ind[a]=ind[a]+(CY.c2TX)[i]*lbs.entries[i][a];
      for (j=0; j<CY.numps; j++)
        for (k=0; k<CY.numps; k++)
           ind[a]=ind[a]+2*(CY.isec)[i][j][k]*lbs.entries[i][a]*lbs.entries[j][a]*lbs.entries[k][a];
    }
    ind[a]=ind[a]/12;
  }
   
}


/******************************************************************************************/
/* compute total index line bundle sum */

int indlbs(struct linebundlesum lbs)
{

  int ind, indlst[MAXRANK];
  int a;
  
  indlstlbs(indlst,lbs);
  ind=0;
  for (a=0; a<lbs.rank; a++) ind=ind+indlst[a];

  return ind;
  
}


/******************************************************************************************/
/* generate random line bundle sum */

struct linebundlesum randomlbs(int rank)
{
  int i, a;
  struct linebundlesum lbs;
  extern struct envsettings set;
  extern struct CYdata CY;
  
  lbs.rank=rank;
  for (i=0; i<CY.numps; i++)
    for (a=0; a<rank; a++)
      lbs.entries[i][a]=randomchoice(set.probdis,set.maxentry-set.minentry+1)+set.minentry;
  
  return lbs;

}


/******************************************************************************************/
/* generate random line bundle sum with c1=0 */

struct linebundlesum randomlbs0(int rank)
{

  int i, inrange, c1zero=0; 
  struct linebundlesum lbs;
  int ch1[MAXPS];
  extern struct envsettings set;
  extern struct CYdata CY;

  while (!c1zero) {
    lbs=randomlbs(rank-1);
    ch1lbs(ch1,lbs); inrange=1; i=0;
    while (inrange && i<CY.numps) {
      if ((-ch1[i]<set.minentry) || (-ch1[i]>set.maxentry)) inrange=0;
      i++;
    }
    c1zero=inrange;
  }
  for (i=0; i<CY.numps; i++) lbs.entries[i][rank-1]=-ch1[i];
  lbs.rank=rank;

  return lbs;
  
}


/******************************************************************************************/
/* test if matrix has positive and negative entries  */

int matposneg(int mat[MAXPS][MAXPS], int size)
{

  int i, a;
  int pos=0, neg=0;

  i=0;
  while ((!pos || !neg) && i<size) {
    for (a=0; a<size; a++) {
      if (mat[i][a]>0) pos=1;
      if (mat[i][a]<0) neg=1;
    }
    i++;
  }

  return (pos && neg) || (!pos && !neg);

}


/******************************************************************************************/
/* test for slope zero of a rank 5 line bundle sum */

int slopezero5(struct linebundlesum lbs)
{

 
  extern struct CYdata CY;
  /*int tcoeff[40][4]={{1,-1,-1,-1},{1,-1,-1,1},{1,-1,1,-1},{1,-1,1,1},{1,1,-1,-1},{1,1,-1,1},
		     {1,1,1,-1},{1,1,1,1},{0,1,-1,-1},{0,1,-1,1},{0,1,1,-1},{0,1,1,1},{1,0,-1,-1},
		     {1,0,-1,1},{1,0,1,-1},{1,0,1,1},{1,-1,0,-1},{1,-1,0,1},{1,1,0,-1},{1,1,0,1},
		     {1,-1,-1,0},{1,-1,1,0},{1,1,-1,0},{1,1,1,0},{1,-1,0,0},{1,1,0,0},{1,0,-1,0},
		     {1,0,1,0},{0,1,-1,0},{0,1,1,0},{0,1,0,-1},{0,1,0,1},{0,0,1,-1},{0,0,1,1},
		     {1,0,0,-1},{1,0,0,1},{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};*/

  int tcoeff[312][4] = {{0,0,0,1},{0,0,0,2},{0,0,1,-2},{0,0,1,-1},{0,0,1,0},{0,0,1,1},{0,0,1,2},{0,0,2,-2},{0,0,2,-1},{0,0,2,0},{0,0,2,1},{0,0,2,2},{0,1,-2,-2},{0,1,-2,-1},{0,1,-2,0},{0,1,-2,1},{0,1,-2,2},{0,1,-1,-2},{0,1,-1,-1},{0,1,-1,0},{0,1,-1,1},{0,1,-1,2},{0,1,0,-2},{0,1,0,-1},{0,1,0,0},{0,1,0,1},{0,1,0,2},{0,1,1,-2},{0,1,1,-1},{0,1,1,0},{0,1,1,1},{0,1,1,2},{0,1,2,-2},{0,1,2,-1},{0,1,2,0},{0,1,2,1},{0,1,2,2},{0,2,-2,-2},{0,2,-2,-1},{0,2,-2,0},{0,2,-2,1},{0,2,-2,2},{0,2,-1,-2},{0,2,-1,-1},{0,2,-1,0},{0,2,-1,1},{0,2,-1,2},{0,2,0,-2},{0,2,0,-1},{0,2,0,0},{0,2,0,1},{0,2,0,2},{0,2,1,-2},{0,2,1,-1},{0,2,1,0},{0,2,1,1},{0,2,1,2},{0,2,2,-2},{0,2,2,-1},{0,2,2,0},{0,2,2,1},{0,2,2,2},{1,-2,-2,-2},{1,-2,-2,-1},{1,-2,-2,0},{1,-2,-2,1},{1,-2,-2,2},{1,-2,-1,-2},{1,-2,-1,-1},{1,-2,-1,0},{1,-2,-1,1},{1,-2,-1,2},{1,-2,0,-2},{1,-2,0,-1},{1,-2,0,0},{1,-2,0,1},{1,-2,0,2},{1,-2,1,-2},{1,-2,1,-1},{1,-2,1,0},{1,-2,1,1},{1,-2,1,2},{1,-2,2,-2},{1,-2,2,-1},{1,-2,2,0},{1,-2,2,1},{1,-2,2,2},{1,-1,-2,-2},{1,-1,-2,-1},{1,-1,-2,0},{1,-1,-2,1},{1,-1,-2,2},{1,-1,-1,-2},{1,-1,-1,-1},{1,-1,-1,0},{1,-1,-1,1},{1,-1,-1,2},{1,-1,0,-2},{1,-1,0,-1},{1,-1,0,0},{1,-1,0,1},{1,-1,0,2},{1,-1,1,-2},{1,-1,1,-1},{1,-1,1,0},{1,-1,1,1},{1,-1,1,2},{1,-1,2,-2},{1,-1,2,-1},{1,-1,2,0},{1,-1,2,1},{1,-1,2,2},{1,0,-2,-2},{1,0,-2,-1},{1,0,-2,0},{1,0,-2,1},{1,0,-2,2},{1,0,-1,-2},{1,0,-1,-1},{1,0,-1,0},{1,0,-1,1},{1,0,-1,2},{1,0,0,-2},{1,0,0,-1},{1,0,0,0},{1,0,0,1},{1,0,0,2},{1,0,1,-2},{1,0,1,-1},{1,0,1,0},{1,0,1,1},{1,0,1,2},{1,0,2,-2},{1,0,2,-1},{1,0,2,0},{1,0,2,1},{1,0,2,2},{1,1,-2,-2},{1,1,-2,-1},{1,1,-2,0},{1,1,-2,1},{1,1,-2,2},{1,1,-1,-2},{1,1,-1,-1},{1,1,-1,0},{1,1,-1,1},{1,1,-1,2},{1,1,0,-2},{1,1,0,-1},{1,1,0,0},{1,1,0,1},{1,1,0,2},{1,1,1,-2},{1,1,1,-1},{1,1,1,0},{1,1,1,1},{1,1,1,2},{1,1,2,-2},{1,1,2,-1},{1,1,2,0},{1,1,2,1},{1,1,2,2},{1,2,-2,-2},{1,2,-2,-1},{1,2,-2,0},{1,2,-2,1},{1,2,-2,2},{1,2,-1,-2},{1,2,-1,-1},{1,2,-1,0},{1,2,-1,1},{1,2,-1,2},{1,2,0,-2},{1,2,0,-1},{1,2,0,0},{1,2,0,1},{1,2,0,2},{1,2,1,-2},{1,2,1,-1},{1,2,1,0},{1,2,1,1},{1,2,1,2},{1,2,2,-2},{1,2,2,-1},{1,2,2,0},{1,2,2,1},{1,2,2,2},{2,-2,-2,-2},{2,-2,-2,-1},{2,-2,-2,0},{2,-2,-2,1},{2,-2,-2,2},{2,-2,-1,-2},{2,-2,-1,-1},{2,-2,-1,0},{2,-2,-1,1},{2,-2,-1,2},{2,-2,0,-2},{2,-2,0,-1},{2,-2,0,0},{2,-2,0,1},{2,-2,0,2},{2,-2,1,-2},{2,-2,1,-1},{2,-2,1,0},{2,-2,1,1},{2,-2,1,2},{2,-2,2,-2},{2,-2,2,-1},{2,-2,2,0},{2,-2,2,1},{2,-2,2,2},{2,-1,-2,-2},{2,-1,-2,-1},{2,-1,-2,0},{2,-1,-2,1},{2,-1,-2,2},{2,-1,-1,-2},{2,-1,-1,-1},{2,-1,-1,0},{2,-1,-1,1},{2,-1,-1,2},{2,-1,0,-2},{2,-1,0,-1},{2,-1,0,0},{2,-1,0,1},{2,-1,0,2},{2,-1,1,-2},{2,-1,1,-1},{2,-1,1,0},{2,-1,1,1},{2,-1,1,2},{2,-1,2,-2},{2,-1,2,-1},{2,-1,2,0},{2,-1,2,1},{2,-1,2,2},{2,0,-2,-2},{2,0,-2,-1},{2,0,-2,0},{2,0,-2,1},{2,0,-2,2},{2,0,-1,-2},{2,0,-1,-1},{2,0,-1,0},{2,0,-1,1},{2,0,-1,2},{2,0,0,-2},{2,0,0,-1},{2,0,0,0},{2,0,0,1},{2,0,0,2},{2,0,1,-2},{2,0,1,-1},{2,0,1,0},{2,0,1,1},{2,0,1,2},{2,0,2,-2},{2,0,2,-1},{2,0,2,0},{2,0,2,1},{2,0,2,2},{2,1,-2,-2},{2,1,-2,-1},{2,1,-2,0},{2,1,-2,1},{2,1,-2,2},{2,1,-1,-2},{2,1,-1,-1},{2,1,-1,0},{2,1,-1,1},{2,1,-1,2},{2,1,0,-2},{2,1,0,-1},{2,1,0,0},{2,1,0,1},{2,1,0,2},{2,1,1,-2},{2,1,1,-1},{2,1,1,0},{2,1,1,1},{2,1,1,2},{2,1,2,-2},{2,1,2,-1},{2,1,2,0},{2,1,2,1},{2,1,2,2},{2,2,-2,-2},{2,2,-2,-1},{2,2,-2,0},{2,2,-2,1},{2,2,-2,2},{2,2,-1,-2},{2,2,-1,-1},{2,2,-1,0},{2,2,-1,1},{2,2,-1,2},{2,2,0,-2},{2,2,0,-1},{2,2,0,0},{2,2,0,1},{2,2,0,2},{2,2,1,-2},{2,2,1,-1},{2,2,1,0},{2,2,1,1},{2,2,1,2},{2,2,2,-2},{2,2,2,-1},{2,2,2,0},{2,2,2,1},{2,2,2,2}};

  int Misec[4][MAXPS][MAXPS], mat[MAXPS][MAXPS];
  int i, j, k, a, b, c, sum, slopemeasure;

  if (lbs.rank!=5) return -1;
  else {
  
    /* compute matrices to be linearly combined */
    for (a=0; a<4; a++)
      for (i=0; i<CY.numps; i++)
        for (j=0; j<CY.numps; j++) {
	  sum=0;
	  for (k=0; k<CY.numps; k++) sum=sum+(CY.isec)[i][j][k]*(lbs.entries)[k][a];
	  Misec[a][i][j]=sum;
        }

    /* check pairs of line bundles */
    /*slopemeasure=0; 
    for (c=0; c<40; c++) {
      for (i=0; i<CY.numps; i++)
        for (j=0; j<CY.numps; j++) {
          mat[i][j]=0;
	  for (a=0; a<4; a++) mat[i][j]=mat[i][j]+tcoeff[c][a]*Misec[a][i][j];
        }*/
    slopemeasure=0; 
    for (c=0; c<312; c++) {
      for (i=0; i<CY.numps; i++)
        for (j=0; j<CY.numps; j++) {
          mat[i][j]=0;
	  for (a=0; a<4; a++) mat[i][j]=mat[i][j]+tcoeff[c][a]*Misec[a][i][j];
        }
      
      /* printf("c = %i, slope = %i\n",c,slopemeasure);
      for (i=0; i<CY.numps; i++) {
	for (j=0; j<CY.numps; j++) printf("%i ",mat[i][j]);
	printf("\n");
      }
      printf("\n"); */
      
      if (!matposneg(mat,CY.numps)) slopemeasure--;
     }
	
    return slopemeasure;
  }
  
} 


 /******************************************************************************************/
/* test if a'th line bundle in a line bundle sum is trivial */

int lbtrivial(struct linebundlesum lbs, int nP, int a)
{

  int i=0, trivial=1;
  
  while (trivial && i<nP) {trivial = trivial && (!lbs.entries[i][a]); i++;}

  return trivial;

}  


 /******************************************************************************************/
/* test if a1'th and fac*a2'th line bundles in a line bundle sum are the same */

int lbssame(struct linebundlesum lbs, int nP, int a1, int a2, int fac)
{

  int i=0, same=1;

  if (a1==a2) return 1;
  else {
    while (same && i<nP) {same = same && (lbs.entries[i][a1]==fac*lbs.entries[i][a2]); i++;}
    return same;
  }

}


/******************************************************************************************/
/* find multiplicities of line bundles in a line bundle sum */

void lbsmul(int mul[MAXRANK], struct linebundlesum lbs, int nP)
{

  int a, b;

  for (a=0; a<lbs.rank; a++) {
    mul[a]=0;
    for (b=0; b<lbs.rank; b++)
      if (lbssame(lbs,nP,a,b,1)) mul[a]++;
  }

}



 /******************************************************************************************/
/* check if line bundle in position a of line bundle sum is greater than line bundle lb2 */

int lbgreater(struct linebundlesum lbs, int a, int lb2[MAXPS], int nP)
{

  int diff=0, i=0;

  if (a<0) return 0;
  else {
    while ((diff==0) && i<nP) {
      diff=lbs.entries[i][a]-lb2[i];
      i++;
    }

    return diff;
  }

}


 /******************************************************************************************/
/* sort line bundles in a line bundle sum */

void sortlbs(struct linebundlesum *lbs)
{

  extern struct CYdata CY;
  int a, ipos, i;
  int activelb[MAXPS];

  for (a=1; a<lbs->rank; a++) {
    for (i=0; i<CY.numps; i++) activelb[i]=(lbs->entries)[i][a];
    ipos=a;
    while (ipos>0 && (lbgreater(*lbs,ipos-1,activelb,CY.numps)>0)) {
      for (i=0; i<CY.numps; i++) (lbs->entries)[i][ipos]=(lbs->entries)[i][ipos-1];
      ipos--;
    }
    for (i=0; i<CY.numps; i++) (lbs->entries)[i][ipos]=activelb[i];
  }
  
}   
    
   
  
 /******************************************************************************************/
/* check if two line bundle sums are equal modulo permutations of line bundles */

int lbsequal(struct linebundlesum lbs1, struct linebundlesum lbs2)
{

  extern struct CYdata CY;
  int a, i=0, diff=0;

  if (lbs1.rank != lbs2.rank) return 0;
  else {
    sortlbs(&lbs1); sortlbs(&lbs2);
    while (i<CY.numps && diff==0) {
      diff=0;
      for (a=0; a<lbs1.rank; a++) diff=diff+abs(lbs1.entries[i][a]-lbs2.entries[i][a]);
      i++;
    }

    return !diff;
    
  }
  
}


 /******************************************************************************************/
/* check if two line bundle sums are equivalent modulo permutations of line bundles and ambient symmetries */

int lbsequiv(struct linebundlesum lbs1, struct linebundlesum lbs2)
{

  extern struct CYdata CY;
  struct linebundlesum lbsperm;
  int i, j, a, k=0, equiv=0;

  if (lbs1.rank != lbs2.rank) return 0;
  else {
    lbsperm.rank=lbs1.rank;
    while (k<CY.numsymm && !equiv) {
      for (i=0; i<CY.numps; i++)
        for (a=0; a<lbs1.rank; a++) {
	  lbsperm.entries[i][a]=0;
	  for (j=0; j<CY.numps; j++) lbsperm.entries[i][a]=lbsperm.entries[i][a]+CY.symm[k][i][j]*lbs1.entries[j][a];
	}
      equiv=lbsequal(lbsperm,lbs2); k++;
    }

    return equiv;
    
  }
      
}


 /******************************************************************************************/
/* perform tensor product between lbs1 and conj*lbs2, where conj=1 or conj=-1 */

struct linebundlesum linetensor(struct linebundlesum lbs1, struct linebundlesum lbs2, int conj)
{

  extern struct CYdata CY;
  struct linebundlesum lbstensor;
  int i, a=0, a1, a2;

  lbstensor.rank=lbs1.rank*lbs2.rank;
  if (lbstensor.rank<=MAXRANK)
    for (a1=0; a1<lbs1.rank; a1++)
      for (a2=0; a2<lbs2.rank; a2++) {
        for (i=0; i<CY.numps; i++) lbstensor.entries[i][a]=lbs1.entries[i][a1]+conj*lbs2.entries[i][a2];
	a++;
      }

  return lbstensor;

}
	

/******************************************************************************************/
/* perform second wedge power of line bundle sum */

struct linebundlesum linewedge2(struct linebundlesum lbs)
{

  extern struct CYdata CY;
  struct linebundlesum lbswedge;
  int i, a=0, a1, a2;

  lbswedge.rank=(lbs.rank*(lbs.rank-1))/2;
   if (lbswedge.rank<=MAXRANK)
     for (a1=0; a1<lbs.rank-1; a1++)
       for (a2=a1+1; a2<lbs.rank; a2++) {
	  for (i=0; i<CY.numps; i++) lbswedge.entries[i][a]=lbs.entries[i][a1]+lbs.entries[i][a2];
	  a++;
       }

   return lbswedge;

}


/******************************************************************************************/
/* compute line bundle cohomology */

int cohline(int coh[4], int lb[MAXPS])
{

  extern struct CYdata CY;

  switch (CY.cicynum) {
  case 7862:
    coh[0]=h07862(lb[0],lb[1],lb[2],lb[3]);
    coh[1]=h17862(lb[0],lb[1],lb[2],lb[3]);
    coh[2]=h17862(-lb[0],-lb[1],-lb[2],-lb[3]);
    coh[3]=h07862(-lb[0],-lb[1],-lb[2],-lb[3]);
    return 1;
  case 7447:
    coh[0]=h07447(lb[0],lb[1],lb[2],lb[3],lb[4]);
    coh[1]=h17447(lb[0],lb[1],lb[2],lb[3],lb[4]);
    coh[2]=h17447(-lb[0],-lb[1],-lb[2],-lb[3],-lb[4]);
    coh[3]=h07447(-lb[0],-lb[1],-lb[2],-lb[3],-lb[4]);
    return 1;
  default: return 0;
  }

}


/******************************************************************************************/
/* compute cohomology of line bundle sum */

int cohlbs(int coh[MAXRANK][4], struct linebundlesum lbs)
{

  extern struct CYdata CY;
  int i, a, ret, co[4], lb[MAXPS];

  for (a=0; a<lbs.rank; a++) {
    for (i=0; i<CY.numps; i++) lb[i]=lbs.entries[i][a];
    ret=cohline(co,lb);
    for (i=0; i<4; i++) coh[a][i]=co[i];
  }

  return ret;

}
  
  
