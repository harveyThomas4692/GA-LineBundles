/******************************************************************************************/
/* lbsmain.c - main functions in line sum environment */
/******************************************************************************************/

/******************************************************************************************/
/* complete a line bundle sum to a model, meth=0: line bundle sum is given and bitlist computed, meth=1: the opposite */

void completemodel(struct lbmodel *lbm, int meth)
{

  int i, j, k, l, p, a, b, c, dc2, rk ,indmirror, indmirror2, cv;
  int vec[IMAXDIM], mul[MAXRANK], indlst2[MAXRANK];
  int ch1[MAXPS];
  struct linebundlesum lbs2;
  extern struct envsettings set;
  extern struct CYdata CY;
  FILE *fp;
  
  /* convert between bitlist and line bundle sum */
  
  /* if meth=1 convert bitlist into line bundle sum */
  if (meth) {
    bits2intvec(vec,lbm->bl,set.numbits,set.minentry);
    for (a=0; a<set.rank-1; a++)
      for (i=0; i<CY.numps; i++) (lbm->lbs).entries[i][a]=vec[a*CY.numps+i];
    rk=(lbm->bl).len/(set.numbits)/(CY.numps); (lbm->lbs).rank=rk;
    ch1lbs(ch1,lbm->lbs); 
    for (i=0; i<CY.numps; i++) (lbm->lbs).entries[i][rk]=-ch1[i]; 
    (lbm->lbs).rank=rk+1;
  }
  /* if meth=0 convert line bundle sum into bitlist */
  else {
    for (a=0; a<set.rank-1; a++)
       for (i=0; i<CY.numps; i++) vec[a*CY.numps+i]=(lbm->lbs).entries[i][a];
    lbm->bl=intvec2bits(vec,(CY.numps)*((lbm->lbs).rank-1),set.numbits,set.minentry);
  } 

  /* compute ch2(v) */
  ch2lbs(lbm->ch2,lbm->lbs);

  /* compute indices of line bundles and total index */
  indlstlbs(lbm->indlst,lbm->lbs);
  lbm->ind=0; indmirror=0;
  for (a=0; a<(lbm->lbs).rank; a++) {
    lbm->ind=(lbm->ind)+(lbm->indlst)[a]; /* total index */
    if ((lbm->indlst)[a]>0) indmirror=indmirror+(lbm->indlst)[a]; /* sum of positive indices -> chiral mirror families */
  }
  
  /* compute indices of line bundles in second wedge power and check for positive values */
  lbs2=linewedge2(lbm->lbs);
  indlstlbs(indlst2,lbs2);
  indmirror2=0;
  for (a=0; a<lbs2.rank; a++) if (indlst2[a]>0) indmirror2=indmirror2+indlst2[a];
  
  /* compute measure for common vanishing slope */
  lbm->slope0=slopezero5(lbm->lbs);




  /* check equivariance */
  /* Check if the symmetry is diagonal */
  if(CY.nonDiagSym == 0){/* If diagonal use the usual equivariance check*/
    lbsmul(mul,lbm->lbs,CY.numps); /* compute multiplicities of line bundle in line bundle sum */
    lbm->equiv=0;
    for (a=0; a<(lbm->lbs).rank; a++) {
       lbm->equiv=lbm->equiv-abs(((lbm->indlst)[a]*mul[a]) % (set.symmorder));
    }
  }
  else if (CY.nonDiagSym == 1){/* If non-diagonal apply P to the list of line bundles */
    lbm->equiv=-50;
    struct linebundlesum lbsTemp = randomlbs((lbm->lbs).rank);

    for(l=0;l < (lbm->lbs).rank; l++){
      for(j=0;j<CY.numps;j++){
        lbsTemp.entries[j][l] = 0;
        for(i=0;i<CY.numps;i++){
          lbsTemp.entries[j][l] += (CY.nonDiagSymGens)[0][j][i] * (lbm->lbs).entries[i][l];
        }
      }
    }
    
    /* Now need to compare the new line bundle witn the old after every possible permutation from symmetry, then check index of each of those */
    fp=fopen("PermData/PERMS5","r");/*Currently just hardcoding the file of interest */
    int permsNumber,rk, par;
    struct linebundlesum lbsTemp2 = randomlbs((lbm->lbs).rank);
    fscanf(fp,"Rank: %i\nPerms: %i\n",&rk,&permsNumber);
    int perm[rk][rk];

    if(rk!= (lbm->lbs).rank){/* Return an error if the file is for the incorrect rank*/
      printf("Error - The rank in the permutation file does not match the rank of the line bundle.");
      exit(0);
    }

    for(p=0; p<permsNumber; p++){/*Loop over all perumutations of line bundles*/
      /*Read in the matrix*/
      fscanf(fp,"Matrix: \n");
      for (j=0; j<rk; j++) {
        for (k=0; k<rk; k++) fscanf(fp,"%i",&(perm[j][k])); 
        fscanf(fp,"\n");
      }
     /* Now we have the matrix, we need to apply it to the temp lineBundleSum and see if it is the same as the original */
      for(l=0;l < CY.numps; l++){
        for(j=0;j<rk;j++){
          lbsTemp2.entries[l][j] = 0;
          for(i=0;i<rk;i++){
            lbsTemp2.entries[l][j] += perm[j][i] * lbsTemp.entries[l][i];
          }
        }
      }
      if(lbsequal(lbm->lbs,lbsTemp2)==1){/*If the two line bundles are the same, we need to look at equivarience checks*/
         fscanf(fp,"Parts: %i\n",&par);/* Find out how many partitions there are */
         int parts[par];
         fscanf(fp,"PartLen: ");
         for(i=0; i<par; i++) fscanf(fp, "%i", &(parts[i])); /*Find out the length of each partition*/
         c=0;
         for(j=0;j<par;j++){/* Compute index for each partition */
           b=0;
           for(i=0;i<parts[j];i++){/* For each partition, find the contribution from each line bundle */
              fscanf(fp,"%i", &a);
              b=b+((lbm->indlst)[a]);
           }
           c = c - (abs(b % (set.symmorder))/par);
           fscanf(fp, "\n");
         }
         fscanf(fp, "\n");
         if(c>lbm->equiv){ /*Check if this is the best check so far, and change fitness as necissary */
           lbm->equiv = c;
           if(c==0)p=permsNumber;
         }
      }
      else{
         b=0;
         while (i<CY.numps) {
           for (c=0; c<rk; c++) b=b-abs(lbm->lbs.entries[i][c]-lbsTemp2.entries[i][c]);
           i++;
         }
         if(b-(set.symmorder)>lbm->equiv){ /*Check if this is the best check so far, and change fitness as necissary */
           lbm->equiv = b-(set.symmorder);
         }
      }
    }
  fclose(fp);
  }

  else { /*With how fidly this is + the small chance of succsess (even worse with more symmetries) currently running with just order two symmetries in mind*/
    lbm->equiv=0;
    printf("Currently only functionality for order 2 non-diagonal symmetries is implemented \n"); 
  }





  /* count number of trivial line bundles in line bundle sum */
  lbm->nOX=0;
   for (a=0; a<(lbm->lbs).rank; a++)
     if (lbtrivial(lbm->lbs,CY.numps,a)) (lbm->nOX)++;

  /* find splits in line bundle sum - enough to check if two line bundles sum to zero */
  lbm->nsplits=0;
  for (a=0; a<(lbm->lbs).rank-1; a++)
    for (b=a+1; b<(lbm->lbs).rank; b++)
      lbm->nsplits=(lbm->nsplits)+lbssame(lbm->lbs,CY.numps,a,b,-1);

  /* if formulae is implemented compute cohomology of line bundle sums */
  if (CY.cohformula) {
    cohlbs(lbm->cohL,lbm->lbs);
    for (i=0; i<4; i++) {
      (lbm->cohV)[i]=0;
      for (a=0; a<(lbm->lbs).rank; a++) (lbm->cohV)[i]=(lbm->cohV)[i]+(lbm->cohL)[a][i];
    }

    /* compute cohomology of second wedge power of line bundle sum */
    cohlbs(lbm->cohLL,lbs2);
    for (i=0; i<4; i++) {
      (lbm->cohV2)[i]=0;
      for (a=0; a<lbs2.rank; a++) (lbm->cohV2)[i]=(lbm->cohV2)[i]+(lbm->cohLL)[a][i];
    }

    /* compute number of Higgs pairs */
    lbm->nHiggs=(lbm->cohV2)[2];
  }
  /* otherwise, set cohomology entries to standard values */
  else {
    for (i=0; i<4; i++) {
      if (i) cv=0; else cv=-1;
      (lbm->cohV)[i]=cv;
      (lbm->cohV2)[i]=cv;
      rk=(lbm->lbs).rank;
      for (a=0; a<rk; a++) (lbm->cohL)[a][i]=cv;
      for (a=0; a<(rk*(rk-1))/2; a++) (lbm->cohLL)[a][i]=cv;
    }
  }

  /* compute contributions to fitness */
  for (i=0; i<set.numweights; i++) (lbm->valuelst)[i]=0;
  /* contribution form anomaly */
  for (i=0; i<CY.numps; i++) {
    dc2=(CY.c2TX)[i]+(lbm->ch2)[i];
    if (dc2<0) (lbm->valuelst)[0]=(lbm->valuelst)[0]+10.*dc2/(CY.numps)/pow(set.maxentry,2)/(lbm->lbs).rank;
  }
  /* contribution from index */
  (lbm->valuelst)[1]=(-100.*(abs((lbm->ind)+3*(set.symmorder))+indmirror+indmirror2))/(CY.numps)/pow(set.maxentry,3)/(lbm->lbs).rank;
  /* contribution from slope */
  (lbm->valuelst)[2]=(1.*(lbm->slope0))/10.;
  /* contribution from equivariance */
  (lbm->valuelst)[3]=1.*(lbm->equiv);
  /* contribution from reduced structure group */
  (lbm->valuelst)[4]=(-1.*((lbm->nOX)+(lbm->nsplits)))/10.;

  /* contributions related to cohomology */
  if (CY.cohformula) {
    /* contribution form non-zero h0 and h3 of line bundle sum */
    (lbm->valuelst)[5]=(-1000.*((lbm->cohV)[0]+(lbm->cohV)[3]+(lbm->cohV2)[0]+(lbm->cohV2)[3]))/
      CY.numps/pow(set.maxentry,3)/pow((lbm->lbs).rank,2);
    /* penalty for mirror generations */
    (lbm->valuelst)[6]=(-100.*lbm->cohV[2])/CY.numps/set.maxentry;
    /* penalty for absence of Higgs */
    if ((lbm->nHiggs)==0) (lbm->valuelst)[7]=-0.1;
    /* penalty for too many Higgs pairs */
    if ((lbm->nHiggs)>2*(set.symmorder))
      (lbm->valuelst)[8]=-5.*((1.*(lbm->nHiggs))-2.)/CY.numps/pow(set.maxentry,3)/(lbm->lbs).rank;
  }
    
  /* compute fitness */
  (lbm->bl).fitness=0.;
  for (i=0; i<set.numweights; i++) (lbm->bl).fitness=(lbm->bl).fitness+(lbm->valuelst)[i]*(set.weights)[i];

  /* check if terminal */
  (lbm->bl).terminal=((lbm->bl).fitness)>=set.termcond;
   
}


/******************************************************************************************/
/* complete a bitlist to a model */

void completebitlist(struct lbmodel *lbm, struct bitlist *bl)
{

  lbm->bl=*bl;
  completemodel(lbm,1);

}

/******************************************************************************************/
/* compute fitness of bitlist for line bundle model */

void fitness(struct bitlist *bl)
{

  struct lbmodel lbm;

  completebitlist(&lbm,bl);
  *bl=lbm.bl;

}

/******************************************************************************************/
/* generate random bitlist for line bundle model */

struct bitlist randomstate()
{
  struct bitlist bl;
  struct lbmodel lbm;
  extern struct envsettings set;

  lbm.lbs=randomlbs(set.rank);
  completemodel(&lbm,0);
  return lbm.bl;

}


/******************************************************************************************/
/* check if two line bundle models are equivalent */

int lbmequiv(struct lbmodel *lbm1, struct lbmodel *lbm2)
{

  extern struct CYdata CY;
  int i, ch2sum1=0, ch2sum2=0;

  for (i=0; i<CY.numps; i++) {
    ch2sum1=ch2sum1+abs((lbm1->ch2)[i]);
    ch2sum2=ch2sum2+abs((lbm2->ch2)[i]);
  }
  if (ch2sum1!=ch2sum2) return 0;
  else return lbsequiv(lbm1->lbs,lbm2->lbs);

}


/******************************************************************************************/
/* compute line bundle models from list of bitlists and remove redundancies */

struct lbmodel * removeredlbm(struct bitlist *bl, int *len)

{

  int i, cnonred, cactive, red, k;
  struct lbmodel *lbm;

  /* allocate memory */
  lbm=calloc(*len,sizeof(struct lbmodel));
  
  /* first compute full line bundle models for all bitlists */
  for (i=0; i<*len; i++) completebitlist(&(lbm[i]),&(bl[i]));

  /* then eliminate redundancies */
  if (*len>1) {
    cnonred=1; cactive=1;
    while (cactive < *len) {
      red=0; k=0;
      while (!red && k<cnonred) {
	if (lbmequiv(&(lbm[k]),&(lbm[cactive]))) red=1;
        k++;
      }
      if (!red) {
	lbm[cnonred]=lbm[cactive];
	cnonred++;}
      cactive++;
    }
    *len=cnonred;
  }

  return lbm;
  
}




  

  
