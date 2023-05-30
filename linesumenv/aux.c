/******************************************************************************************/
/*aux.c - auxiliary functions for line sum environment - mostly imput/output */
/******************************************************************************************/


/******************************************************************************************/
/* load in global settings and CY data from file */

void confenv(FILE *fp)
{
  
  FILE *fps;
  int i, j, k, p;
  extern struct envsettings set;
  extern struct CYdata CY;
  int coh[4], lb[MAXPS];
  int ord, symmok;

  /* initialise random number generator */
  srand(time(0));
  
  /* read in global environment settings from file options */
  fps=fopen("options","r");
  fscanf(fps,"SlopeZeroMeth: %i\nRank: %i\nSymmOrder: %i\n",&(set.slopezerometh),&(set.rank),&(set.symmorder));
  fscanf(fps,"WilsonLine: %i %i\n",&((set.wilsonline)[0]),&((set.wilsonline)[1]));
  fscanf(fps,"DownStairs: %i\nMinEntry: %i\nMaxEntry: %i\nNumBits: %i\nValueWeights: %i\n",
	 &(set.downstairs),&(set.minentry),&(set.maxentry),&(set.numbits),&(set.numweights));
  for (i=0; i<set.numweights; i++) fscanf(fps,"%i ",&((set.weights)[i]));
  fscanf(fps,"ProbDis: ");
  for (i=0; i<=set.maxentry-set.minentry; i++) fscanf(fps,"%f ",&((set.probdis)[i]));
  fscanf(fps,"\nTermCond: %f",&set.termcond);
  fclose(fps);

  /* read in CY data */
  fscanf(fp,"CicyNum: %i\nNumPs: %i\nNumPols: %i\nConf:\n",&(CY.cicynum),&(CY.numps),&(CY.numpols));
  /* configuration matric */
  for (i=0; i<CY.numps; i++) {
    for (p=0; p<CY.numpols; p++) fscanf(fp,"%i",&((CY.conf)[i][p]));
    fscanf(fp,"\n");
  }
  fscanf(fp,"c2TX: ");
  for (i=0; i<CY.numps; i++) fscanf(fp,"%i",&((CY.c2TX)[i]));
  fscanf(fp,"\nSymmOrd: %i \n",&(CY.numord)); symmok=0;
  for (i=0; i<CY.numord; i++) {
    fscanf(fp,"%i ",&ord); (CY.symmord)[i]=ord;
    if (ord==set.symmorder) symmok=1;
  }
  if (!symmok) {
    set.symmorder=(CY.symmord)[0];
    printf("Warning: SymmOrd was incomatible with space. Set to %i!\n",set.symmorder);
  }
  fscanf(fp,"\nIntersecNumbers: \n");
  for (i=0; i<CY.numps; i++) {
    for (j=0; j<CY.numps; j++) {
      for (k=0; k<CY.numps; k++) fscanf(fp,"%i",&((CY.isec)[i][j][k]));
      fscanf(fp,"\n");
    }
    fscanf(fp,"\n");
  }
  fscanf(fp,"KahlerCone: %i\n",&(CY.numkahlerconegen)); 
  for (i=0; i<CY.numkahlerconegen; i++) {
    for (j=0; j<CY.numps; j++) fscanf(fp,"%i",&((CY.kahlerconegen)[i][j]));
    fscanf(fp,"\n");
  }
  fscanf(fp,"SKahlerCone: %i\n",&(CY.numskahlerconegen));
  for (i=0; i<CY.numskahlerconegen; i++) {
    for (j=0; j<CY.numps; j++) fscanf(fp,"%i",&((CY.skahlerconegen)[i][j]));
    fscanf(fp,"\n");
  }
  fscanf(fp,"CYSymmetry: %i\n",&(CY.numsymm));
  for (i=0; i<CY.numsymm; i++) {
    for (j=0; j<CY.numps; j++) {
      for (k=0; k<CY.numps; k++) fscanf(fp,"%i",&((CY.symm)[i][j][k]));
      fscanf(fp,"\n");
    }
    fscanf(fp,"\n");
  }

  fscanf(fp,"NonDiagSyms: %i\n",&(CY.nonDiagSym));
  for (i=0; i<CY.nonDiagSym; i++) {
    for (j=0; j<CY.numps; j++) {
      for (k=0; k<CY.numps; k++) fscanf(fp,"%i",&((CY.nonDiagSymGens)[i][j][k]));
      fscanf(fp,"\n");
    }
    fscanf(fp,"\n");
  }

  for (i=0; i<CY.numps; i++) lb[i]=0;
  CY.cohformula=cohline(coh,lb);
  
}

/******************************************************************************************/
/* print global settings */

void printsettings()
{

  int i;
  extern struct envsettings set;

  printf("SopeZeroMeth: %i\nRank: %i\nSymmOrder: %i\n",set.slopezerometh,set.rank,set.symmorder);
  printf("WilsonLines: %i %i\n",(set.wilsonline)[0],(set.wilsonline)[1]);
  printf("Downstairs: %i\nMinEntry: %i\nMaxEntry: %i\nNumBits: %i\nValueWeights: %i\n",
	 set.downstairs,set.minentry,set.maxentry,set.numbits,set.numweights);
  for (i=0; i<set.numweights; i++) printf("%i ",(set.weights)[i]);
  printf("\nProbDis: ");
  for (i=0; i<=set.maxentry-set.minentry; i++) printf("%4.2f ",set.probdis[i]);
  printf("\nTermCond: %f",set.termcond);
  printf("\n");
}

/******************************************************************************************/
/* print CY data */

void printCYdata()
{

  int i, j, k, p;
  extern struct CYdata CY;

  printf("CicyNum: %i\nNumPs: %i\nNumPols: %i\nConf:\n",CY.cicynum,CY.numps,CY.numpols);
  for (i=0; i<CY.numps; i++) {
    for (p=0; p<CY.numpols; p++) printf("%i ",(CY.conf)[i][p]);
    printf("\n");
  }
  printf("c2TX: ");
  for (i=0; i<CY.numps; i++) printf("%i ",(CY.c2TX)[i]);
  printf("\nSymmOrd: %i \n",CY.numord);
  for (i=0; i<CY.numord; i++) printf("%i ",(CY.symmord)[i]);
  printf("\nIntersecNumbers: \n");
  for (i=0; i<CY.numps; i++) {
    for (j=0; j<CY.numps; j++) {
      for (k=0; k<CY.numps; k++) printf("%i ",(CY.isec)[i][j][k]);
      printf("\n");
    }
    printf("\n");
  }
  printf("Kahlercone: %i\n",CY.numkahlerconegen);
  for (i=0; i<CY.numkahlerconegen; i++) {
    for (j=0; j<CY.numps; j++) printf("%i ",(CY.kahlerconegen)[i][j]);
    printf("\n");
  }
  printf("SKahlerCone: %i\n",CY.numskahlerconegen);
  for (i=0; i<CY.numskahlerconegen; i++) {
    for (j=0; j<CY.numps; j++) printf("%i ",(CY.skahlerconegen)[i][j]);
    printf("\n");
  }
  printf("CYSymmetry: %i\n",CY.numsymm);
  for (i=0; i<CY.numsymm; i++) {
    for (j=0; j<CY.numps; j++) {
      for (k=0; k<CY.numps; k++) printf("%i ",(CY.symm)[i][j][k]);
      printf("\n");
    }
    printf("\n");
  }
  printf("CohFormula: %i\n",CY.cohformula);
  
  printf("NonDiagSyms: %i\n",CY.nonDiagSym);
  for (i=0; i<CY.nonDiagSym; i++) {
    for (j=0; j<CY.numps; j++) {
      for (k=0; k<CY.numps; k++) printf("%i ",(CY.nonDiagSymGens)[i][j][k]);
      printf("\n");
    }
    printf("\n");
  }
  
}  
  

/******************************************************************************************/
/* print line bundle sum */

void printlbs(struct linebundlesum lbs)
{

  int i, a;
  extern struct CYdata CY;

  printf("{");
  for (i=0; i<CY.numps-1; i++) printlist(lbs.entries[i],lbs.rank,',');
  printlist(lbs.entries[CY.numps-1],lbs.rank,'}');
  
}
    

/******************************************************************************************/
/* write line bundle sum to file */

void fprintlbs(FILE *fp, struct linebundlesum lbs)
{

  int i, a;
  extern struct CYdata CY;

  fprintf(fp,"{");
  for (i=0; i<CY.numps-1; i++) fprintlist(fp,lbs.entries[i],lbs.rank,',');
  fprintlist(fp,lbs.entries[CY.numps-1],lbs.rank,'}');
  
}


/******************************************************************************************/
/* print a line bundle model */

void printlbm(struct lbmodel lbm)
{

  int i, a, rk2;
  extern struct envsettings set;
  extern struct CYdata CY;
  
  printf("<|\"LBS\"->");
  printlbs(lbm.lbs);
  printf(",\"BL\"->");
  printbitlist(lbm.bl,' ');
  printf(",\"ch2(V)\"->");
  printlist(lbm.ch2,CY.numps,',');
  printf("\"indLa\"->");
  printlist(lbm.indlst,lbm.lbs.rank,',');
  printf("\"ind\"->%i,\"Slope0\"->%i,\"Equiv\"->%i,\"n(OX)\"->%i,\"Splits\"->%i,\"h(L)\"->{",
	 lbm.ind,lbm.slope0,lbm.equiv,lbm.nOX,lbm.nsplits);
  for (a=0; a<lbm.lbs.rank-1; a++) printlist(lbm.cohL[a],4,',');
  printlist(lbm.cohL[lbm.lbs.rank-1],4,'}');
  printf(",\"h(V)\"->"); printlist(lbm.cohV,4,',');
  printf("\"h(LL)\"->{");
  rk2=(lbm.lbs.rank*(lbm.lbs.rank-1))/2;
  for (a=0; a<rk2-1; a++) printlist(lbm.cohLL[a],4,',');
  printlist(lbm.cohLL[rk2-1],4,'}');
  printf(",\"h(V2)\"->"); printlist(lbm.cohV2,4,',');
  printf("\"ValueLst\"->{");
  for (i=0; i<(set.numweights)-1; i++) printf("%5f,",lbm.valuelst[i]);
  printf("%5f}",lbm.valuelst[(set.numweights)-1]);
  printf("|>\n");

}


/******************************************************************************************/
/* write a line bundle model to a file */

void fprintlbm(FILE *fp, struct lbmodel lbm)
{

  int i, a, rk2;
  extern struct envsettings set;
  extern struct CYdata CY;
  
  fprintf(fp,"<|\"LBS\"->");
  fprintlbs(fp,lbm.lbs);
  fprintf(fp,",\"BL\"->");
  fprintbitlist(fp,lbm.bl,' ');
  fprintf(fp,",\"ch2(V)\"->");
  fprintlist(fp,lbm.ch2,CY.numps,',');
  fprintf(fp,"\"indLa\"->");
  fprintlist(fp,lbm.indlst,lbm.lbs.rank,',');
  fprintf(fp,"\"ind\"->%i,\"Slope0\"->%i,\"Equiv\"->%i,\"n(OX)\"->%i,\"Splits\"->%i,\"h(L)\"->{",
	 lbm.ind,lbm.slope0,lbm.equiv,lbm.nOX,lbm.nsplits);
  for (a=0; a<lbm.lbs.rank-1; a++) fprintlist(fp,lbm.cohL[a],4,',');
  fprintlist(fp,lbm.cohL[lbm.lbs.rank-1],4,'}');
  fprintf(fp,",\"cohV\"->"); fprintlist(fp,lbm.cohV,4,',');
  fprintf(fp,"\"cohLL\"->{");
  rk2=(lbm.lbs.rank*(lbm.lbs.rank-1))/2;
  for (a=0; a<rk2-1; a++) fprintlist(fp,lbm.cohLL[a],4,',');
  fprintlist(fp,lbm.cohLL[rk2-1],4,'}');
  fprintf(fp,",\"cohV2\"->"); fprintlist(fp,lbm.cohV2,4,',');
  fprintf(fp,"\"ValueLst\"->{");
  for (i=0; i<(set.numweights)-1; i++) fprintf(fp,"%5f,",lbm.valuelst[i]);
  fprintf(fp,"%5f}",lbm.valuelst[(set.numweights)-1]);
  fprintf(fp,"|>\n");

}  
