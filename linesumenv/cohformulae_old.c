
int h07862(int k1, int k2, int k3, int k4) {

  extern struct CYdata CY;
  int m[53][MAXPS][MAXPS]={{{-3, -2, 0, 0}, {4, 3, 0, 0}, {12, 6, 1, 0}, {12, 6, 0, 1}},
        {{-3, 0, -2, 0}, {12, 1, 6, 0}, {4, 0, 3, 0}, {12, 0, 6, 1}},
        {{-3, 0, 0, -2}, {12, 1, 0, 6}, {12, 0, 1, 6}, {4, 0, 0, 3}},
        {{-1, -6, -2, 0}, {2, 15, 6, 0}, {2, 10, 3, 0}, {2, 18, 6, 1}},
        {{-1, -6, 0, -2}, {2, 15, 0, 6}, {2, 18, 1, 6}, {2, 10, 0, 3}},
        {{-1, -2, -6, 0}, {2, 3, 10, 0}, {2, 6, 15, 0}, {2, 6, 18, 1}},
        {{-1, -2, 0, -6}, {2, 3, 0, 10}, {2, 6, 1, 18}, {2, 6, 0, 15}},
        {{-1, -2, 0, 0}, {2, 3, 0, 0}, {2, 6, 1, 0}, {2, 6, 0, 1}},
        {{-1, 0, -6, -2}, {2, 1, 18, 6}, {2, 0, 15, 6}, {2, 0, 10, 3}},
        {{-1, 0, -2, -6}, {2, 1, 6, 18}, {2, 0, 3, 10}, {2, 0, 6, 15}},
        {{-1, 0, -2, 0}, {2, 1, 6, 0}, {2, 0, 3, 0}, {2, 0, 6, 1}},
        {{-1, 0, 0, -2}, {2, 1, 0, 6}, {2, 0, 1, 6}, {2, 0, 0, 3}},
        {{-1, 0, 0, 0}, {2, 1, 0, 0}, {2, 0, 1, 0}, {2, 0, 0, 1}},
        {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}},
        {{1, 0, 0, 2}, {0, 1, 0, 2}, {0, 0, 1, 2}, {0, 0, 0, -1}},
        {{1, 0, 2, 0}, {0, 1, 2, 0}, {0, 0, -1, 0}, {0, 0, 2, 1}},
        {{1, 0, 2, 6}, {0, 1, 2, 6}, {0, 0, -1, -2}, {0, 0, 2, 3}},
        {{1, 0, 6, 2}, {0, 1, 6, 2}, {0, 0, 3, 2}, {0, 0, -2, -1}},
        {{1, 0, 6, 12}, {0, 1, 6, 12}, {0, 0, 3, 4}, {0, 0, -2, -3}},
        {{1, 0, 12, 6}, {0, 1, 12, 6}, {0, 0, -3, -2}, {0, 0, 4, 3}},
        {{1, 2, 0, 0}, {0, -1, 0, 0}, {0, 2, 1, 0}, {0, 2, 0, 1}},
        {{1, 2, 0, 6}, {0, -1, 0, -2}, {0, 2, 1, 6}, {0, 2, 0, 3}},
        {{1, 2, 6, 0}, {0, -1, -2, 0}, {0, 2, 3, 0}, {0, 2, 6, 1}},
        {{1, 2, 6, 18}, {0, -1, -2, -6}, {0, 2, 3, 10}, {0, 2, 6, 15}},
        {{1, 2, 18, 6}, {0, -1, -6, -2}, {0, 2, 15, 6}, {0, 2, 10, 3}},
        {{1, 6, 0, 2}, {0, 3, 0, 2}, {0, 6, 1, 2}, {0, -2, 0, -1}},
        {{1, 6, 0, 12}, {0, 3, 0, 4}, {0, 6, 1, 12}, {0, -2, 0, -3}},
        {{1, 6, 2, 0}, {0, 3, 2, 0}, {0, -2, -1, 0}, {0, 6, 2, 1}},
        {{1, 6, 2, 18}, {0, 3, 2, 10}, {0, -2, -1, -6}, {0, 6, 2, 15}},
        {{1, 6, 12, 0}, {0, 3, 4, 0}, {0, -2, -3, 0}, {0, 6, 12, 1}},
        {{1, 6, 18, 2}, {0, 3, 10, 2}, {0, 6, 15, 2}, {0, -2, -6, -1}},
        {{1, 12, 0, 6}, {0, -3, 0, -2}, {0, 12, 1, 6}, {0, 4, 0, 3}},
        {{1, 12, 6, 0}, {0, -3, -2, 0}, {0, 4, 3, 0}, {0, 12, 6, 1}},
        {{1, 18, 2, 6}, {0, 15, 2, 6}, {0, -6, -1, -2}, {0, 10, 2, 3}},
        {{1, 18, 6, 2}, {0, 15, 6, 2}, {0, 10, 3, 2}, {0, -6, -2, -1}},
        {{3, 0, 0, 2}, {6, 1, 0, 2}, {6, 0, 1, 2}, {-2, 0, 0, -1}},
        {{3, 0, 0, 4}, {6, 1, 0, 12}, {6, 0, 1, 12}, {-2, 0, 0, -3}},
        {{3, 0, 2, 0}, {6, 1, 2, 0}, {-2, 0, -1, 0}, {6, 0, 2, 1}},
        {{3, 0, 2, 10}, {6, 1, 2, 18}, {-2, 0, -1, -6}, {6, 0, 2, 15}},
        {{3, 0, 4, 0}, {6, 1, 12, 0}, {-2, 0, -3, 0}, {6, 0, 12, 1}},
        {{3, 0, 10, 2}, {6, 1, 18, 2}, {6, 0, 15, 2}, {-2, 0, -6, -1}},
        {{3, 2, 0, 0}, {-2, -1, 0, 0}, {6, 2, 1, 0}, {6, 2, 0, 1}},
        {{3, 2, 0, 10}, {-2, -1, 0, -6}, {6, 2, 1, 18}, {6, 2, 0, 15}},
        {{3, 2, 10, 0}, {-2, -1, -6, 0}, {6, 2, 15, 0}, {6, 2, 18, 1}},
        {{3, 4, 0, 0}, {-2, -3, 0, 0}, {6, 12, 1, 0}, {6, 12, 0, 1}},
        {{3, 10, 0, 2}, {6, 15, 0, 2}, {6, 18, 1, 2}, {-2, -6, 0, -1}},
        {{3, 10, 2, 0}, {6, 15, 2, 0}, {-2, -6, -1, 0}, {6, 18, 2, 1}},
        {{15, 0, 2, 6}, {18, 1, 2, 6}, {-6, 0, -1, -2}, {10, 0, 2, 3}},
        {{15, 0, 6, 2}, {18, 1, 6, 2}, {10, 0, 3, 2}, {-6, 0, -2, -1}},
        {{15, 2, 0, 6}, {-6, -1, 0, -2}, {18, 2, 1, 6}, {10, 2, 0, 3}},
        {{15, 2, 6, 0}, {-6, -1, -2, 0}, {10, 2, 3, 0}, {18, 2, 6, 1}},
        {{15, 6, 0, 2}, {10, 3, 0, 2}, {18, 6, 1, 2}, {-6, -2, 0, -1}},
        {{15, 6, 2, 0}, {10, 3, 2, 0}, {-6, -2, -1, 0}, {18, 6, 2, 1}}};
    
    
    int coh=-1;
    int kvec[MAXPS];
    int nP;

    nP=CY.numps;
    
    if ((k1>=0)&&(k2>=0)&&(k3>=0)&&(k4>=0)) {
        kvec[0]=k1;kvec[1]=k2;kvec[2]=k3;kvec[3]=k4;
    }
    else {
        kvec[0]=k1;kvec[1]=k2;kvec[2]=k3;kvec[3]=k4;
        int ok=0;
        int i=0;
        int j,k;
        int newkvec[nP];
        
        int numMatrices=53;
        while ((i<numMatrices) && (ok==0)) {
            ok=1;
            newkvec[0]=0;newkvec[1]=0;newkvec[2]=0;newkvec[3]=0;
            j=0;
            while ((j<nP)&&(ok==1)) {
                for (k=0; k<nP; k++) {
                    newkvec[j]=newkvec[j]+(m[i][j][k])*(kvec[k]);
                };
                if (newkvec[j]<0) {
                    ok=0;
                    newkvec[0]=0;newkvec[1]=0;newkvec[2]=0;newkvec[3]=0;
                };
                j++;
            };
            i++;
        };
        if (ok==0) {
            coh=0;
        }
        else {
            kvec[0]=newkvec[0];kvec[1]=newkvec[1];kvec[2]=newkvec[2];kvec[3]=newkvec[3];
        };
    };
    
    if (coh==-1) {
        if ((kvec[0]>=0) && (kvec[1]>=0) && (kvec[2]==0) && (kvec[3]==0)) {
            coh = (1 + kvec[0])*(1 + kvec[1]);
        }
        else {
            if ((kvec[0]>=0) && (kvec[1]==0) && (kvec[2]>=0) && (kvec[3]==0)) {
                coh = (1 + kvec[0])*(1 + kvec[2]);
            }
            else {
                if ((kvec[0]>=0) && (kvec[1]==0) && (kvec[2]==0) && (kvec[3]>=0)) {
                    coh = (1 + kvec[0])*(1 + kvec[3]);
                }
                else {
                    if ((kvec[0]==0) && (kvec[1]>=0) && (kvec[2]>=0) && (kvec[3]==0)) {
                        coh = (1 + kvec[1])*(1 + kvec[2]);
                    }
                    else {
                        if ((kvec[0]==0) && (kvec[1]>=0) && (kvec[2]==0) && (kvec[3]>=0)) {
                            coh = (1 + kvec[1])*(1 + kvec[3]);
                        }
                        else {
                            if ((kvec[0]==0) && (kvec[1]==0) && (kvec[2]>=0) && (kvec[3]>=0)) {
                                coh = (1 + kvec[2])*(1 + kvec[3]);
                            }
                            else
                                coh = 2*(kvec[0] + kvec[1] + kvec[2] + kvec[3] + kvec[0]*kvec[1]*kvec[2] + kvec[0]*kvec[1]*kvec[3] + kvec[0]*kvec[2]*kvec[3] + kvec[1]*kvec[2]*kvec[3]);
                        }
                    }
                }
            }
        }
    };
    
    
    return coh;
    
};


int h17862(int k1, int k2, int k3, int k4) {
    int coh0=h07862(k1,k2,k3,k4);
    int coh1;
    int coh3=h07862(-k1,-k2,-k3,-k4);
    int ind = 2*(k1 + k2 + k3 + k4 + k1*k2*k3 + k1*k2*k4 + k1*k3*k4 + k2*k3*k4);
    
    int zeroscount=0;
    int nonzeroscount=0;
    int i, nP;
    extern struct CYdata CY;
    
    int zeros[MAXPS];
    int nonzeros[MAXPS];

    nP=CY.numps;
    
    if (k1==0) {
        zeros[zeroscount]=k1;
        zeroscount++;
    }
    else {
        nonzeros[nonzeroscount]=k1;
        nonzeroscount++;
    };
    
    if (k2==0) {
        zeros[zeroscount]=k2;
        zeroscount++;
    }
    else {
        nonzeros[nonzeroscount]=k2;
        nonzeroscount++;
    };
    
    if (k3==0) {
        zeros[zeroscount]=k3;
        zeroscount++;
    }
    else {
        nonzeros[nonzeroscount]=k3;
        nonzeroscount++;
    };
    
    if (k4==0) {
        zeros[zeroscount]=k4;
        zeroscount++;
    }
    else {
        nonzeros[nonzeroscount]=k4;
        nonzeroscount++;
    };
    
    
    if (zeroscount==2) {
        if ((nonzeros[0]*nonzeros[1]<0)&&(abs(nonzeros[0])>1)&&(abs(nonzeros[1])>1)) {
            coh1 = (coh0 - coh3 - ind)/2 - (1 + nonzeros[0]*nonzeros[1]);
        }
        else coh1 = coh0 - coh3 - ind;
    }
    else coh1 = coh0 - coh3 - ind;

    if (coh1>0) {
        return coh1;
    }
    else return 0;

};


/******************************************************/
