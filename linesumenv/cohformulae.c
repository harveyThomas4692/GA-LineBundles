
int h07862(int k1, int k2, int k3, int k4) {

    extern struct CYdata CY;
    extern const int m7862[484][4][4];
    int coh=-1;
    int kvec[MAXPS];
    int nP;
    
    nP=CY.numps;
    
    
    
    if ((k1>=0)&&(k2>=0)&&(k3>=0)&&(k4>=0)) {
        kvec[0]=k1;kvec[1]=k2;kvec[2]=k3;kvec[3]=k4;
    }
    else {
        int maxk,numMatrices;
        maxk=k1;
        if (k2>maxk) {
            maxk=k2;
        };
        if (k3>maxk) {
            maxk=k3;
        };
        if (k4>maxk) {
            maxk=k4;
        };
        
        if (maxk<6) {
            numMatrices=4;
        }
        else {
            if (maxk<12) {
                numMatrices=16;
            }
            else {
                if (maxk<20) {
                    numMatrices=52;
                }
                else {
                    if (maxk<30) {
                        numMatrices=160;
                    }
                    else {
                        numMatrices=484;
                    }
                }
            }
        };
        
        kvec[0]=k1;kvec[1]=k2;kvec[2]=k3;kvec[3]=k4;
        int ok=0;
        int i=0;
        int j,k;
        int newkvec[nP];
        
        while ((i<numMatrices) && (ok==0)) {
            ok=1;
            newkvec[0]=0;newkvec[1]=0;newkvec[2]=0;newkvec[3]=0;
            j=0;
            while ((j<nP)&&(ok==1)) {
                for (k=0; k<nP; k++) {
                    newkvec[j]=newkvec[j]+(m7862[i][j][k])*(kvec[k]);
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
/******************************************************/
/******************************************************/
/********************* 7447  **************************/
/******************************************************/
/******************************************************/
/******************************************************/


int h07447(int k1, int k2, int k3, int k4, int k5) {
    
    extern struct CYdata CY;
    extern const int m7447[1145][5][5];
    
    
    int coh=-1;
    int kvec[MAXPS];
    int nP;
    
    nP=CY.numps;
    
    
    
    if ((k1>=0)&&(k2>=0)&&(k3>=0)&&(k4>=0)&&(k5>=0)) {
        kvec[0]=k1;kvec[1]=k2;kvec[2]=k3;kvec[3]=k4;kvec[4]=k5;
    }
    else {
        int maxk,numMatrices;
        maxk=k1;
        if (k2>maxk) {
            maxk=k2;
        };
        if (k3>maxk) {
            maxk=k3;
        };
        if (k4>maxk) {
            maxk=k4;
        };
        if (k5>maxk) {
            maxk=k5;
        };
        
        if (maxk<2) {
            numMatrices=5;
        }
        else {
            if (maxk<4) {
                numMatrices=95;
            }
            else {
                if (maxk<6) {
                    numMatrices=335;
                }
                else { /* the next change at maxk = 9 */
                    numMatrices=1145;
                }
            }
        };
        
        kvec[0]=k1;kvec[1]=k2;kvec[2]=k3;kvec[3]=k4;kvec[4]=k5;
        int ok=0;
        int i=0;
        int j,k;
        int newkvec[nP];
        
        while ((i<numMatrices) && (ok==0)) {
            ok=1;
            newkvec[0]=0;newkvec[1]=0;newkvec[2]=0;newkvec[3]=0;newkvec[4]=0;
            j=0;
            while ((j<nP)&&(ok==1)) {
                for (k=0; k<nP; k++) {
                    newkvec[j]=newkvec[j]+(m7447[i][j][k])*(kvec[k]);
                };
                if (newkvec[j]<0) {
                    ok=0;
                    newkvec[0]=0;newkvec[1]=0;newkvec[2]=0;newkvec[3]=0;newkvec[4]=0;
                };
                j++;
            };
            i++;
        };
        if (ok==0) {
            coh=0;
        }
        else {
            kvec[0]=newkvec[0];kvec[1]=newkvec[1];kvec[2]=newkvec[2];kvec[3]=newkvec[3];kvec[4]=newkvec[4];
        };
    };
    
    if (coh==-1) {
        int zeroscount=0;
        int nonzeroscount=0;
        
        int zeros[nP];
        int nonzeros[nP];
        
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
        
        if (k5==0) {
            zeros[zeroscount]=k5;
            zeroscount++;
        }
        else {
            nonzeros[nonzeroscount]=k5;
            nonzeroscount++;
        };
        
        if (zeroscount>=2) {
            coh=(1+kvec[0])*(1+kvec[1])*(1+kvec[2])*(1+kvec[3])*(1+kvec[4]);
        }
        else
            coh=2*(kvec[0] + kvec[1] + kvec[2] + kvec[3] + kvec[4] + kvec[0]*kvec[3]*kvec[4] + kvec[1]*kvec[3]*kvec[4] + kvec[2]*kvec[3]*kvec[4] + kvec[0]*kvec[2]*(kvec[3] + kvec[4]) + kvec[1]*kvec[2]*(kvec[3] + kvec[4]) + kvec[0]*kvec[1]*(kvec[2] + kvec[3] + kvec[4]));
        
    };
    
    return coh;
    
};


int h17447(int k1, int k2, int k3, int k4, int k5) {
    int coh0=h07447(k1,k2,k3,k4,k5);
    int coh1;
    int coh3=h07447(-k1,-k2,-k3,-k4,-k5);
    int ind = 2*(k1 + k2 + k3 + k4 + k5 + k1*k2*(k3 + k4 + k5) + k1*k3*(k4+k5) + k1*k4*k5 + k2*k3*(k4+k5) + k2*k4*k5 + k3*k4*k5);
    
    
    int zeroscount=0;
    int nonzeroscount=0;
    int i;
    
    int zeros[5];
    int nonzeros[5];
    
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
    
    if (k5==0) {
        zeros[zeroscount]=k5;
        zeroscount++;
    }
    else {
        nonzeros[nonzeroscount]=k5;
        nonzeroscount++;
    };
    
    if (zeroscount==3) {
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
