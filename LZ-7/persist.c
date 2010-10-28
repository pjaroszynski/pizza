 
// Implements the persistence of the LZIndex

        // writes index to filename fname.*

int save_index(void *index, char *filename)
 { 
    lzindex I = *(lzindex *)index;
    char fnamext[1024];
    FILE *f;
    // save lztrie
    sprintf(fnamext,"%s.lzt",filename);
    f = fopen(fnamext,"w");
    if (f == NULL) { 
       fprintf(stderr,"Error: cannot create file %s\n",fnamext);
       exit(1);
    }
    fwrite(&I.u, sizeof(ulong), 1,f); // writes the text length
    saveLZTrie(I.fwdtrie,f);
    savePosition(f, I.TPos);
    if (fclose(f) != 0) { 
       fprintf(stderr,"Error: cannot write file %s\n",fnamext);
       exit(1);
    }
    // save revtrie
    sprintf(fnamext,"%s.rvt",filename);
    f = fopen(fnamext,"w");
    if (f == NULL) { 
       fprintf(stderr,"Error: cannot create file %s\n",fnamext);
       exit(1);
    }
    saveRevTrie(I.bwdtrie,f);
    // save mapping R
    if (fwrite(I.R,sizeof(uint),(I.fwdtrie->n*I.fwdtrie->nbits+W-1)/W,f) !=
              (I.fwdtrie->n*I.fwdtrie->nbits+W-1)/W) {
       fprintf(stderr,"Error: Cannot write R mapping on file\n");
       exit(1);
    }
    if (fclose(f) != 0) { 
       fprintf(stderr,"Error: cannot write file %s\n",fnamext);
       exit(1);
    }

    return 0;
 }

        //loads index from filename fname.*

extern uint nbits_SB, nbits_Offs;

int load_index(char *filename, void **index) 
 { 
    lzindex *I;
    char fnamext[1024];
    FILE *f;
    uint n; 
    //long pos;
    unsigned long long aux;
    
    I = malloc(sizeof(lzindex));
    // load lztrie
    sprintf(fnamext,"%s.lzt",filename);
    f = fopen(fnamext,"r");
    if (f == NULL) { 
       fprintf(stderr,"Error: cannot open file %s\n",fnamext);
       exit(1);
    }
    fread(&I->u, sizeof(ulong), 1,f); // reads the text length
    I->fwdtrie = loadLZTrie(f);
    I->TPos = loadPosition(f,I->u); 
    nbits_SB = I->TPos->nbitsSB;
    nbits_Offs = I->TPos->nbitsOffs;
    fclose(f);
    n = I->fwdtrie->n;
    // load revtrie
    sprintf(fnamext,"%s.rvt",filename);
    f = fopen(fnamext,"r");
    if (f == NULL) { 
       fprintf(stderr,"Error: cannot open file %s\n",fnamext);
       exit(1);
    }
    I->bwdtrie = loadRevTrie(f,I->fwdtrie);    
    // load the R mapping
    aux = (((unsigned long long)I->fwdtrie->n*bits(I->fwdtrie->n-1)+W-1)/W);
    I->R = malloc(aux*sizeof(uint));
    if (fread(I->R,sizeof(uint),aux,f) != aux) {
       fprintf(stderr,"Error: Cannot read R mapping from file\n");
       exit(1);
    }
    
    I->fwdtrie->R = I->R;
    I->bwdtrie->R = I->R; 
   
    fclose(f);
    
    *index = I;
    return 0;
 }
 

