 
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
    if (fclose(f) != 0) { 
       fprintf(stderr,"Error: cannot write file %s\n",fnamext);
       exit(1);
    }

    return 0;
 }

extern uint nbits_perm;

        //loads index from filename fname.*

int load_index(char *filename, void **index) 
 { 
    lzindex *I;
    char fnamext[1024];
    FILE *f;
    uint depth,i,j,n;
    long pos;
    
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
    I->TPos = loadPosition(f, I->u);
    fclose(f);
    n = I->fwdtrie->n;
    nbits_perm = I->fwdtrie->ids->nbits;
    // load revtrie
    sprintf(fnamext,"%s.rvt",filename);
    f = fopen(fnamext,"r");
    if (f == NULL) { 
       fprintf(stderr,"Error: cannot open file %s\n",fnamext);
       exit(1);
    }
    I->bwdtrie = loadRevTrie(f,I->fwdtrie);    
    fclose(f);
    *index = I;
    return 0;
 }
 

