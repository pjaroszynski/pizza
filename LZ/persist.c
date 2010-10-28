 
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

        //loads index from filename fname.*

int load_index(char *filename, void **index) 
 { 
    lzindex *I;
    char fnamext[1024];
    FILE *f;
    uint depth,*map,i,j,n;
    long pos;
    
    I = malloc(sizeof(lzindex));
    // load lztrie
    sprintf(fnamext,"%s.lzt",filename);
    f = fopen(fnamext,"r");
    if (f == NULL) { 
       fprintf(stderr,"Error: cannot open file %s\n",fnamext);
       exit(1);
    }
    fread(&I->u, sizeof(ulong), 1,f); // writes the text length
    I->fwdtrie = loadLZTrie(f);
    fclose(f);
    // produce map
    n = I->fwdtrie->n;
    map = malloc(n*sizeof(uint));
    map[0] = ROOT;
    i = ROOT;
    for (j=1;j<n;j++) { 
       i = nextLZTrie(I->fwdtrie,i,&depth);
       map[rthLZTrie(I->fwdtrie,j)] = i;
    }
    I->map = createNodemap(map,n,n);
    free (map);
    // load revtrie
    sprintf(fnamext,"%s.rvt",filename);
    f = fopen(fnamext,"r");
    if (f == NULL) { 
       fprintf(stderr,"Error: cannot open file %s\n",fnamext);
       exit(1);
    }
    I->bwdtrie = loadRevTrie(f,I->fwdtrie,I->map);
    fclose (f);
    // produce rmap
    n = I->bwdtrie->n;
    map = malloc (n*sizeof(uint));
    map[0] = ROOT;
    i = ROOT;
    for (j=1;j<n;j++)
       { i = nextRevTrie (I->bwdtrie,i);
         map[rthRevTrie(I->bwdtrie,j)] = i;
                               // when equality, the innermost gets the mapping
       }
    I->rmap = createNodemap (map,n,n);
    free (map);
    I->TPos = createPosition(I->fwdtrie,I->u,I->map);
    *index = I;
    return 0;
 }
 

