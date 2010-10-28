   
 // Indexing module

#include "trie.h"
#include "lztrie.h"
#include "nodemap.h"
#include "revtrie.h"
#include "lzindex.h"
#include <math.h>

	// creates lztrie over a null-terminated text
	// it also creates *ids

#ifdef INDEXREPORT
struct tms time;
clock_t t1,t2;
uint ticks;
#endif

lztrie buildLZTrie(byte *text, byte s)
 { 
    trie T;
    uint n;
    uint *parent, *ids;
    byte *letters;
    lztrie LZT;
    unsigned long long aux;
    // first creates a full trie T
#ifdef INDEXREPORT
    ticks= sysconf(_SC_CLK_TCK);
    times(&time); t1 = time.tms_utime;
    printf("  Building LZTrie...\n"); fflush(stdout);
    printf("    Building normal trie...\n"); fflush(stdout);
#endif
    T = createTrie();
    do {
       text = insertTrie(T,text);
    }   
    while (text[-1] != s);
    // now compresses it
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf("    Representing with parentheses, letters and ids...\n"); fflush(stdout);
#endif
    n       = T->nid;
    aux = ((unsigned long long)2*n+W-1)/W;
    parent  = malloc(aux*sizeof(uint));
    letters = malloc(n*sizeof(byte));
    aux = ((unsigned long long)n*bits(n-1)+W-1)/W;
    ids    = malloc(aux*sizeof(uint));
    //malloc(n*sizeof(uint));
    representTrie(T,parent,letters,ids,NULL,bits(n-1));
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf("    Freing trie...\n"); fflush(stdout);
#endif
    destroyTrie(T);
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf("    Creating compressed trie...\n"); fflush(stdout);
#endif
    LZT = createLZTrie(parent,letters,ids,n);
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf("  End of LZTrie\n"); fflush(stdout);
#endif
    return LZT;
 }


	// builds reverse trie from LZTrie, Map, and maximum LZTrie depth
	// returns reverse ids

revtrie buildRevTrie(lztrie T, uint maxdepth, uint **rids)
 { 
    byte *str;
    uint n,rn,depth,j;
    trieNode i;
    trie RT;
    uint *parent, *emptybmap;
    revtrie CRT;
    unsigned long long aux;
    // first create a full trie RT
#ifdef INDEXREPORT
    times(&time); t1 = time.tms_utime;
    printf ("  Building RevTrie...\n"); fflush(stdout);
    printf ("    Creating full trie...\n"); fflush(stdout);
#endif
    str = malloc(maxdepth*sizeof(byte));
    RT  = createTrie(); 
    i   = ROOT; depth = 0;
    for (j=1;j<T->n;j++) { 
       i = nextLZTrie(T,i,&depth);
       str[maxdepth-depth] = letterLZTrie(T,i);
       insertstringTrie(RT,str+maxdepth-depth,depth,idLZTrie(T,i));
    }
    free(str);
    // now compresses it
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1        = t2;
    printf("    Representing with parentheses and ids...\n"); fflush(stdout);
#endif
    n         = T->n; 
    rn        = RT->nid;
    aux       = ((unsigned long long)2*rn+W-1)/W;
    parent    = malloc(aux*sizeof(uint)); // 2*rn bits
    emptybmap = calloc(((rn+W-1)/W),sizeof(uint));   // rn bits
    aux       = ((unsigned long long)n*bits(n-1)+W-1)/W;
    *rids     = malloc(aux*sizeof(uint)); // the rids array has n entries
                                         // (only for the non-empty nodes)
    representTrie(RT,parent,NULL,*rids,emptybmap,bits(n-1));
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf("    Freeing trie...\n"); fflush(stdout);
#endif
    destroyTrie(RT);
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf("    Creating compressed trie...\n"); fflush(stdout);
#endif
    CRT = createRevTrie(parent,T,emptybmap,rn);
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf("  End of RevTrie...\n"); fflush(stdout);
#endif
    return CRT;
 }

        // builds Map from LZTrie, for temporary use

uint *buildMap(lztrie T)
 {
    uint *map;
    trieNode i;
    uint n,j,depth;
#ifdef INDEXREPORT
    times(&time); t1 = time.tms_utime;
    printf ("  Building Map...\n"); fflush(stdout);
    printf ("    Computing indexes...\n"); fflush(stdout);
#endif
    n = T->n;
    map = malloc (n*sizeof(uint));
    map[0] = ROOT; depth = 0;
    i = ROOT;
    for (j = 1; j < n; j++) {
       i = nextLZTrie (T,i,&depth);
       map[rthLZTrie(T,j)] = i;
    }
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf ("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf ("  End of BuildMap...\n"); fflush(stdout);
#endif
    return map;
 }
 

 
void buildMaps(lztrie T, revtrie RT, uint *rids, nodemap *RMap, nodemap *Rev)
 { 
    uint *rmap,*rev,*map,nbits;
    trieNode i;
    uint n,rn,j;
#ifdef INDEXREPORT
    times(&time); t1 = time.tms_utime;
    printf("  Building RMap and Rev...\n"); fflush(stdout);
    printf("    Computing indexes...\n"); fflush(stdout);
#endif
    map     = buildMap(T);
    n       = T->n; // number of LZTrie nodes 
    rn      = RT->n; // number of RevTrie nodes
    rmap    = malloc(n*sizeof(uint));
    rev     = malloc(n*sizeof(uint));
    rmap[0] = ROOT; // RevTrie root
    rev[0]  = ROOT; // LZTrie root
    i       = ROOT; 
    nbits   = T->nbits;
    for (j = 1; j < n;) {
       i = nextRevTrie(RT,i);
       if (!isemptyRevTrie(RT,i)) {
          uint id = bitget(rids, j*nbits, nbits);
	  rmap[id] = i;
          rev[j] = map[id];
          j++;
       }
    }
    free(rids);
    free(map);
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf("    Creating nodemap...\n"); fflush(stdout);
#endif
    *RMap = createNodemap(rmap,n,rn);
    free (rmap);
    *Rev = createNodemap(rev,n,n);
    free(rev);
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf("  End of RMap\n"); fflush(stdout);
#endif
 }

byte selectSymbol(byte *text, ulong length)
 {
    ulong i;
    byte s;
    bool *A = calloc(256, sizeof(bool));;
    
    for (i=0;i<length;i++) A[text[i]]= true;
    for (s=0;s<256;s++)
       if (!A[s]) break;
    return s;
 } 

        // creates lzindex over a null-terminated text
	// frees text
int build_index(byte *text, ulong length, char *build_options, void **index)
 { 
    lzindex *I;
    uint *rids,maxdepth;
    
    I = malloc(sizeof(lzindex));
    text[length] = selectSymbol(text, length);
    // build index
    I->fwdtrie = buildLZTrie(text, text[length]);
    maxdepth   = maxdepthLZTrie(I->fwdtrie);
    I->bwdtrie = buildRevTrie(I->fwdtrie,maxdepth, &rids);
    buildMaps(I->fwdtrie,I->bwdtrie,rids,&I->rmap,&I->Rev);
    setRevTrieRev(I->bwdtrie,I->Rev);
    I->TPos    = createPosition(*I, length); 
    I->u       = length;
    *index = I; // return index
    return 0; // no errors yet
 }
