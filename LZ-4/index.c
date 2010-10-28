
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

extern uint PARAMETER_T_IDS, PARAMETER_T_RIDS;


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
    aux     = (2*((unsigned long long)n)+W-1)/W;
    parent  = malloc(aux*sizeof(uint));
    letters = malloc(n*sizeof(byte));
    aux     = (((unsigned long long)n)*bits(n-1)+W-1)/W;
    ids     = malloc(aux*sizeof(uint));
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

revtrie buildRevTrie(lztrie T, uint maxdepth)
 { 
    byte *str;
    uint n,rn,depth,j;
    trieNode i;
    trie RT;
    uint *parent, *emptybmap, *ids;
    revtrie CRT;
    unsigned long long aux;
    // first create a full trie RT
#ifdef INDEXREPORT
    times(&time); t1 = time.tms_utime;
    printf ("  Building RevTrie...\n"); fflush(stdout);
    printf ("    Creating full trie...\n"); fflush(stdout);
#endif
    str = malloc(maxdepth*sizeof(byte));
    RT = createTrie(); 
    i = ROOT; depth = 0;
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
    aux       = (2*(unsigned long long)rn+W-1)/W;
    parent    = malloc(aux*sizeof(uint)); // 2*rn bits
    emptybmap = calloc(((rn+W-1)/W),sizeof(uint));   // rn bits
    aux       = (((unsigned long long)n)*bits(n-1)+W-1)/W;
    ids       = malloc(aux*sizeof(uint)); // the rids array has n entries
                                         // (only for the non-empty nodes)
    representTrie(RT,parent,NULL,ids,emptybmap,bits(n-1));
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
    CRT = createRevTrie(parent,ids,T,emptybmap,rn);
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf("  End of RevTrie...\n"); fflush(stdout);
#endif
    return CRT;
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

extern uint nbits_perm;

        // creates lzindex over a null-terminated text
	// frees text
int build_index(byte *text, ulong length, char *build_options, void **index)
 { 
    lzindex *I;
    uint *ids,maxdepth;
    char RIDS[15], IDS[15];
    uint i, j;
    
    // set default values for the parameters
    PARAMETER_T_IDS  = 4;
    PARAMETER_T_RIDS = 4;
    
    if (build_options)
       PARAMETER_T_IDS = PARAMETER_T_RIDS = atoi(build_options);       
    I = malloc(sizeof(lzindex));
    text[length] = selectSymbol(text, length);
    // build index
    I->fwdtrie = buildLZTrie(text, text[length]);
    nbits_perm = I->fwdtrie->ids->nbits;
    maxdepth   = maxdepthLZTrie(I->fwdtrie);
    I->bwdtrie = buildRevTrie(I->fwdtrie,maxdepth);
    I->TPos    = createPosition(I->fwdtrie, length); 
    I->u       = length;
    *index = I; // return index
    return 0; // no errors yet
 }

