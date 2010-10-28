
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

lztrie buildLZTrie(byte *text, uint **ids, byte s)
 { 
    trie T;
    uint n;
    uint *parent;
    byte *letters;
    lztrie LZT;
    // first creates a full trie T
#ifdef INDEXREPORT
    ticks= sysconf(_SC_CLK_TCK);
    times(&time); t1 = time.tms_utime;
    printf ("  Building LZTrie...\n"); fflush(stdout);
    printf ("    Building normal trie...\n"); fflush(stdout);
#endif
    T = createTrie();
    do {
       text = insertTrie(T,text);
    }   
    while (text[-1]!=s);
    // now compresses it
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf ("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf ("    Representing with parentheses, letters and ids...\n"); fflush(stdout);
#endif
    n           = T->nid;
    parent = malloc (((2*n+W-1)/W)*sizeof(uint));
    letters = malloc (n*sizeof(byte));
    *ids = malloc (n*sizeof(uint));
    representTrie (T,parent,letters,*ids);
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf ("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf ("    Freing trie...\n"); fflush(stdout);
#endif
    destroyTrie(T);
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf ("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf ("    Creating compressed trie...\n"); fflush(stdout);
#endif
    LZT = createLZTrie (parent,letters,*ids,n);
#ifdef INDEXREPORT
    times(&time); t2 = time.tms_utime;
    printf ("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
    t1 = t2;
    printf ("  End of LZTrie\n"); fflush(stdout);
#endif
    return LZT;
 }

	// builds Map from LZTrie and ids, which gets freed
	// it also writes the maximum depth of the trie

nodemap buildMap (lztrie T, uint *ids, uint *maxdepth)

   { nodemap M;
     uint *map;
     trieNode i;
     uint n,j,depth,mdepth;
#ifdef INDEXREPORT
     times(&time); t1 = time.tms_utime;
     printf ("  Building Map...\n"); fflush(stdout);
     printf ("    Computing indexes...\n"); fflush(stdout);
#endif
     n = T->n;
     map = malloc (n*sizeof(uint));
     map[0] = ROOT; depth = mdepth = 0;
     i = ROOT;
     for (j=1;j<n;j++)
	{ i = nextLZTrie (T,i,&depth);
	  if (depth > mdepth) mdepth = depth;
	  map[ids[j]] = i;
	}
     free (ids);
#ifdef INDEXREPORT
     times(&time); t2 = time.tms_utime;
     printf ("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
     t1 = t2;
     printf ("    Creating nodemap...\n"); fflush(stdout);
#endif
     M = createNodemap (map,n,n);
     free (map);
     *maxdepth = mdepth;
#ifdef INDEXREPORT
     times(&time); t2 = time.tms_utime;
     printf ("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
     t1 = t2;
     printf ("  End of Map\n"); fflush(stdout);
#endif
     return M;
   }

	// builds reverse trie from LZTrie, Map, and maximum LZTrie depth
	// returns reverse ids

revtrie buildRevTrie (lztrie T, nodemap M, uint maxdepth, uint **ids)

   { byte *str;
     uint n,depth,j;
     trieNode i;
     trie RT;
     uint *parent, *emptybmap;
     revtrie CRT;
	// first create a full trie RT
#ifdef INDEXREPORT
     times(&time); t1 = time.tms_utime;
     printf ("  Building RevTrie...\n"); fflush(stdout);
     printf ("    Creating full trie...\n"); fflush(stdout);
#endif
     str = malloc (maxdepth*sizeof(byte));
     RT = createTrie(); 
     i = ROOT; depth = 0;
     for (j=1;j<T->n;j++)
	{ i = nextLZTrie (T,i,&depth);
	  str[maxdepth-depth] = letterLZTrie (T,i);
          insertstringTrie (RT,str+maxdepth-depth,depth,idLZTrie(T,i));
	}
     free (str);
        // now compresses it
#ifdef INDEXREPORT
     times(&time); t2 = time.tms_utime;
     printf ("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
     t1 = t2;
     printf ("    Representing with parentheses and ids...\n"); fflush(stdout);
#endif 
     n = RT->nid;
     parent = malloc (((2*n+W-1)/W)*sizeof(uint)); 
     *ids = malloc (n*sizeof(uint));              
     representTrie (RT,parent,NULL,*ids);
#ifdef INDEXREPORT
     times(&time); t2 = time.tms_utime;
     printf ("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
     t1 = t2;
     printf ("    Freeing trie...\n"); fflush(stdout);
#endif
     destroyTrie(RT);
#ifdef INDEXREPORT
     times(&time); t2 = time.tms_utime;
     printf ("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
     t1 = t2;
     printf ("    Creating compressed trie...\n"); fflush(stdout);
#endif
     CRT = createRevTrie(parent,T,M,*ids,n);
#ifdef INDEXREPORT
     times(&time); t2 = time.tms_utime;
     printf ("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
     t1 = t2;
     printf ("  End of RevTrie...\n"); fflush(stdout);
#endif
     return CRT;
   }

	// builds Map from RevTrie and ids, which gets freed

nodemap buildRMap (revtrie T, uint *ids)

   { nodemap M;
     uint *map;
     trieNode i;
     uint j, n;
#ifdef INDEXREPORT
     times(&time); t1 = time.tms_utime;
     printf ("  Building RMap...\n"); fflush(stdout);
     printf ("    Computing indexes...\n"); fflush(stdout);
#endif
     n = T->n;
     map = malloc (n*sizeof(uint));
     map[0] = ROOT;
     i = ROOT;
     for (j=1;j<n;j++)
        { i = nextRevTrie (T,i);
          map[ids[j]] = i; // when equality, the innermost gets the mapping
        }
     free (ids);
#ifdef INDEXREPORT
     times(&time); t2 = time.tms_utime;
     printf ("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
     t1 = t2;
     printf ("    Creating nodemap...\n"); fflush(stdout);
#endif
     M = createNodemap (map,n,n);
     free (map);
#ifdef INDEXREPORT
     times(&time); t2 = time.tms_utime;
     printf ("    User time: %f secs\n",(t2-t1)/(float)ticks); fflush(stdout);
     t1 = t2;
     printf ("  End of Map\n"); fflush(stdout);
#endif
     return M;
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

int build_index(byte *text, ulong length, char *build_options, void **index)
 { 
    lzindex *I;
    uint *ids,maxdepth;
    
    I = malloc(sizeof(lzindex));
    text[length] = selectSymbol(text, length);
    // build index
    I->fwdtrie = buildLZTrie(text,&ids,text[length]);
    I->map     = buildMap(I->fwdtrie,ids,&maxdepth);
    I->bwdtrie = buildRevTrie(I->fwdtrie,I->map,maxdepth,&ids);
    I->rmap    = buildRMap(I->bwdtrie,ids);
    I->TPos    = createPosition(I->fwdtrie, length, I->map); 
    I->u       = length;
    *index = I; // return index
    return 0; // no errors yet
 }

