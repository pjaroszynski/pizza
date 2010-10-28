
// Implements the Map data structure, which maps block ids to lztrie nodes

#include "nodemap.h"

	// creates a nodemap structure from a mapping array, not owning it
	// n is number of blocks
	// max is the number of trie nodes

nodemap createNodemap(uint *map, uint n, uint max)
 { 
    nodemap M;
    uint i;
    unsigned long long aux;
    M = malloc(sizeof(struct snodemap));
    M->nbits = bits(2*max-1);
    aux    = ((unsigned long long)n*M->nbits+W-1)/W;
    M->map = malloc(aux*sizeof(uint));
    M->n   = n;
    for (i=0;i<n;i++)
       bitput(M->map,i*M->nbits,M->nbits,map[i]);
    return M;
 }

	// frees revtrie structure, including the owned data

void destroyNodemap(nodemap M)
 { 
    free(M->map);
    free(M);
 }

	// mapping

trieNode mapto(nodemap M, uint id)
 { 
    return bitget(M->map,id*M->nbits,M->nbits);
 }

         // saves nodemap to file f

void saveNodemap(nodemap M, FILE *f)
 { 
    unsigned long long aux;
    if (fwrite(&M->n,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot write nodemap on file\n");
       exit(1);
    }
    if (fwrite(&M->nbits,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot write nodemap on file\n");
       exit(1);
    }
    aux = ((unsigned long long)M->n*M->nbits+W-1)/W;
    if (fwrite(M->map,sizeof(uint),aux,f) != aux) {
       fprintf(stderr,"Error: Cannot write nodemap on file\n");
       exit(1);
    }
 }

        // loads nodemap from f

nodemap loadNodemap(FILE *f)
 { 
    nodemap M = malloc (sizeof(struct snodemap));
    unsigned long long aux;
    if (fread(&M->n,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot read nodemap from file\n");
       exit(1);
    }
    if (fread(&M->nbits,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot read nodemap from file\n");
       exit(1);
    }
    aux = ((unsigned long long)M->n*M->nbits+W-1)/W;
    M->map = malloc(aux*sizeof(uint));
    if (fread(M->map,sizeof(uint),aux,f) != aux) {
       fprintf(stderr,"Error: Cannot read nodemap from file\n");
       exit(1);
    }
    return M;
 }
