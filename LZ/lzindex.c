
#include "lzindex.h"
#include "search.c"
#include "index.c"
#include "persist.c"

        // frees lzindex

int free_index(void *index)
 { 
    lzindex I = *(lzindex *)index;
    destroyLZTrie(I.fwdtrie);
    destroyRevTrie(I.bwdtrie);
    destroyNodemap(I.map);
    destroyNodemap(I.rmap);
    destroyPosition(I.TPos);
    return 0; // no errors
 }

int index_size(void *index, ulong *size)
 {
    lzindex I = *(lzindex *)index;
    *size = sizeof(lzindex) +
           sizeofLZTrie(I.fwdtrie) +
           sizeofRevTrie(I.bwdtrie, (I.fwdtrie->n)) +
           sizeof(struct snodemap)+((I.fwdtrie->n/W)*I.map->nbits+1)
                  *sizeof(uint)+ //Node
           sizeof(struct snodemap)+((I.fwdtrie->n/W)*I.rmap->nbits+1)
                  *sizeof(uint)+ //RNode
           sizeofPosition(I.TPos); 
    return 0; // no errors
 }

 
int get_length(void *index, ulong *text_length) 
 {
    *text_length = ((lzindex *)index)->u;
    return 0; // no errors
 }

char *error_index (int e) 
 {
    printf("char *error_index(int e): Function not yet implemented...\n");
    exit(1);
 }
