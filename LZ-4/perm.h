
#ifndef PERMINCLUDED
#define PERMINCLUDED

#include "basics.h"
#include "bitmap.h"


typedef struct sperm {
   uint *elems;   // elements of the permutation
   uint nelems;   // # of elements
   bitmap bmap;   // bitmap allowing rank() queries in O(1) time
   uint *bwdptrs; // array of backward pointers
   uint nbits;    // log(nelems)
   uint nbwdptrs; // # of backward pointers
   uint t;
#ifdef QUERYREPORT
   uint cont_perm;
   uint cont_invperm;
#endif
} *perm;

typedef struct {
   uint key;
   uint pointer;
} auxbwd;
 
perm createPerm(uint *elems, uint nelems, uint t);

uint getelemPerm(perm P, uint i);

void destroyPerm(perm P);

uint inversePerm(perm P, uint i);

uint savePerm(perm P, FILE *f);

perm loadPerm(FILE *f);

uint sizeofPerm(perm P);

#endif
