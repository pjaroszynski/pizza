

#include "perm.h"
#include "bitmap.h"
#include "basics.h"
#include <math.h>


int compare(const void *p1, const void *p2)
 {
    return  ((auxbwd *)p1)->key - ((auxbwd *)p2)->key;
 }


perm createPerm(uint *elems, uint nelems, uint t)
 {
    perm P;
    uint *b, *baux, nextelem, i, j, bptr, 
         aux, antbptr,nbwdptrs, elem,nbits, firstelem, cyclesize;
    auxbwd *auxbwdptr;

    P = malloc(sizeof(struct sperm));
    P->elems  = elems;
    P->nelems = nelems;
    nbits = P->nbits  = bits(nelems-1);
    P->t = t;
    if (t==1) {
       P->bwdptrs = malloc((((unsigned long long)nelems*nbits+W-1)/W)*sizeof(uint));
       P->nbwdptrs = nelems;
       for (i=0; i<nelems; i++)
          bitput(P->bwdptrs, bitget(elems, i*nbits, nbits)*nbits, nbits, i);
       P->bmap = NULL;  
    }
    else {
       auxbwdptr = malloc(sizeof(auxbwd)*(t+((int)ceil((double)nelems/t))));
       b = calloc(((nelems+W-1)/W), sizeof(uint));
       baux = calloc(((nelems+W-1)/W), sizeof(uint));
       nbwdptrs = 0; 
       for (i = 0; i < nelems; i++) {
          if (bitget1(baux,i) == 0) {
             nextelem = j = bptr = antbptr = i; 
             aux = 0;
             bitset(baux, j);
             cyclesize = 0;
             firstelem = j;
             while ((elem=bitget(elems,j*nbits,nbits)) != nextelem) {
                j = elem;
                bitset(baux, j);
                aux++;
                if (aux >= t) {
                   auxbwdptr[nbwdptrs].key = j;
                   auxbwdptr[nbwdptrs++].pointer = bptr;
                   antbptr = bptr;
                   bptr    = j;
                   aux     = 0;
                   bitset(b, j);
                }
                cyclesize++;
             }
             if (cyclesize >= t) {
                auxbwdptr[nbwdptrs].key = nextelem;
                auxbwdptr[nbwdptrs++].pointer = bptr;
                bitset(b, nextelem);
             }
          }
       }
       qsort(auxbwdptr, nbwdptrs, sizeof(auxbwd), &compare);
       aux = ((unsigned long long)nbwdptrs*P->nbits+W-1)/W;
       P->bwdptrs = malloc(sizeof(uint)*aux); 
       P->nbwdptrs = nbwdptrs;
       for (i = 0; i < nbwdptrs; i++) 
          bitput(P->bwdptrs, i*nbits, nbits, auxbwdptr[i].pointer); 
       P->bmap = createBitmap(b, nelems, false);
       free(baux);
       free(auxbwdptr);
    }
#ifdef QUERYREPORT
    P->cont_invperm = 0;
    P->cont_perm = 0;
#endif    
    return P;
 }

void destroyPerm(perm P)
 {
    free(P->elems);
    if (P->bmap) destroyBitmap(P->bmap);
    free(P->bwdptrs);
    free(P);
 }


uint nbits_perm, *e_perm;

inline uint bitget_perm(uint p)
 {
     register uint i=(p>>5), j=p&0x1F, answ;
     if (j+nbits_perm <= W)
        answ = (e_perm[i] << (W-j-nbits_perm)) >> (W-nbits_perm);
     else
        answ = (e_perm[i] >> j) | ((e_perm[i+1] << (W-j-nbits_perm)) >> (W-nbits_perm));
     return answ;
 }




// Computes P-1[i] 
uint inversePerm(perm P, uint i) 
 {
    uint j, elem, *data;
    bitmap bmap;
    
    if (P->t==1) {
       e_perm = P->bwdptrs;
       j = bitget_perm(i*nbits_perm);
    } 
    else {
       j     = i;
       e_perm = P->elems;
       bmap  = P->bmap;
       data  = bmap->data;
       while (((elem=bitget_perm(j*nbits_perm)) != i)&&(bitget1(data,j)==0)) 
          j = elem;
       
       if (elem != i) { 
          // follows the backward pointer
          j = bitget(P->bwdptrs, rank(bmap,j)*nbits_perm, nbits_perm);
	  while ((elem = bitget_perm(j*nbits_perm))!= i)  
             j = elem;
       }
    }
#ifdef QUERYREPORT
    P->cont_invperm++;
#endif
    return j;
 }


        // gets the ith element of a perm P

uint getelemPerm(perm P, uint i)
 {
#ifdef QUERYREPORT
    P->cont_perm++;
#endif
    return bitget(P->elems, i*P->nbits, P->nbits);
 }


uint savePerm(perm P, FILE *f) 
 {
    unsigned long long aux;
    uint v;
    
    if (fwrite(&P->nelems,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot write Permutation on file\n");
       exit(1);
    }

    aux = (((unsigned long long)P->nelems)*P->nbits+W-1)/W;
    if (fwrite(P->elems,sizeof(uint),aux,f) != aux) {
       fprintf(stderr,"Error: Cannot write Permutation on file\n");
       exit(1);
    }
    
    aux = ((P->nelems+W-1)/W);
    
    if (P->bmap) {
       v=1;
       if (fwrite(&v,sizeof(uint),1,f) != 1) {
          fprintf(stderr,"Error: Cannot write Permutation on file\n");
          exit(1);
       }
       if (fwrite(P->bmap->data,sizeof(uint),aux,f) != aux) {
          fprintf(stderr,"Error: Cannot write Permutation on file\n");
          exit(1);
       }
       saveBitmap(f, P->bmap);
    }
    else {
       v=0;
       if (fwrite(&v,sizeof(uint),1,f) != 1) {
          fprintf(stderr,"Error: Cannot write Permutation on file\n");
          exit(1);
       }
    }
    
    if (fwrite(&P->nbwdptrs,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot write Permutation on file\n");
       exit(1);
    }
    
    aux = ((unsigned long long)P->nbwdptrs*P->nbits+W-1)/W;
    if (fwrite(P->bwdptrs,sizeof(uint),aux,f) != aux) {
       fprintf(stderr,"Error: Cannot write Permutation on file\n");
       exit(1);
    }
    if (fwrite(&P->t,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot write Permutation on file\n");
       exit(1);       
    } 
 } 


perm loadPerm(FILE *f) 
 {
    unsigned long long aux;
    perm P;
    uint *data, v;
    
    P = malloc(sizeof(struct sperm));
    
    if (fread(&P->nelems,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot read Permutation from file\n");
       exit(1);
    }
    P->nbits = bits(P->nelems-1);
    aux = (((unsigned long long)P->nelems)*P->nbits+W-1)/W;
    P->elems = malloc(aux*sizeof(uint));
    
    if (fread(P->elems,sizeof(uint),aux,f) != aux) {
       fprintf(stderr,"Error: Cannot read Permutation from file\n");
       exit(1);
    }
    
    if (fread(&v,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot read Permutation from file\n");
       exit(1);
    }
    
    if (v) {
       aux = (P->nelems+W-1)/W;
       data = malloc(aux*sizeof(uint));
       if (fread(data,sizeof(uint),aux,f) != aux) {
          fprintf(stderr,"Error: Cannot read Permutation from file\n");
          exit(1);
       }
         
       P->bmap = loadBitmap(f,P->nelems,data);
    }
    else P->bmap = NULL;
    
    if (fread(&P->nbwdptrs,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot read Permutation from file\n");
       exit(1);
    }
    
    aux = ((unsigned long long)P->nbwdptrs*P->nbits+W-1)/W;
    P->bwdptrs = malloc(aux*sizeof(uint));
    
    if (fread(P->bwdptrs,sizeof(uint),aux,f) != aux) {
       fprintf(stderr,"Error: Cannot read Permutation from file\n");
       exit(1);
    }
    
    if (fread(&P->t,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot read Permutation from file\n");
       exit(1);
    }
    
    return P;
 } 

uint sizeofPerm(perm P)
 {
    return sizeof(struct sperm) +
           ((((unsigned long long)P->nelems)*P->nbits+W-1)/W)*sizeof(uint) +
           ((P->bmap)?sizeofBitmap(P->bmap):0) +
           (((P->nbwdptrs*P->nbits+W-1)/W)*sizeof(uint));
 }

