
// Closed hash table

#include "hash.h"

  // creates a table to store up to n values with guaranteed load factor.
  // vbits = # of bits per entry, ENTRIES CANNOT HAVE ZERO VALUE

hash createHash (uint n, uint vbits, float factor)
   
   { hash H = malloc (sizeof(struct shash));
     int i,N;
     if (n == 0) return NULL;
     N = n*factor; if (N <= n) N = n+1;
     H->size = (1 << bits(N-1)) - 1;
     H->bits = vbits;
     H->table = malloc ((((H->size+1)*vbits+W-1)/W)*sizeof(uint));
#ifdef INDEXREPORT
       printf ("     Also created hash table of %i bits\n",
		  (((H->size+1)*vbits+W-1)/W)*W);
#endif
     for (i=(((H->size+1)*vbits+W-1)/W)-1;i>=0;i--) H->table[i] = 0;
     return H;
   }

  // frees the structure

void destroyHash (hash H)

   { if (H == NULL) return;
     free (H->table);
     free (H);
   }

void saveHash(FILE *f, hash H)
 {
    uint indicator;
    if (H) {
       indicator = 1;
       fwrite(&indicator, sizeof(uint), 1, f);
       fwrite(&H->size, sizeof(uint), 1, f);
       fwrite(&H->bits, sizeof(uint), 1, f);
       fwrite(H->table, sizeof(uint), ((H->size+1)*H->bits+W-1)/W, f);
    }
    else {
       indicator = 0;
       fwrite(&indicator, sizeof(uint), 1, f);
    }
 }   
   

hash loadHash(FILE *f)
 {
    hash H;
    uint indicator,i;
   
    fread(&indicator, sizeof(uint), 1, f);
    if (indicator) {
       H = malloc(sizeof(struct shash));
       fread(&H->size, sizeof(uint), 1, f);
       fread(&H->bits, sizeof(uint), 1, f);
       H->table = malloc((((H->size+1)*H->bits+W-1)/W)*sizeof(uint));
       fread(H->table, sizeof(uint), ((H->size+1)*H->bits+W-1)/W, f);      
    }
    else H = NULL;
    return H;
 }
   
  // inserts an entry, not prepared for overflow

void insertHash (hash H, uint key, uint value)

   { uint pos = (key*PRIME1) & H->size;
     if (bitget(H->table,pos*H->bits,H->bits) != 0)
	{ do pos = (pos + PRIME2) & H->size;
          while (bitget(H->table,pos*H->bits,H->bits) != 0);
	}
     bitput(H->table,pos*H->bits,H->bits,value);
   }

  // looks for a key, returns first value (zero => no values)
  // writes in pos a handle to get next values

uint searchHash (hash H, uint key, handle *h)

   { *h = (key*PRIME1) & H->size;
     return bitget(H->table,*h*H->bits,H->bits);
   }

  // gets following values using handle *pos, which is rewritten
  // returns next value (zero => no more values)

uint nextHash (hash H, handle *h)

   { *h = (*h +PRIME2) & H->size;
     return bitget(H->table,*h*H->bits,H->bits);
   }


uint sizeofHash(hash H)
 {
    if (H) {
       unsigned long long aux;
       aux = (((unsigned long long)(H->size+1))*H->bits+W-1)/W;
       return sizeof(struct shash) + aux*sizeof(uint);
    }
    else return 0;
 }

