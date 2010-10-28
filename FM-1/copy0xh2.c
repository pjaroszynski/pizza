/*  
 * algoritmo seward-like modificato per il calcolo
 * del vero suffix array (non quello ciclico).
 * 
 * ver 1.0 (12-apr-00)  
 * ver 1.1 (26-mag-00) usabile anche per file contenenti 0
 * ver 1.2 (22-giu-00) usa una hash table per i long matches 
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "sacopy.h"

extern int Verbose;

// ---------------- typedef's 
typedef int          Int32;
typedef unsigned int UInt32;
typedef short          Int16;
typedef unsigned short UInt16;
typedef char          Char;
typedef unsigned char UChar;
typedef unsigned char Bool;

// -------------------- global variables 
Int32  nBlock;
UChar  *block;   // [M_BLOCK + N_OVERSHOOT];
Int32  *ptr;       // [M_BLOCK];
Int32  ftab [65537];
Int32  runningOrder[256];

void   calc_running_order ( void );

// ------------------ constants and macros
#define True   ((Bool)1)
#define False  ((Bool)0)
#define BIGFREQ(b) (ftab[((b)+1) << 8] - ftab[(b) << 8])
#define SETMASK (1 << 30)
#define CLEARMASK (~(SETMASK))

static Bool local_debug = False;

static __inline__ Int32 cmp_aux14(Int32, Int32);

/* ---------------------------------------------------------------- */
Int32 Hash_table_diff[HASH_SIZE];
Int32 Hash_table_start[HASH_SIZE];
Int32 Hash_table_lenNris[HASH_SIZE];

Int32 Num_lm; 
Int32 Curr_len;


int *suffixsort5n(UChar *x, int n)
{
  void d_copyEQ_u12( void );
  Int32 i;

  block=x;
  nBlock=n;
  ptr = malloc(n*sizeof(int *));
  if(ptr==NULL) {
    fprintf(stderr,"malloc failed -suffix_sort-\n");
    exit(1);
  }
  // ---- init hash table with long matches
  Num_lm=0;
  for(i=0;i<HASH_SIZE;i++)
    Hash_table_diff[i] = 0;
  d_copyEQ_u12();  
  // fprintf(stderr,"Size of Lm_list: %d\n",Num_lm);
  return ptr;
}


/* ---------------------------------------------------------------
   compare the two strings starting at
   b1 and b2. If the comparison ends because a mismatch is found,
   returns + or - and stores in Curr_len the number of matching chars.
   If limit chars are compared returns 0
   --------------------------------------------------------------- */
#if (N_OVERSHOOT<16)
#error
#endif
static __inline__
Int32 cmp_aux14( Int32 b1, Int32 b2)
{
  Int32 b1_start;
  UChar c1, c2;

   b1_start=b1;
   b1 += 2; b2 += 2;

#if (N_OVERSHOOT<16)
#error
#endif

   // execute blocks of 14 comparisons untill a difference
   // is found or we run out of the string 
   do {
     // 3
     c1 = block[b1]; c2 = block[b2];
     if (c1 != c2) {Curr_len=b1-b1_start; return ((UInt32)c1 - (UInt32)c2);}
     b1++; b2++; 

     // 4
     c1 = block[b1]; c2 = block[b2];
     if (c1 != c2) {Curr_len=b1-b1_start; return ((UInt32)c1 - (UInt32)c2);}
     b1++; b2++; 

     // 5
     c1 = block[b1]; c2 = block[b2];
     if (c1 != c2) {Curr_len=b1-b1_start; return ((UInt32)c1 - (UInt32)c2);}
     b1++; b2++; 

     // 6
     c1 = block[b1]; c2 = block[b2];
     if (c1 != c2) {Curr_len=b1-b1_start; return ((UInt32)c1 - (UInt32)c2);}
     b1++; b2++; 

     // 7
     c1 = block[b1]; c2 = block[b2];
     if (c1 != c2) {Curr_len=b1-b1_start; return ((UInt32)c1 - (UInt32)c2);}
     b1++; b2++; 

     // 8
     c1 = block[b1]; c2 = block[b2];
     if (c1 != c2) {Curr_len=b1-b1_start; return ((UInt32)c1 - (UInt32)c2);}
     b1++; b2++; 
   
     // 9
     c1 = block[b1]; c2 = block[b2];
     if (c1 != c2) {Curr_len=b1-b1_start; return ((UInt32)c1 - (UInt32)c2);}
     b1++; b2++; 

     // 10
     c1 = block[b1]; c2 = block[b2];
     if (c1 != c2) {Curr_len=b1-b1_start; return ((UInt32)c1 - (UInt32)c2);}
     b1++; b2++; 

     // 11
     c1 = block[b1]; c2 = block[b2];
     if (c1 != c2) {Curr_len=b1-b1_start; return ((UInt32)c1 - (UInt32)c2);}
     b1++; b2++; 

     // 12
     c1 = block[b1]; c2 = block[b2];
     if (c1 != c2) {Curr_len=b1-b1_start; return ((UInt32)c1 - (UInt32)c2);}
     b1++; b2++; 

     // 13
     c1 = block[b1]; c2 = block[b2];
     if (c1 != c2) {Curr_len=b1-b1_start; return ((UInt32)c1 - (UInt32)c2);}
     b1++; b2++; 

     // 14
     c1 = block[b1]; c2 = block[b2];
     if (c1 != c2) {Curr_len=b1-b1_start; return ((UInt32)c1 - (UInt32)c2);}
     b1++; b2++; 

     // 15
     c1 = block[b1]; c2 = block[b2];
     if (c1 != c2) {Curr_len=b1-b1_start; return ((UInt32)c1 - (UInt32)c2);}
     b1++; b2++; 

     // 16
     c1 = block[b1]; c2 = block[b2];
     if (c1 != c2) {Curr_len=b1-b1_start; return ((UInt32)c1 - (UInt32)c2);}
     b1++; b2++; 

   } while( b1<nBlock && b2<nBlock);

   Curr_len=b1-b1_start;
   return b2-b1;   // we have  b2>b1 <=> s2<s1
} 


/* ------------------------------------------------------------
   insert the data for a long match in Lm_list.
   the Lm_list is sorted for increasing diff. elements of equal diff
   are sorted for increasing start
   ------------------------------------------------------------ */ 
void lm_insert(Int32 diff, Int32 start, Int32 len, Bool ris)
{
  Int32 i,h1,h2;

  if(Verbose>2) fprintf(stderr,"#");
  /* ----- try to backward extend a block ------ */ 
  for(i=start-1;i>=0;i--) 
    if(block[i]!=block[i+diff]) break;
  len += (start-i-1);
  start = i+1;

  if(++Num_lm==LM_LIST_SIZE) {
    fprintf(stderr,"Too many long matches!\n");
    exit(1);
  }
  // search for an empty position in the hash table 
  h1 = HASH1(diff); h2=HASH2(diff);
  while(Hash_table_diff[h1]!=0) {
    h1 += h2; h1 &= (HASH_SIZE-1);
  }   
  Hash_table_diff[h1]=diff;
  Hash_table_start[h1]=start;
  Hash_table_lenNris[h1] = (len<<1) | ((Int32) ris);
}

/* -------------------------------------------------------------------
   search for a long match with difference diff. This function return
   either:
     1) an integer h1>=0 with identify a long match lm such that:
        lm.diff==diff,  lm.start<=start and  lm.start + lm.len > start
     2) -1
   NOTE: this function is no loger used: it has been in-lined in mycmp  
   -------------------------------------------------------------------- */ 
__inline__
Int32 lm_search(Int32 diff, Int32 start)
{
  Int32 h1,h2,d;

  h1=HASH1(diff);
  if(Hash_table_diff[h1]==0) /* no matches with difference=diff*/
    return -1;           
  h2=HASH2(diff);            /* scan the hash table */
  while((d=Hash_table_diff[h1])!=0) {
    if(d==diff) {
      if(Hash_table_start[h1] <= start && 
	 Hash_table_start[h1]+(Hash_table_lenNris[h1]>>1)>start)
	return h1;
    }
    h1 += h2; h1 &= (HASH_SIZE-1);
  }
  return -1;
}


static __inline__
Int32 my_cmp(Int32 b1, Int32 b2)
{
  void lm_insert(Int32 diff, Int32 start, Int32 len, Bool ris);
  Int32 first, second, diff, ris, pmris, xch, h1, h2, d;

  /* ---------- try the easy way -------- */
  if(b1==b2) return 0;
  if(Num_lm==0) {
    pmris = cmp_aux14(b1,b2);
    if(Curr_len<=LM_MIN_LEN) return pmris;
  }

  /* ---- reorder b1 and b2 ------- */
  if(b1>b2) {
    first = b2; second=b1; xch = 1;
  }
  else {
    first = b1; second=b2; xch = 0;
  }
  diff = second-first;

  /* --------- check if there are relevant long matches -------- */
  /* ---- this used to be the function lm_search() ------------- */
  h1=HASH1(diff);
  if((d=Hash_table_diff[h1])!=0) {        /* there are matches with difference=diff*/
    h2=HASH2(diff);                       /* scan the hash table */
    do {
      if(d==diff) {
	if(Hash_table_start[h1] <= first && 
	   Hash_table_start[h1]+(Hash_table_lenNris[h1]>>1)>first)
	  break;
      }
      h1 += h2; h1 &= (HASH_SIZE-1);    // update h1
      d=Hash_table_diff[h1];            // get new diff
      if(d==0) 
	goto make_comparison;
    } while(1);
    ris = (Hash_table_lenNris[h1] & 1);   /* get result of comparison */
    pmris =  2*ris -1;                    /* translate 0/1 ris to +=1 ris */
  }
  else {                            /* nothing useful in the hash table */ 
  make_comparison:
    pmris = cmp_aux14(first,second);
    assert(pmris!=0);
    if(Curr_len>LM_MIN_LEN) { 
      pmris = pmris>0 ? 1 : -1;
      lm_insert(diff,first,Curr_len,(pmris+1)/2);
    }
  }
  if(xch) pmris = -pmris;
  return pmris;
}




#if !MK
#define  TEMPLATE  d_copyEQ_u12
#define  CMP       my_cmp
#define  FMAP      ptr
#include "qsort3.template"
#endif


void calc_running_order ( void )
{
   Int32 i, j;
   for (i = 0; i <= 255; i++) runningOrder[i] = i;

   {
      Int32 vv;
      Int32 h = 1;
      do h = 3 * h + 1; while (h <= 256);
      do {
         h = h / 3;
         for (i = h; i <= 255; i++) {
            vv = runningOrder[i];
            j = i;
            while ( BIGFREQ(runningOrder[j-h]) > BIGFREQ(vv) ) {
               runningOrder[j] = runningOrder[j-h];
               j = j - h;
               if (j <= (h - 1)) goto zero;
            }
            zero:
            runningOrder[j] = vv;
         }
      } while (h != 1);
   }
}



void d_copyEQ_u12 ( void )
{
   void ssort1main(UChar *, int *, int);
   Int32  i, j, k, ss, sb;
   UChar  c1, c2;
   Bool   bigDone[256];
   Int32  copyStart[256];
   Int32  copyEnd  [256];
   Int32  numQSorted = 0;

   for (i = 0; i < N_OVERSHOOT; i++) // init tail
      block[nBlock + i] = 0;

   for (i = 0; i <= 65536; i++) ftab[i] = 0;

   c1 = block[0];
   for (i = 1; i <= nBlock; i++) {
      c2 = block[i];
      ftab[(c1 << 8) + c2]++;
      c1 = c2;
   } 

   for (i = 1; i <= 65536; i++) ftab[i] += ftab[i-1];

   c1 = block[0];
   for (i = 0; i < nBlock; i++) {
      c2 = block[i+1];
      j = (c1 << 8) + c2;
      c1 = c2;
      ftab[j]--;
      ptr[ftab[j]] = i;
   }

   /* decide on the running order */
   calc_running_order();
   for (i = 0; i < 256; i++) bigDone[i] = False;

   /* Really do the sorting */
   for (i = 0; i <= 255; i++) {
  
      /*--
         Process big buckets, starting with the least full.
      --*/
      ss = runningOrder[i];

      /*--
         Complete the big bucket [ss] by quicksorting
         any unsorted small buckets [ss, j].  Hopefully
         previous pointer-scanning phases have already
         completed many of the small buckets [ss, j], so
         we don't have to sort them at all.
      --*/
      for (j = 0; j <= 255; j++) {
         if (j != ss) {
            sb = (ss << 8) + j;
            if ( ! (ftab[sb] & SETMASK) ) {
               Int32 lo = ftab[sb]   & CLEARMASK;
               Int32 hi = (ftab[sb+1] & CLEARMASK) - 1;
               if (hi > lo) {
                  if (local_debug)
                     fprintf ( stderr,
                               "        qsort [0x%x, 0x%x]   done %d"
                               "   this %d\n",
                               ss, j, numQSorted, hi - lo + 1 );
                  #if MK
                  ssort1main(block, ptr+lo, hi-lo+1);
                  #else
                  d_copyEQ_u12_qsort3 ( lo, hi );
                  #endif
                  numQSorted += ( hi - lo + 1 );
               }
            }
            ftab[sb] |= SETMASK;
         }
      }

      assert (!bigDone[ss]);

      {
         for (j = 0; j <= 255; j++) {
            copyStart[j] =  ftab[(j << 8) + ss]     & CLEARMASK;
            copyEnd  [j] = (ftab[(j << 8) + ss + 1] & CLEARMASK) - 1;
         }
	 // take care of the virtual -1 char in position nBlock+1
	 if(ss==0) {
	   k=nBlock-1;
	   c1 = block[k];
	   if (!bigDone[c1])
	     ptr[ copyStart[c1]++ ] = k;
	 }
         for (j = ftab[ss << 8] & CLEARMASK; j < copyStart[ss]; j++) {
            k = ptr[j]-1; if (k < 0) continue;   // k += nBlock;
            c1 = block[k];
            if (!bigDone[c1])
               ptr[ copyStart[c1]++ ] = k;
         }
         for (j = (ftab[(ss+1) << 8] & CLEARMASK) - 1; j > copyEnd[ss]; j--) {
            k = ptr[j]-1; if (k < 0) continue; // k += nBlock;
            c1 = block[k];
            if (!bigDone[c1]) 
               ptr[ copyEnd[c1]-- ] = k;
         }
      }

      assert (copyStart[ss] - 1 == copyEnd[ss]);

      for (j = 0; j <= 255; j++) ftab[(j << 8) + ss] |= SETMASK;
      bigDone[ss] = True;

   }

   if (local_debug)
   fprintf ( stderr, "        %d pointers, %d sorted, %d scanned\n",
                     nBlock, numQSorted, nBlock - numQSorted );

}


#if 0
// -----------------------------------------------------------------
// this is a macro equivalent to cmp_aux14()
// it introduces a saving of roughly 5% in the running time
// -----------------------------------------------------------------
#define cmp_aux14_macro(b1,b2,ris)          \
   b1_start=b1;                     \
   b1 += 2; b2 += 2;                \
   do {                             \
     c1 = block[b1]; c2 = block[b2];\
     if (c1 != c2) {Curr_len=b1-b1_start; ris=((UInt32)c1 - (UInt32)c2);goto done_cmp;}\
     b1++; b2++;                    \
     c1 = block[b1]; c2 = block[b2];\
     if (c1 != c2) {Curr_len=b1-b1_start; ris=((UInt32)c1 - (UInt32)c2);goto done_cmp;}\
     b1++; b2++;                    \
     c1 = block[b1]; c2 = block[b2];\
     if (c1 != c2) {Curr_len=b1-b1_start; ris=((UInt32)c1 - (UInt32)c2);goto done_cmp;}\
     b1++; b2++;                    \
     c1 = block[b1]; c2 = block[b2];\
     if (c1 != c2) {Curr_len=b1-b1_start; ris=((UInt32)c1 - (UInt32)c2);goto done_cmp;}\
     b1++; b2++;                    \
     c1 = block[b1]; c2 = block[b2];\
     if (c1 != c2) {Curr_len=b1-b1_start; ris=((UInt32)c1 - (UInt32)c2);goto done_cmp;}\
     b1++; b2++;                    \
     c1 = block[b1]; c2 = block[b2];\
     if (c1 != c2) {Curr_len=b1-b1_start; ris=((UInt32)c1 - (UInt32)c2);goto done_cmp;}\
     b1++; b2++;                    \
     c1 = block[b1]; c2 = block[b2];\
     if (c1 != c2) {Curr_len=b1-b1_start; ris=((UInt32)c1 - (UInt32)c2);goto done_cmp;}\
     b1++; b2++;                    \
     c1 = block[b1]; c2 = block[b2];\
     if (c1 != c2) {Curr_len=b1-b1_start; ris=((UInt32)c1 - (UInt32)c2);goto done_cmp;}\
     b1++; b2++;                    \
     c1 = block[b1]; c2 = block[b2];\
     if (c1 != c2) {Curr_len=b1-b1_start; ris=((UInt32)c1 - (UInt32)c2);goto done_cmp;}\
     b1++; b2++;                    \
     c1 = block[b1]; c2 = block[b2];\
     if (c1 != c2) {Curr_len=b1-b1_start; ris=((UInt32)c1 - (UInt32)c2);goto done_cmp;}\
     b1++; b2++;                    \
     c1 = block[b1]; c2 = block[b2];\
     if (c1 != c2) {Curr_len=b1-b1_start; ris=((UInt32)c1 - (UInt32)c2);goto done_cmp;}\
     b1++; b2++;                    \
     c1 = block[b1]; c2 = block[b2];\
     if (c1 != c2) {Curr_len=b1-b1_start; ris=((UInt32)c1 - (UInt32)c2);goto done_cmp;}\
     b1++; b2++;                    \
     c1 = block[b1]; c2 = block[b2];\
     if (c1 != c2) {Curr_len=b1-b1_start; ris=((UInt32)c1 - (UInt32)c2);goto done_cmp;}\
     b1++; b2++;                    \
     c1 = block[b1]; c2 = block[b2];\
     if (c1 != c2) {Curr_len=b1-b1_start; ris=((UInt32)c1 - (UInt32)c2);goto done_cmp;}\
     b1++; b2++;                    \
   } while( b1<nBlock && b2<nBlock);\
   Curr_len=b1-b1_start;            \
   ris=b2-b1;   /* we have  b2>b1 <=> s2<s1 */ \
 done_cmp:      
#endif




