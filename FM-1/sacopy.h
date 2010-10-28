// constant used by the seward-like sa construction algorithms

#define N_OVERSHOOT (16)      // length of overshoot for loop unrolling 
#define LM_MIN_LEN 1000       // minimum length of long matches
#define LM_LIST_SIZE 4000     // max # of long matches stores
#define HASH_SIZE 8192        // size of hash table for long matches

// ---- make sure the hash table stays half empty -----
#if((HASH_SIZE/2) < LM_LIST_SIZE)
#error
#endif

// --- hashing functions 
#define FIINV 5063 // this is 8192*0.618033988
#define HASH1(k)  ( ((k)*FIINV) & (HASH_SIZE-1) ) 
#define HASH2(k) (( ((k>>7)*FIINV) & (HASH_SIZE-1) ) | 1) 

