
// Implements the revtrie data structure

#include "revtrie.h"

 // The idea is to use the revtrie as follows. Since we have to search for
 // every substring of P in lztrie, let's assume it has been already done, 
 // which is easier to handle. Note that every substring s searched in lztrie 
 // corresponds to the same block number for s^R in revtrie, so those cells 
 // (i,j) of C that have an answer in lztrie have their answer here. Indeed,
 // we need to search only for prefixes in revtrie, so to search for (0,j) we
 // start from the smallest i such that (i,j) has an answer in lztrie.
 // We find (i,j) and have to find out our way to (0,i-1). We want to know
 // which child goes down by a letter c, if any. For each child, we go to 
 // lztrie, go up the part that corresponds to P^R (i..j), and see which is 
 // the parent. After we get the correct child, we continue from there until
 // we use up all the string that connects the current node with the correct
 // child. This is known by going up the child and also the largest descendant
 // of it: when they differ by a letter, the common string has been read. At
 // that point we finally go down to the child and continue until finding, if
 // at all, (0,j).

        // creates a revtrie structure from a parentheses bitstring,
        // a corresponding lztrie, the mapping from block id to lztrie,
        // and an id array in preorder, which are not owned
        // n is the total number of nodes (n ids, 2n parentheses)

revtrie createRevTrie (uint *string, lztrie trie, nodemap map, uint *id, 
                       uint n)

   { revtrie T;
     uint i;
     T = malloc (sizeof(struct srevtrie));
     T->data = string;
     T->pdata = createParentheses (string,2*n,false);
     T->n = n;
     T->nbits = bits(n-1);
     T->trie = trie;
     T->map = map;
     T->id = malloc (((n*T->nbits+W-1)/W)*sizeof(uint));
     for (i=0;i<n;i++)
	bitput (T->id,i*T->nbits,T->nbits,id[i]);
     return T;
   }

	// frees revtrie structure, including the owned data

void destroyRevTrie (revtrie T)
   { destroyParentheses (T->pdata);
     free (T->id);
     free (T);
   }

	// stores the revtrie on file f

void saveRevTrie (revtrie T, FILE *f)

   { if (fwrite (&T->n,sizeof(uint),1,f) != 1)
        { fprintf (stderr,"Error: Cannot write RevTrie on file\n");
          exit(1);
        }
     if (fwrite (T->data,sizeof(uint),(2*T->n+W-1)/W,f) != (2*T->n+W-1)/W)
        { fprintf (stderr,"Error: Cannot write RevTrie on file\n");
          exit(1);
        }                                                                       
     if (fwrite (T->id,sizeof(uint),(T->n*T->nbits+W-1)/W,f) !=
                 (T->n*T->nbits+W-1)/W)
        { fprintf (stderr,"Error: Cannot write RevTrie on file\n");
          exit(1);
        }
   }                                                                            

	// loads revtrie from file f

revtrie loadRevTrie (FILE *f, lztrie trie, nodemap map)

   { revtrie T;
     T = malloc (sizeof(struct srevtrie));
     if (fread (&T->n,sizeof(uint),1,f) != 1)
        { fprintf (stderr,"Error: Cannot read RevTrie from file\n");
          exit(1);
        }
     T->data = malloc (((2*T->n+W-1)/W)*sizeof(uint));
     if (fread (T->data,sizeof(uint),(2*T->n+W-1)/W,f) != (2*T->n+W-1)/W)
        { fprintf (stderr,"Error: Cannot read RevTrie from file\n");
          exit(1);
        }                                                                       
     T->pdata = createParentheses (T->data,2*T->n,false);
     T->nbits = bits(T->n-1);
     T->trie = trie;
     T->map = map;
     T->id = malloc (((T->n*T->nbits+W-1)/W)*sizeof(uint));
     if (fread (T->id,sizeof(uint),(T->n*T->nbits+W-1)/W,f) !=
                 (T->n*T->nbits+W-1)/W)
        { fprintf (stderr,"Error: Cannot read RevTrie from file\n");
          exit(1);
        }
     return T;
   }

	// give children by (a string starting in) letter c, if possible.
	// leave in lzl,lzr the corresponding lztrie nodes (leftmost and
	// rightmost) below the node returned, if it is empty, otherwise
	// lzl=lztrie node corresponding to the node returned,lzr=null.
	// depth is the length of the string represented by node i

#ifdef QUERYREPORT
int BCALLS = 0; // calls to childRevTrie
int DEPTH = 0;  // parent operations due to depth in the trie + search
int RJUMPS = 0; // number of nodes considered to go by child
int DEPTH1 = 0; // like depth but once per call
#endif

trieNode childRevTrie (revtrie T, trieNode i, uint depth, byte c, 
			trieNode *lzl, trieNode *lzr)

   { trieNode j;
     trieNode jlz;
     uint d,pos,jid;
     byte tc;
#ifdef QUERYREPORT
     BCALLS++;
#endif
     j = i+1;
#ifdef QUERYREPORT
		DEPTH1 += depth;
#endif
     while (!bitget1(T->data,j)) // there is a child here
	{ pos = leftrankRevTrie(T,j); 
	  jid = rthRevTrie(T,pos); // leftmost id in subtree
	  jlz = mapto(T->map,jid);
	  for (d=depth;d;d--) jlz = parentLZTrie(T->trie,jlz);
#ifdef QUERYREPORT
		RJUMPS++;
		DEPTH += depth;
#endif
	  tc = letterLZTrie(T->trie,jlz);
	  if (tc > c) break;
	  if (tc == c)
	     { *lzl = jlz; 
	       if ((pos < T->n-1) && (rthRevTrie(T,pos+1) == jid)) // empty
	          { pos = rightrankRevTrie(T,j); 
	            jlz = mapto(T->map,rthRevTrie(T,pos)); // rightmost
	            for (d=depth;d;d--) jlz = parentLZTrie(T->trie,jlz);
	            *lzr = jlz;
		  }
	       else *lzr = NULLT;
	       return j; 
	     }
	  j = findclose(T->pdata,j)+1;
	}
     *lzl = *lzr = NULLT;
     return NULLT; // no child to go down by c
   }

	// subtree size

uint subtreesizeRevTrie (revtrie T, trieNode i)

   { trieNode j;
     j = findclose (T->pdata,i);
     return (j-i+1)/2;
   }

	// smallest rank in subtree

uint leftrankRevTrie (revtrie T, trieNode i)

   { return i-rank(T->pdata->bdata,i);
   }

	// largest rank in subtree

uint rightrankRevTrie (revtrie T, trieNode i)

   { trieNode j;
     j = findclose (T->pdata,i);
     return j-rank(T->pdata->bdata,j)-1;
   }

	// id of node

uint idRevTrie (revtrie T, trieNode i)

   { uint pos = leftrankRevTrie (T,i);
     return rthRevTrie (T,pos);
   }

	// rth of position

uint rthRevTrie (revtrie T, uint pos)

   { return bitget (T->id,pos*T->nbits,T->nbits);
   }

        // is node i ancestor of node j?

bool ancestorRevTrie (revtrie T, trieNode i, trieNode j)

   { return (i <= j) && (j < findclose (T->pdata,i));
   } 

        // next node from i, in preorder
        // assumes it *can* go on!

trieNode nextRevTrie (revtrie T, trieNode i)

   { uint *data = T->data; 
     i++; 
     while (bitget1(data,i)) i++;
     return i;
   }



uint sizeofRevTrie(revtrie T, uint n)
 {
    return sizeof(struct srevtrie) +
           sizeofParentheses(T->pdata) +
           ((n/W)*T->nbits+1)*sizeof(uint); // rids     
 }
