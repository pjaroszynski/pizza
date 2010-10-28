
 // LZ78 trie data structure

#include "trie.h"

	// creates trie

trie createTrie (void)

   { trie pTrie;
     uint i;
     pTrie = malloc (sizeof(struct strie));
     pTrie->nid = 0;
     pTrie->trie.id = pTrie->nid++;
     pTrie->trie.nchildren = 0;
     pTrie->trie.children = NULL;
     pTrie->heaps[0] = createHeap(sizeof(triebody));
     for (i=1;i<256;i++)
	pTrie->heaps[i] = createHeap(i*sizeof(struct schild));
     return pTrie;
   }

	// frees the trie

void destroyTrie (trie pTrie)

   { uint i;
     for (i=0;i<256;i++) destroyHeap (pTrie->heaps[i]);
     free (pTrie);
   }

	// inserts word[0...] into pTrie and returns new text ptr
	// insertion proceeds until we get a new trie node

byte *insertTrie (trie pTrie, byte *word)

   { triebody *t = &pTrie->trie;
     triebody *nt;
     struct schild *newc;
     int i,j;
     int m = 0;
	// traverse pTrie with word[0...]
     while (true)
       { i = 0;
         while (i < t->nchildren)
	    { if (t->children[i].car >= word[m]) break;
	      i++;
	    }
	 if ((i == t->nchildren) || (t->children[i].car > word[m]))
	    break;  // not found, get out
         t = t->children[i].trie;
	 m++;
       }
	// at this point we fell off the trie, which is guaranteed to occur
	// since the text finishes with the unique character 0
     newc = mallocHeap(pTrie->heaps[t->nchildren+1]);
     memcpy (newc,t->children,i*sizeof(struct schild));
     memcpy (newc+i+1,t->children+i,(t->nchildren-i)*sizeof(struct schild));
     freeHeap (pTrie->heaps[t->nchildren],t->children);
     t->children = newc;
     t->children[i].car = word[m];
     nt = mallocHeap (pTrie->heaps[0]);
     t->children[i].trie = nt;
     t->nchildren++;
	// new node created	
     nt->id = pTrie->nid++;
     nt->nchildren = 0;
     nt->children = NULL;
     	// return rest of text
     return word+m+1;
   }

	// inserts word[0..len-1] into pTrie, with id = id
	// assumes that no two equal strings are ever inserted

void insertstringTrie (trie pTrie, byte *word, uint len, uint id)

   { triebody *t,*nt;
     uint i,j,m;
     struct schild *newc;
	// traverse pTrie with word[0...]
     t = &pTrie->trie;
     m = 0;
     while (m < len)
       { i = 0;
         while (i < t->nchildren)
	    { if (t->children[i].car >= word[m]) break;
	      i++;
	    }
	 if ((i == t->nchildren) || (t->children[i].car > word[m]))
	    break;  // not found, get out
         t = t->children[i].trie;
	 m++;
       }
	// if we fell off the trie, we create more (unary and empty) nodes
     while (m < len)
       { newc = mallocHeap(pTrie->heaps[t->nchildren+1]);
         memcpy (newc,t->children,i*sizeof(struct schild));
         memcpy (newc+i+1,t->children+i,(t->nchildren-i)*sizeof(struct schild));
         freeHeap (pTrie->heaps[t->nchildren],t->children);
         t->children = newc;
	 if ((t->id == ~0) && (t->nchildren == 1)) pTrie->nid++; //not mute now
         t->children[i].car = word[m];
         nt = mallocHeap (pTrie->heaps[0]);
	 nt->id = ~0; // empty node, at least for now
         nt->nchildren = 0;
         nt->children = NULL;
         t->children[i].trie = nt;
         t->nchildren++;
	 t = nt;
         m++; i = 0;
       }
	// new node created or existing node with id added
     t->id = id;
     if (t->nchildren <= 1) pTrie->nid++; //not mute now
   }

        // represents pTrie with parentheses, letters and ids

	// also returns the leftmost id
static uint traverse (triebody *t, uint *parent, byte *letters, uint *ids,
		      uint *pi, uint *pli)

   { uint i,chid,oldpli,myid;
     myid = t->id;
     if ((myid == ~0) && (t->nchildren == 1)) // mute node
	return traverse (t->children[0].trie,parent,letters,ids,pi,pli);
	// open parenthesis
     bitclean (parent,*pi); 
     (*pi)++;
     oldpli = *pli;
	// traverse children
     for (i=0;i<t->nchildren;i++)
        { (*pli)++;
	  if (letters) letters[*pli] = t->children[i].car;
	  chid = traverse (t->children[i].trie,parent,letters,ids,pi,pli);
	  if (myid == ~0) myid = chid; // if I hadn't, get leftmost as mine
        }
	// write id 
     ids[oldpli] = myid;
	// close parenthesis
     bitset (parent,*pi); 
     (*pi)++; 
	// return leftmost id
     return myid;
   }

void representTrie (trie pTrie, uint *parent, byte *letters, uint *ids)


   { uint pi,pli;
     if (letters) letters[0] = 0; // dummy value
     pi = 0; pli = 0;
     traverse (&pTrie->trie,parent,letters,ids,&pi,&pli);
   }

