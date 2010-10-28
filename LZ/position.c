
// This code implements the data structure to get
// "real" text positions at search time. 

#include "position.h"
#include <math.h>

 
position createPosition(lztrie T, uint text_length, nodemap map)
 {
    position P;
    uint *Offset;
    ulong n, M, superblock_size, current_superblock,
          starting_pos, i, depth, pos;
    trieNode node;
    
    P = malloc(sizeof(struct tpos));
    n = T->n;
    P->SBlock_size = bits(n-1); // superblock size is log n
    P->nbitsSB     = bits(text_length-1);
    P->nSuperBlock = (ulong)ceil((double)n/P->SBlock_size); // number of superblocks
    P->SuperBlock  = malloc(((P->nbitsSB*P->nSuperBlock+W-1)/W)*sizeof(uint));
    P->Tlength     = text_length;
    M = 0; // maximum superblock size (in number of characters)
    
    current_superblock = 0;
    superblock_size    = 0;
    starting_pos       = 0;
    P->nOffset         = n;
    Offset = malloc(n*sizeof(uint)); // array for temporary use
    
    bitput(P->SuperBlock, 0, P->nbitsSB, 0); 
    
    for (i = 0, pos = 0; i < n; i++) {
       node = mapto(map, i);
       depth = depthLZTrie(T, node);
       pos += depth;
       
       superblock_size++;
       
       if ((superblock_size > P->SBlock_size) || (i==0)) {
          if (i) 
             bitput(P->SuperBlock,P->nbitsSB*(++current_superblock),P->nbitsSB,pos);
          if (i!=0 && Offset[i-1]>M) M = Offset[i-1];
          superblock_size = 1;
          starting_pos = pos;
       }
       
       Offset[i] = pos-starting_pos; 
    }
    
    if (Offset[i-1]>M) M = Offset[i-1];
        
    P->nbitsOffs = bits(M-1);
    P->Offset    = malloc(((P->nOffset*P->nbitsOffs+W-1)/W)*sizeof(uint));
    for (i = 0; i < n; i++)
       bitput(P->Offset, i*P->nbitsOffs, P->nbitsOffs, Offset[i]);
    
    free(Offset);
    
    return P;
 }
 
 
// given a phrase id, gets the corresponding text position

ulong getPosition(position P, uint id)
 {
     if (id > P->nOffset) return P->Tlength; 
     else {
       ulong posSB = (ulong)ceil((double)id/P->SBlock_size)-1;
       return bitget(P->SuperBlock,posSB*P->nbitsSB,P->nbitsSB)
            + bitget(P->Offset,(id-1)*P->nbitsOffs,P->nbitsOffs);
     }
 }

 
// given a text position, gets the identifier of the
// LZ78 phrase containing that position
 
ulong getphrasePosition(position P, ulong text_pos)
 {
    ulong li, ls, nbits, elem, med, SB, temp, i;
    
    li = 0;
    ls = P->nSuperBlock-1;
    nbits = P->nbitsSB;
    while ((ls-li+1) > 0) {
       med  = (li+ls)/2;
       elem = bitget(P->SuperBlock,med*nbits,nbits);
       if (elem == text_pos) break;
       if (elem < text_pos) li = med+1;
       else ls = med - 1;
    }
    if (elem > text_pos && med) med--;
    SB = bitget(P->SuperBlock,med*nbits,nbits);
    
    temp = li;
    li = med*P->SBlock_size;
    
    if (med+1 == P->nSuperBlock) 
       ls = P->nOffset-1;
    else 
       ls = (med+1)*P->SBlock_size - 1;
    
    nbits = P->nbitsOffs;

    while ((ls-li+1) > 8) {
       med  = (li+ls)/2;
       elem = bitget(P->Offset,med*nbits,nbits)+SB;
       if (elem == text_pos) break;
       if (elem < text_pos) li = med+1;
       else ls = med-1;
    }
            
    for (i = li; i <= ls; i++)
      if ((elem=bitget(P->Offset,i*nbits,nbits) + SB) >= text_pos) break;
    
    if (elem > text_pos) i--;
    
    return i+(elem>text_pos);
 }
 
  
ulong sizeofPosition(position P)
 {
    return sizeof(struct tpos)
          + (P->nbitsSB*(P->nSuperBlock/W)+1)*sizeof(uint)
          + ((P->nOffset/W)*P->nbitsOffs+1)*sizeof(uint); 
 }
 
void destroyPosition(position P)
 {
    free(P->SuperBlock);
    free(P->Offset);
 }
