
#include <math.h>

 
position createPosition(lzindex I, uint text_length)
 {
    position P;
    uint *Offset;
    ulong n, M, superblock_size, current_superblock,
          starting_pos, i, depth, pos;
    trieNode node;
    lztrie T = I.fwdtrie;
    revtrie RT = I.bwdtrie; 
    nodemap rmap = I.rmap;
    
    
    P = malloc(sizeof(struct tpos));
    n = T->n;
    P->SBlock_size = 32;//bits(n-1); // superblock size
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
       node  = NODE(I, i);
       depth = depthLZTrie(T, node);
       pos  += depth;
       
       superblock_size++;
       
       if ((superblock_size > P->SBlock_size) || (i==0)) {
          if (i) 
             bitput(P->SuperBlock,P->nbitsSB*(++current_superblock),P->nbitsSB,pos);
          if (i!=0 && Offset[i-1]>M) M = Offset[i-1];
          superblock_size = 1;
          starting_pos    = pos;
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
       ulong posSB = (id-1)>>5;
       return bitget(P->SuperBlock,posSB*P->nbitsSB,P->nbitsSB)
            + bitget(P->Offset,(id-1)*P->nbitsOffs,P->nbitsOffs);
    }
 }

 
// given a text position, gets the identifier of the
// LZ78 phrase containing that position
 
ulong getphrasePosition(position P, ulong text_pos)
 {
    ulong li, ls, nbits, elem, med, SB, i;
    
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
    
    
    li = med*P->SBlock_size;
    if (med+1 == P->nSuperBlock) ls = P->nOffset-1;
    else ls = (med+1)*P->SBlock_size - 1;
    nbits = P->nbitsOffs;
    
    while ((ls-li+1) > 6) {
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
 

position loadPosition(FILE *f, uint text_length)
 {
    position P;
    uint aux;

    P = malloc (sizeof(struct tpos));

    if (fread(&P->nSuperBlock,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot read LZTrie from file\n");
       exit(1);
    }

    P->nbitsSB = bits(text_length-1);
    aux = (((unsigned long long) P->nbitsSB*P->nSuperBlock+W-1)/W);
    P->SuperBlock  = malloc(aux*sizeof(uint));

    if (fread(P->SuperBlock,sizeof(uint),aux,f) != aux) {
       fprintf(stderr,"Error: Cannot read LZTrie from file\n");
       exit(1);
    }

    P->SBlock_size = 32;

    P->Tlength = text_length;

    if (fread(&P->nOffset,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot read LZTrie from file\n");
       exit(1);
    }

    if (fread(&P->nbitsOffs,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot read LZTrie from file\n");
       exit(1);
    }

    aux = (((unsigned long long)P->nOffset*P->nbitsOffs+W-1)/W);
    P->Offset    = malloc(aux*sizeof(uint));
    if (fread(P->Offset,sizeof(uint),aux,f) != aux) {
       fprintf(stderr,"Error: Cannot read LZTrie from file\n");
       exit(1);
    }

    return P;
 }


void savePosition(FILE *f, position P)
 {
    uint aux;

    if (fwrite(&P->nSuperBlock,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot write LZTrie on file\n");
       exit(1);
    }

    aux = (((unsigned long long) P->nbitsSB*P->nSuperBlock+W-1)/W);

    if (fwrite(P->SuperBlock,sizeof(uint),aux,f) != aux) {
       fprintf(stderr,"Error: Cannot write LZTrie on file\n");
       exit(1);
    }

    if (fwrite(&P->nOffset,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot write LZTrie on file\n");
       exit(1);
    }

    if (fwrite(&P->nbitsOffs,sizeof(uint),1,f) != 1) {
       fprintf(stderr,"Error: Cannot write LZTrie on file\n");
       exit(1);
    }

    aux = (((unsigned long long)P->nOffset*P->nbitsOffs+W-1)/W);

    if (fwrite(P->Offset,sizeof(uint),aux,f) != aux) {
       fprintf(stderr,"Error: Cannot write LZTrie on file\n");
       exit(1);
    }

 }

 
 
 
 
 
   
ulong sizeofPosition(position P)
 {
    return sizeof(struct tpos)
          + ((P->nbitsSB*P->nSuperBlock+W-1)/W)*sizeof(uint)
          + ((P->nOffset*P->nbitsOffs+W-1)/W)*sizeof(uint); 
 }
 
void destroyPosition(position P)
 {
    free(P->SuperBlock);
    free(P->Offset);
 }

