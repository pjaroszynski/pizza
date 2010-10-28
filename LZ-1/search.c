 
// Search module

#ifdef QUERYREPORT

#include <sys/time.h>

struct tms time1,time2;
double fwdtime;
double bwdtime;
double type1time;
double type2time;
double type3time;
double ticks;

#endif


#define INIT_OCCARRAY 5000
#define ALPHA 0.7

static uint Count, SIZE_OCCARRAY=INIT_OCCARRAY;
static ulong *OccArray=NULL, *OffsetArray=NULL;
position TPos;


uint nbits_aux;


inline uint bitget_aux(uint *e, uint p)
 {
    register uint i=(p>>5), j=p&0x1F, answ;

    if (j+nbits_aux <= W)
       answ = (e[i] << (W-j-nbits_aux)) >> (W-nbits_aux);
    else
       answ = (e[i] >> j) | ((e[i+1]<<(W-j-nbits_aux)) >> (W-nbits_aux));
    return answ;
 }




uint p_aux3, nbits_aux3, *e_aux3;

inline uint bitget_aux3()
 {
    register uint i=(p_aux3>>5), j=p_aux3&0x1F, answ;

    if (j+nbits_aux3 <= W)
        answ = (e_aux3[i] << (W-j-nbits_aux3)) >> (W-nbits_aux3);
    else
       answ = (e_aux3[i] >> j) | ((e_aux3[i+1]<<(W-j-nbits_aux3))>>(W-nbits_aux3));
    return answ;
 }


uint *e_aux4, p_aux4, nbits_aux4;

inline uint bitget_aux4()
 {
    register uint i=(p_aux4>>5), j=p_aux4&0x1F, answ;

    if (j+nbits_aux4 <= W)
       answ = (e_aux4[i] << (W-j-nbits_aux4)) >> (W-nbits_aux4);
    else
       answ = (e_aux4[i] >> j) | ((e_aux4[i+1]<<(W-j-nbits_aux4))>> (W-nbits_aux4));
    return answ;
 }



inline static void report(uint id,uint delta)
 { 
    OccArray[Count]      = id;
    OffsetArray[Count++] = delta;
    if (Count == SIZE_OCCARRAY) {
       SIZE_OCCARRAY = (SIZE_OCCARRAY+1)/ALPHA;
       OccArray      = realloc(OccArray, SIZE_OCCARRAY*sizeof(ulong));
       OffsetArray   = realloc(OffsetArray, SIZE_OCCARRAY*sizeof(ulong));
    }    
 }

        // returns answ[i,j] = lztrie node corresponding to pat[i..j], or
        //    NULLT if it does not exist. answ is really used as answ[i*m+j]
        // time is O(m^2s) in the worst case, but probably closer to O(m^2)

#ifdef QUERYREPORT
int QUERIES = 0; // number of queries
int BCELLS = 0;  // number of bwd cells filled by hand (tot-ffill-empty)
int FFILL = 0;   // number of bwd cells filled using fwd info
#endif

static trieNode *fwdmatrix(lzindex I, byte *pat, uint m, uint **fwdid)
 { 
    trieNode *answ, cur, new;
    uint i,j,ptr, *id;
    lztrie T = I.fwdtrie;
#ifdef QUERYREPORT
    QUERIES++;
#endif
    answ = malloc(m*m*sizeof(trieNode));
    id   = malloc(m*m*sizeof(uint));
    ptr  = 0; // ptr = i*m+j
    for (i=0; i<m; i++) {
       cur = ROOT; ptr += i;
       for (j=i; j<m; j++) {
          new = childLZTrie(T,cur,pat[j]);
          if (new != NULLT) {
             answ[ptr] = new;
             id [ptr]  = idLZTrie(T,new);
             cur       = new;
             ptr++;
          }
          else // no children, so no more entries in answ
             for (;j<m;j++) {
                answ[ptr] = NULLT; 
                id[ptr]   = ~0; 
                ptr++;
             }
       }
    }
    *fwdid = id;
    return answ;
 }

        // returns answ[j] = revtrie node corresponding to pat[0..j]^R, or
        //    NULLT if it does not exist.
        //    note that it needs the lztrie LZT, the answers computed by 
        //    fwdmatrix (fwd) and the corresponding block ids (fedid)
        // time is O(m^3s) in the worst case, but probably closer to O(m^2)

static trieNode *bwdmatrix(lzindex I, byte *pat, uint m, uint *fwdid)
 { 
    int i,j,k;
    trieNode cur,cur2,lzcur,lzcur2,*answ;
    byte c1,c2;
    lztrie LZT = I.fwdtrie;
    revtrie T = I.bwdtrie;
    nodemap rmap = I.rmap;
    
    answ = malloc(m*sizeof(trieNode));
    // answ[j] = node(pat[0..j]);
    for (j=0; j<m; j++) { 
       for (i=0;i<=j;i++) if (fwdid[m*i+j] != ~0) break;
       if (i <= j) {// i is last nonempty node in the path
          cur = RNODE(I,fwdid[m*i+j]);
#ifdef QUERYREPORT
          FFILL += j-i+1;
#endif
       }
       else cur = ROOT;
       i--; // unresolved 0..i
       while (i >= 0) {// once per (empty) node traversed
          cur = childRevTrie(T,cur,j-i,pat[i],&lzcur,&lzcur2);
#ifdef QUERYREPORT
          BCELLS++;
#endif
          if (cur == NULLT) break; // no further prefixes exist
          while (--i >= 0) {// once per letter
             lzcur = parentLZTrie(LZT,lzcur); //cannot be NULLT
#ifdef QUERYREPORT
             BCELLS++;
#endif
             if (lzcur == ROOT) break; // arrived at cur
             c1 = letterLZTrie(LZT,lzcur);
             if (lzcur2 != NULLT) {
                lzcur2 = parentLZTrie(LZT,lzcur2);
                c2     = letterLZTrie(LZT,lzcur2);
                if (c1 != c2) break; // end of common path
             }
             if (c1 != pat[i]) // no further prefixes exist
             { cur = NULLT; i = -1; break; }
          }
       }
       answ[j] = cur;
    }
    return answ;
 }

        // reports occurrences of type 1
        // time is O(occ1),

static void reportType1(lzindex I, trieNode *fwd, trieNode *bwd, uint *fwdid, 
                        uint m, bool counting)
 { 
    uint from,to,id,oldid,siz,pos,delta;
    trieNode lzcur,new, oldlzcur;
    lztrie fwdtrie = I.fwdtrie;
    revtrie bwdtrie = I.bwdtrie;
    uint /**Rev,*/ nbits;
    
    if (bwd[m-1] == NULLT) return; // there is not any LZ78 phrase ending with P
    from  = leftrankRevTrie(bwdtrie,bwd[m-1]);
    to    = rightrankRevTrie(bwdtrie,bwd[m-1]);
    oldid = 0;
    //Rev = bwdtrie->Rev;
    nbits = fwdtrie->nbits;
    e_aux3 = fwdtrie->ids;
    nbits_aux3 = nbits;
    while (from <= to) {
       lzcur = mapto(bwdtrie->Rev, getposRevTrie(bwdtrie,from++));
       pos = leftrankLZTrie(fwdtrie,lzcur); 
       p_aux3 = pos*nbits;
       id  = bitget_aux3();
       if (oldid != id) {// replicas due to empty nodes in revtrie
          oldid = id;
          siz = subtreesizeLZTrie(fwdtrie,lzcur);
          if (counting) Count += siz; // just count, faster
          else { 
             pos   = leftrankLZTrie(fwdtrie,lzcur); 
             delta = m;
             while (true) {
                report(id,delta);
                if (!--siz) break;
                p_aux3 = (++pos)*nbits;
		id    = bitget_aux3();//rthLZTrie(fwdtrie,++pos);
                lzcur = nextLZTrie(fwdtrie,lzcur,&delta);
             }
          }
       }
    }
 }

        // reports occurrences of type 2
        // time is O(min range among the 2 one-dimensional ones)

#ifdef QUERYREPORT
int WORK2 = 0; // amount of checks at level 2
#endif 

static void reportType2(lzindex I, trieNode *fwd, trieNode *bwd, uint m, 
                        bool counting)
 { 
    uint i,ffrom,fto,bfrom,bto;
    uint id,*data,rvtpos;
    trieNode cur,close,fopen, bopen;
    lztrie fwdtrie = I.fwdtrie;
    revtrie bwdtrie = I.bwdtrie;
    nodemap rmap = I.rmap, Rev = I.Rev;
    
    data = bwdtrie->B->data;
    e_aux3 = fwdtrie->ids;
    nbits_aux3 = fwdtrie->nbits;
    for (i=1; i<m; i++) {
       fopen = fwd[i*m+m-1];
       bopen = bwd[i-1];
       if ((fopen != NULLT) && (bopen != NULLT)) {
          // compute corresponding LZTrie preorder interval [ffrom, fto]
          ffrom = leftrankLZTrie(fwdtrie,fopen);
          fto   = rightrankLZTrie(fwdtrie,fopen);
          
          // compute corresponding RevTrie preorder interval [bfrom, bto] 
          bfrom = leftrankRevTrie(bwdtrie,bopen);
          bto   = rightrankRevTrie(bwdtrie,bopen);
          
          // Traverse "smaller" interval (traversing RevTrie is more expensive)
          if (fto-ffrom <= 6*(bto-bfrom)) {
             // LZTrie preorder interval is smaller
             close = findclose(bwdtrie->pdata,bopen);
             while (ffrom <= fto) {
                p_aux3 = (ffrom++)*nbits_aux3;
		id     = bitget_aux3();
                cur    = mapto(rmap,id-1);// RNODE(I,id-1);
                if ((bopen <= cur) && (cur <= close))
                   if (counting) Count++;
                   else report(id-1,i);
             }
          }
          else {
             // RevTrie preorder interval is smaller
             close = findclose(fwdtrie->pdata,fopen);
             rvtpos = getposRevTrie(bwdtrie, bfrom);
	     while (bfrom <= bto) {
                if (bitget1(data, bfrom)) {
                   p_aux3 = leftrankLZTrie(fwdtrie, mapto(Rev, rvtpos))*nbits_aux3;
		   id = bitget_aux3();
		   rvtpos++;
                   cur = mapto(Rev, getposRevTrie(bwdtrie,leftrankRevTrie(bwdtrie,mapto(rmap,id+1))));
                   if ((fopen <= cur) && (cur <= close))
                      if (counting) Count++;
                      else report(id,i);
		}
		bfrom++;
             }
          }
       }
    }
 }

        // reports occurrences of type 3
        // time is O(m^3s) worst case, but probably closer to O(m^2)

        // hashing

typedef struct 
 { 
    uint k;
    short i,j;
 } helem;

static helem *hcreate(uint m, uint *size)
 { 
    helem *table;
    int i;
    *size = 1 << bits(m*m); // so factor is at least 2.0
    table = malloc (*size*sizeof(helem));
    (*size)--;
    for (i=*size;i>=0;i--) table[i].k = ~0;
    return table;
 }

static void hinsert(helem *table, uint size, uint k, uint i, uint j)
 { 
    uint key = (k*i) & size;
    while (table[key].k != ~0) key = (key + PRIME2) & size;
    table[key].k = k; table[key].i = i; table[key].j = j;
 }

static int hsearch(helem *table, uint size, uint k, uint i)
 { 
    uint key = (k*i) & size;
    uint sk;
    while ((sk = table[key].k) != ~0) {
       if ((sk == k) && (table[key].i == i)) return table[key].j;
       key = (key + PRIME2) & size;
    }
    return -1;
 }

#ifdef QUERYREPORT
int EXTEND3 = 0; // how much we worked to extend cells at level 3
int CHECK3 = 0; // how many cells were checked for matches
#endif

static void reportType3(lzindex I, trieNode *fwd, uint *fwdid, trieNode *bwd,
                        byte *pat, uint m, bool counting)
 {
    uint i,j,k,f,t,nt,ok;
    trieNode node;
    helem *table;
    uint tsize;
    lztrie fwdtrie = I.fwdtrie;
    revtrie bwdtrie = I.bwdtrie;
    nodemap rmap = I.rmap;
    
    // find and store all blocks contained in pat
    table = hcreate(m,&tsize);
    for (i=1; i<m; i++)
       for (j=i; j<m; j++)
          if (fwdid[i*m+j] != ~0)
             hinsert(table,tsize,fwdid[i*m+j],i,j);
    // now find maximal segments
    for (i=0;i<m;i++)
       for (j=i; j<m; j++)
          if ((k = fwdid[i*m+j]) != ~0) {
             f = i; t = j; ok = k;
#ifdef QUERYREPORT
             EXTEND3++;
#endif
             while ((nt = hsearch(table,tsize,k+1,t+1)) != -1) {
                fwdid[(t+1)*m+nt] = ~0; t = nt; k++; 
#ifdef QUERYREPORT
                EXTEND3++;
#endif
             }
             // now we know that pat[f,t] is a maximal sequence
             if ((k-ok+1) + (f > 0) + (t < m-1) < 3) continue;
             // ok, it spans 3 blocks at least
#ifdef QUERYREPORT
             CHECK3++;
#endif
             if (t < m-1) {// see if block k+1 descends from pat[t+1..m-1]
                if (k+1 == fwdtrie->n) continue; // this was the last block
                node = fwd[(t+1)*m+m-1];
                if (node == NULLT) continue; // rest does not exist 
                if (!ancestorLZTrie(fwdtrie,node,NODE(I,k+1))) continue;
             }
             // ok, test on block k+1 passed
             if (f == 0)  // left test not needed
                if (counting) Count++;
                else report(ok,j-i+1);
             else {
                node = bwd[f-1];
                if (node == NULLT) continue; // rest does not exist 
                if (!ancestorRevTrie(bwdtrie,node,RNODE(I,ok-1))) 
                   continue;
                if (counting) Count++;
                else report(ok-1,f);
             }
          }
    free(table);
 }


int count(void *index, byte *pattern, ulong length, ulong* numocc)
 {
    lzindex I = *(lzindex *)index;
    trieNode *fwd,*bwd;
    uint *fwdid;
    
    Count       = 0;    
    
#ifdef QUERYREPORT
    ticks= (double)sysconf(_SC_CLK_TCK);
    times(&time1);  
#endif

    fwd = fwdmatrix(I,pattern,length,&fwdid);
    
#ifdef QUERYREPORT
    times(&time2);
    fwdtime += (time2.tms_utime - time1.tms_utime)/ticks;
    times(&time1);
#endif    

    bwd = bwdmatrix(I,pattern,length,fwdid);
    
#ifdef QUERYREPORT
    times(&time2);
    bwdtime += (time2.tms_utime - time1.tms_utime)/ticks;
    times(&time1);
#endif    

    reportType1(I,fwd,bwd,fwdid,length,true);
    
#ifdef QUERYREPORT
    times(&time2);
    type1time += (time2.tms_utime - time1.tms_utime)/ticks;
    times(&time1);
#endif    
    
    if (length > 1) {
       reportType2(I,fwd,bwd,length,true);
    
#ifdef QUERYREPORT
       times(&time2);
       type2time += (time2.tms_utime - time1.tms_utime)/ticks;
       times(&time1);
#endif               
    }
    
    if (length > 2) {
       reportType3(I,fwd,fwdid,bwd,pattern,length,true);
    
#ifdef QUERYREPORT
       times(&time2);
       type3time += (time2.tms_utime - time1.tms_utime)/ticks;
       times(&time1);
#endif        
    }
    
    free(fwd); 
    free(bwd); 
    free(fwdid);
    *numocc     = Count;
    return 0; // no errors yet
 }
 

int locate(void *index, byte *pattern, ulong length, 
           ulong **occ, ulong *numocc)
 {
    lzindex I = *(lzindex *)index;
    trieNode *fwd, *bwd;
    uint *fwdid, i;
    ulong *occaux;
    uint nbits_SB, nbits_Offs, id, nOffset, Tlength;;
    
#ifdef QUERYREPORT
    ticks= (double)sysconf(_SC_CLK_TCK);
#endif

    TPos = I.TPos;  
    nOffset       = TPos->nOffset;
    Tlength       = TPos->Tlength;
    Count         = 0;
    OccArray      = malloc(sizeof(ulong)*SIZE_OCCARRAY);
    if (OffsetArray == NULL)
       OffsetArray   = malloc(sizeof(ulong)*SIZE_OCCARRAY);
        
#ifdef QUERYREPORT
    times(&time1);
#endif    

    fwd = fwdmatrix(I,pattern,length,&fwdid);
    
#ifdef QUERYREPORT
    times(&time2);
    fwdtime += (time2.tms_utime - time1.tms_utime)/ticks;
    times(&time1);
#endif    

    bwd = bwdmatrix(I,pattern,length,fwdid);
    
#ifdef QUERYREPORT
    times(&time2);
    bwdtime += (time2.tms_utime - time1.tms_utime)/ticks;
    times(&time1);
#endif        

    reportType1(I,fwd,bwd,fwdid,length,false);
    
#ifdef QUERYREPORT
    times(&time2);
    type1time += (time2.tms_utime - time1.tms_utime)/ticks;
    times(&time1);
#endif        
    
    if (length > 1) {
       reportType2(I,fwd,bwd,length,false);

#ifdef QUERYREPORT
       times(&time2);
       type2time += (time2.tms_utime - time1.tms_utime)/ticks;
       times(&time1);
#endif    
    }
    

    if (length > 2) {
       reportType3(I,fwd,fwdid,bwd,pattern,length,false);

#ifdef QUERYREPORT
       times(&time2);
       type3time += (time2.tms_utime - time1.tms_utime)/ticks;
       times(&time1);
#endif            
    }

    free(fwd); 
    free(bwd); 
    free(fwdid);
    *numocc = Count;
    
    e_aux4        = TPos->SuperBlock;
    nbits_aux4    = TPos->nbitsSB;
    e_aux3        = TPos->Offset;
    nbits_aux3    = TPos->nbitsOffs;
    nbits_SB = TPos->nbitsSB;
    nbits_Offs = TPos->nbitsOffs;
    
    if (Count > 0) {
    for (i=Count-1; i; i--) {
       id = OccArray[i];
       if ((id+1) > nOffset) OccArray[i] = Tlength - OffsetArray[i];
       else {
          p_aux4 = (id>>5)*nbits_SB;
          p_aux3 = id*nbits_Offs;
          OccArray[i] = (bitget_aux4() + bitget_aux3()) - OffsetArray[i];
       }
    }
    id = OccArray[0];
    if ((id+1) > nOffset) OccArray[0] = Tlength - OffsetArray[0];
    else {
       p_aux4 = (id>>5)*nbits_SB;
       p_aux3 = id*nbits_Offs;
       OccArray[0] = (bitget_aux4() + bitget_aux3()) - OffsetArray[0];
    }
    OccArray = realloc(OccArray, Count*sizeof(uint));
    *occ = OccArray;
    }
    return 0; // no errors yet
 } 
 
 


void extract_text(void *index, ulong from, ulong to, byte *snippet, 
                  ulong *snippet_length)
 {
    ulong snpp_pos, posaux, idaux;
    lzindex I = *(lzindex *)index;
    uint idfrom,idto;
    trieNode p;
    
    idfrom = getphrasePosition(I.TPos, from);
    idto   = getphrasePosition(I.TPos, to+1);
    *snippet_length = to-from+1;
    
    posaux = getPosition(I.TPos, idto+1)-1;
    p = NODE(I, idto);
    while (p&&(posaux > to)) {
       p = parentLZTrie(I.fwdtrie, p);
       posaux--;
    }
    snpp_pos = (*snippet_length)-1;
    for (idaux = idto; idaux != idfrom;) {
       while (p) {
          snippet[snpp_pos--] = letterLZTrie(I.fwdtrie, p);
          p = parentLZTrie(I.fwdtrie, p);   
       }
       p = NODE(I, --idaux);
    }
    if (idfrom != idto) posaux = getPosition(I.TPos, idfrom+1)-(from!=0);
    while (p && (posaux >= from)) {
       snippet[snpp_pos--] = letterLZTrie(I.fwdtrie, p);
       p = parentLZTrie(I.fwdtrie, p);
       posaux--;
    }
 }



int extract(void *index, ulong from, ulong to, byte **snippet, 
            ulong *snippet_length)
 {
    ulong text_length;
    lzindex I = *(lzindex *)index;
    
    get_length(index, &text_length);
    if (to > (text_length-1)) to = text_length-1;
    if (from > (text_length-1)) from = text_length-1;

    *snippet_length = to-from+1;
    *snippet = malloc(*snippet_length);
    extract_text(&I,from,to,*snippet,snippet_length);
    return 0; // no errors yet
 }
 
 
   
int display(void *index, byte *pattern, ulong plength, ulong numc,
            ulong *numocc, byte **snippet_text, ulong **snippet_lengths)
 {
    ulong *occ, i, k, from, to, text_length;
    lzindex I = *(lzindex *)index;
    
    get_length(index, &text_length);
    locate(&I,pattern,plength,(ulong **)&occ,numocc);
    k = *numocc;
    (*snippet_text)    = malloc(sizeof(byte)*k*(plength+2*numc));
    for (i = 0; i < k; i++) {
       from = occ[i] - numc;
       to   = occ[i] + plength -1 + numc;
       if (numc >= occ[i]) from = 0;
       if (to > (text_length-1)) to = text_length-1;
       if (from > (text_length-1)) from = text_length-1;
       extract_text(&I,from,to,&((*snippet_text)[i*(plength+2*numc)]),&(occ[i]));
    }
    *snippet_lengths = occ;
    return 0; // no errors yet
 }

