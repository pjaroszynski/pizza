 
// Search module


#ifdef QUERYREPORT
int OCC1 = 0,OCC2 = 0,OCC3 = 0;
#endif


// reports occurrences,
// if node is not computed (==NULLT), report will take care by doing Map
// time is O(1) for levels 0 and 1,

// 0 -> count, 1 -> show positions
static uint reportLevel = 1; 
static uint Count, SIZE_OCCARRAY;
static ulong *OccArray;
static byte blockend = '\n';
position TPos;

#define INIT_OCCARRAY 5000
#define ALPHA 0.7


static void report(lztrie trie, nodemap map, trieNode node, uint id, 
                   uint delta, uint m)
 { 
    Count++;
    if (reportLevel == 0) return;
    // level >= 1
    if (reportLevel >= 1) {
       OccArray[Count-1] = getPosition(TPos, id+1)-delta;
       if (Count == SIZE_OCCARRAY) {
          SIZE_OCCARRAY = (SIZE_OCCARRAY+1)/ALPHA;
          OccArray       = realloc(OccArray, SIZE_OCCARRAY*sizeof(ulong));
       }    
       return;
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

static trieNode *fwdmatrix(lztrie T, byte *pat, uint m, uint **fwdid)
 { 
    trieNode *answ, cur, new;
    uint i,j,ptr, *id;
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

static trieNode *bwdmatrix(revtrie T, nodemap rmap, byte *pat, uint m, 
                           uint *fwdid, lztrie LZT)
 { 
    int i,j,k;
    trieNode cur,cur2,lzcur,lzcur2,*answ;
    byte c1,c2;
    answ = malloc(m*sizeof(trieNode));
    // answ[j] = node(pat[0..j]);
    for (j=0; j<m; j++) { 
       for (i=0;i<=j;i++) if (fwdid[m*i+j] != ~0) break;
       if (i <= j) {// i is last nonempty node in the path
          cur = mapto(rmap,fwdid[m*i+j]);
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

static void reportType1(lztrie fwdtrie, revtrie bwdtrie, nodemap map,
                        trieNode *fwd, trieNode *bwd, uint *fwdid, uint m)
 { 
    uint from,to,id,oldid,siz,pos,delta;
    trieNode lzcur,new;
    
    if (bwd[m-1] == NULLT) return; // there is not any LZ78 phrase ending with P
    from  = leftrankRevTrie(bwdtrie,bwd[m-1]);
    to    = rightrankRevTrie(bwdtrie,bwd[m-1]);
    oldid = 0;
    while (from <= to) {
       id = rthRevTrie(bwdtrie,from++);
       if (oldid != id) {// replicas due to empty nodes in revtrie
          oldid = id;
          lzcur = mapto(map,id);
          siz   = subtreesizeLZTrie(fwdtrie,lzcur);
          if (reportLevel == 0) Count += siz; // just count, faster
          else { 
             pos   = leftrankLZTrie(fwdtrie,lzcur); 
             delta = m;
             while (true) {
                report(fwdtrie,map,lzcur,id,delta,m);
                if (!--siz) break;
                id    = rthLZTrie(fwdtrie,++pos);
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

static void reportType2(lztrie fwdtrie, revtrie bwdtrie, nodemap map,
                        nodemap rmap, trieNode *fwd, trieNode *bwd, uint m)
 { 
    uint i,ffrom,fto,bfrom,bto;
    uint id,pid;
    trieNode cur,open,close;
    
    for (i=1; i<m; i++)
       if ((fwd[i*m+m-1] != NULLT) && (bwd[i-1] != NULLT)) {
          ffrom = leftrankLZTrie(fwdtrie,fwd[i*m+m-1]);
          fto   = rightrankLZTrie(fwdtrie,fwd[i*m+m-1]);
          bfrom = leftrankRevTrie(bwdtrie,bwd[i-1]);
          bto   = rightrankRevTrie(bwdtrie,bwd[i-1]);
          if (fto-ffrom <= bto-bfrom) {
             open = bwd[i-1];
             close = findclose (bwdtrie->pdata,open);
#ifdef QUERYREPORT
             WORK2 += fto-ffrom+1;
#endif 
             while (ffrom <= fto) {
                id  = rthLZTrie(fwdtrie,ffrom++);
                cur = mapto(rmap,id-1);
                if ((open <= cur) && (cur <= close))
                   report(fwdtrie,map,NULLT,id-1,i,m);
             }
          }
          else {
             open  = fwd[i*m+m-1];
             close = findclose(fwdtrie->pdata,open);
             pid   = ~0;
#ifdef QUERYREPORT
             WORK2 += bto-bfrom+1;
#endif 
             while (bfrom <= bto) {
                id  = rthRevTrie(bwdtrie,bfrom++);
                if (id == pid) continue; pid = id;
                cur = mapto (map,id+1);
                if ((open <= cur) && (cur <= close))
                   report(fwdtrie,map,NULLT,id,i,m);
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

static void reportType3(lztrie fwdtrie, revtrie bwdtrie,
                        nodemap map, nodemap rmap,
                        trieNode *fwd, uint *fwdid, trieNode *bwd, 
                        byte *pat, uint m)
 {
    uint i,j,k,f,t,nt,ok;
    trieNode node;
    helem *table;
    uint tsize;
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
                if (!ancestorLZTrie(fwdtrie,node,mapto(map,k+1))) continue;
             }
             // ok, test on block k+1 passed
             if (f == 0)  // left test not needed
                report(fwdtrie,map,fwd[i*m+j],ok,j-i+1,m);
             else {
                node = bwd[f-1];
                if (node == NULLT) continue; // rest does not exist 
                if (!ancestorRevTrie(bwdtrie,node,mapto(rmap,ok-1))) 
                   continue;
                report(fwdtrie,map,NULLT,ok-1,f,m);
             }
          }
    free(table);
 }


int count(void *index, byte *pattern, ulong length, ulong* numocc)
 {
    lzindex I = *(lzindex *)index;
    trieNode *fwd,*bwd;
    uint *fwdid;
    
    reportLevel = 0;
    Count       = 0;    
    fwd = fwdmatrix(I.fwdtrie,pattern,length,&fwdid);
    bwd = bwdmatrix(I.bwdtrie,I.rmap,pattern,length,fwdid,I.fwdtrie);
    reportType1(I.fwdtrie,I.bwdtrie,I.map,fwd,bwd,fwdid,length);
    if (length > 1)
       reportType2(I.fwdtrie,I.bwdtrie,I.map,I.rmap,fwd,bwd,length);
    if (length > 2)
       reportType3(I.fwdtrie,I.bwdtrie,I.map,I.rmap,fwd,fwdid,bwd,pattern,length);
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
    
    reportLevel   = 1;
    TPos          = I.TPos;
    Count         = 0;
    SIZE_OCCARRAY = INIT_OCCARRAY;
    OccArray      = malloc(sizeof(ulong)*SIZE_OCCARRAY);
    fwd = fwdmatrix(I.fwdtrie,pattern,length,&fwdid);
    bwd = bwdmatrix(I.bwdtrie,I.rmap,pattern,length,fwdid,I.fwdtrie);
    reportType1(I.fwdtrie,I.bwdtrie,I.map,fwd,bwd,fwdid,length);
    if (length > 1)
       reportType2(I.fwdtrie,I.bwdtrie,I.map,I.rmap,fwd,bwd,length);
    if (length > 2)
       reportType3(I.fwdtrie,I.bwdtrie,I.map,I.rmap,fwd,fwdid,bwd,pattern,length);
    free(fwd); 
    free(bwd); 
    free(fwdid);
    *numocc = Count;
    *occ    = malloc(sizeof(ulong)*(*numocc));
    occaux = *occ;
    for (i=0; i<Count; i++) occaux[i] = OccArray[i];
    free(OccArray);
    return 0; // no errors yet
 } 
 
 


void extract_text(void *index, ulong from, ulong to, byte *snippet, 
                  ulong *snippet_length)
 {
    ulong snpp_pos, posaux, idaux;
    lzindex I = *(lzindex *)index;
    uint idfrom,idto;
    trieNode p;
    
    idfrom = getphrasePosition(TPos, from);
    idto   = getphrasePosition(TPos, to);
    *snippet_length = to-from+1;
    
    posaux = getPosition(TPos, idto+1)-1;
    p = mapto(I.map, idto);
    while (p && (posaux > to)) {
       p = parentLZTrie(I.fwdtrie, p);
       posaux--;
    }
    snpp_pos = (*snippet_length)-1;
    for (idaux = idto; idaux != idfrom;) {
       while (p) {
          snippet[snpp_pos--] = letterLZTrie(I.fwdtrie, p);
          p = parentLZTrie(I.fwdtrie, p);   
       }
       p = mapto(I.map, --idaux);
    }
    if (idfrom != idto) posaux = getPosition(TPos, idfrom+1)-(from!=0);
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
    TPos     = I.TPos;
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
