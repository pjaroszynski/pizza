
#include "mappings.h"

inline uint Rev(lzindex I, uint rpos)
 {
    return mapto(I.Rev,rpos);
 }

inline uint NODE(lzindex I, uint id)
 {
    if (!id) return 0;
    return Rev(I, getposRevTrie(I.bwdtrie,leftrankRevTrie(I.bwdtrie, RNODE(I,id))));
 }

inline uint RNODE(lzindex I, uint id)
 {
    return mapto(I.rmap,id);
 }

inline uint IDS(lzindex I, uint pos)
 {
    return rthLZTrie(I.fwdtrie, pos);
 }

inline uint RIDS(lzindex I, uint rpos)
 {
    return IDS(I, leftrankLZTrie(I.fwdtrie,Rev(I,rpos)));
 }
