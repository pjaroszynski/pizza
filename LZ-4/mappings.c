
#include "mappings.h"


inline uint NODE(lztrie T, uint id)
 {
    if (!id) return 0;
    return select_0(T->pdata->bdata, inversePerm(T->ids, id)+1);
 }

inline uint RNODE(revtrie RT, uint id)
 {
    if (!id) return 0;
    return getnodeRevTrie(RT, inversePerm(RT->rids, id));
 }

inline uint IDS(lztrie T, uint pos)
 {
    return getelemPerm(T->ids, pos);
 }

inline uint RIDS(revtrie RT, uint rpos)
 {
    return getelemPerm(RT->rids, rpos);
 }

