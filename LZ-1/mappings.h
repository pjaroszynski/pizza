
// implementation of the different mappings of LZ-index

#include "bitmap.h"
//#include "revtrie.h"
//#include "lztrie.h"
#include "lzindex.h"

inline uint NODE(lzindex I, uint id);

inline uint RNODE(lzindex I, uint id);

inline uint IDS(lzindex I, uint pos);

inline uint RIDS(lzindex I, uint rpos);