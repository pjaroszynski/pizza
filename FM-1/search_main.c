/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
   search_main.c   
   P. Ferragina & G. Manzini, 29-sept-00 

   these are the main routines for searching in bwi files
   which are called by bwsearch, bwxml, bwrand and bwhuffw
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "common.h"

// -------------------------------------------------------------------
// ----------- Global variables & Macros for searching   -------------
// -------------------------------------------------------------------
#define EOF_shift(n) (n < s->bwt_eof_pos) ? n+1 :  n

static uchar Bool_map_sb[256]; // info alphabet for superbucket to be read
uchar Inv_map_sb[256];         // inverse map for the current superbucket
int Alpha_size_sb;             // actual size of alphabet in superbucket
static uchar Bool_map_b[256];  // info alphabet for the bucket to be read
static uchar Inv_map_b[256];   // inverse map for the current superbucket
static int Alpha_size_b;       // actual size of alphabet in bucket 
static uchar Mtf[256]; // stores MTF-picture of bucket to be decompressed


// -------------------------------------------------------------------
// -------------              Caching System          ----------------
// -------------------------------------------------------------------
static uchar **Cache_of_buckets;  // Array of pointers to uncompressed buckets
static int Num_bucs_in_cache;     // # buckets into cache
static int *LRU_queue_bucs;  // LRU ordered array of indices of uncmp. buckets
static int LRU_index;             // index of the least recently used bucket
static int Max_cached_buckets;    // Maximum number of cacheable buffers
static int Use_bwi_cache=0;       // !=0 if some cache is used
double Cache_percentage;          // size of cache (% wrt uncompressed size)

/* **********************************************************
   This is the bwsearch procedure as described in the paper.  Given a
   patter p[] returns the number of occurrences in the uncompressed
   file, and stores in sp and ep the start/end positions of the
   pattern occurrences in the BWI.  Note that the EOF_shift macro
   defined below is required since in the BWT we remove the EOF symbol
   storing its position in bwt_eof_pos
   Note: 
   we are assuming that read_basic_prologue() has been already called!
   ********************************************************** */
int bwsearch(bwi_out *s, uchar *p, int len, int *sp, int *ep)
{
  int occ(bwi_out *,int,uchar),uint_read(void);
  int i;
  uchar c;

  /* ---- remap pattern ------ */
  assert(len>0);
  for(i=0;i<len;i++) {
    if(s->bool_char_map[p[i]]==0) return 0;  /* char not in file */
    p[i]=s->char_map[p[i]];                  /* remap char */
  }

  /* ---- get initial sp and ep values ----- */ 
  i=len-1;
  c=p[i];
  *sp=s->bwt_occ[c];
  if(c==s->alpha_size-1)
    *ep=s->text_size-1;
  else 
    *ep=s->bwt_occ[c+1]-1;   

  /* ----- main loop (see paper) ---- */
  while((*sp <= *ep) && (i>0)) {
    c=p[--i];
    *sp = s->bwt_occ[c]+occ(s, EOF_shift(*sp - 1), c);
    *ep = s->bwt_occ[c]+occ(s, EOF_shift(*ep), c)-1;
  }

  /* ----- inverse remap pattern -------- */
  for(i=0;i<len;i++) 
    p[i]=s->inv_char_map[p[i]];
  /* ----- return number of occurrences ----- */
  if(*ep < *sp) return 0;
  return(*ep - *sp + 1);
}

/* *********************************************************************
   read len chars before the position given by row using the LF mapping
   stop if the beginning of the file is encountered
   return the number of chars actually read 
   ********************************************************************* */
int go_back(int row, int len, char *dest, bwi_out *s)
{
  int written,curr_row,n;
  uchar c,c_sb;
  int occ_sb[256],occ_b[256];
  
  curr_row=row;
  for(written=0;written<len; ) {
    // fetches info from the header of the superbucket
    get_info_sb(s,curr_row,occ_sb);  
    // fetches occ into occ_b properly remapped and returns
    // the remapped code for occ_b of the char in the  specified position
    c = get_info_b(s,NULL_CHAR,curr_row,occ_b,WHAT_CHAR_IS);  
    assert(c < Alpha_size_sb);
    c_sb = Inv_map_sb[c];
    assert(c_sb < s->alpha_size);
    printf("%c\n",s->inv_char_map[c_sb]);

    dest[written++] = s->inv_char_map[c_sb]; // store char    
    n = occ_sb[c_sb] + occ_b[c] - 1;         // # of occ before curr_row
    curr_row = s->bwt_occ[c_sb] + n;         // get next row
    if(curr_row==s->bwt_eof_pos) break;    
    curr_row = EOF_shift(curr_row);        
  }
  return written;
}

/* *******************************************************
   write in dest[] the first "len" chars of "row" (unless
   the EOF is encountered).
   return the number of chars actually read 
   ******************************************************* */
int go_forw(int row, int len, char *dest, bwi_out *s)
{
  uchar get_firstcolumn_char(int row, bwi_out *s);
  int fl_map(int,uchar,bwi_out *);
  int written;
  uchar c;

  for(written=0;written<len; ) {
    c = get_firstcolumn_char(row,s);
    assert(c<s->alpha_size);
    dest[written++] = s->inv_char_map[c];
    // compute the first to last mapping
    row = fl_map(row,c,s);
    // adjust row to take the EOF symbol into account
    if(row<=s->bwt_eof_pos) row -= 1;
    if(row<0) break;   // end of file     
  }
  return written;
}

/* *********************************************************************
   read num_word words before the position given by row using the LF mapping
   stop if the beginning of the file is encountered
   return the number of chars actually read 
   ********************************************************************* */
int go_back_word(int row, int num_word, char *dest, bwi_out *s, int bsize)
{
  int written,curr_row,n,curr_word;
  uchar c,c_sb;
  int occ_sb[256],occ_b[256];
  
  curr_row=row;
  for(written=0,curr_word=0; curr_word<num_word; ) {
    // fetches info from the header of the superbucket
    get_info_sb(s,curr_row,occ_sb);  
    // fetches occ into occ_b properly remapped and returns
    // the remapped code for occ_b of the char in the  specified position
    c = get_info_b(s,NULL_CHAR,curr_row,occ_b,WHAT_CHAR_IS);  
    assert(c < Alpha_size_sb);
    c_sb = Inv_map_sb[c];
    assert(c_sb < s->alpha_size);
    dest[written] = s->inv_char_map[c_sb]; // store char
    if (dest[written] == ' ' || dest[written]=='\n') 
      curr_word++;                           // end of word found  
    if(++written == bsize)
      break;                                 // buffer is full. 
    n = occ_sb[c_sb] + occ_b[c] - 1;         // # of occ before curr_row
    curr_row = s->bwt_occ[c_sb] + n;         // get next row
    if(curr_row==s->bwt_eof_pos) break;    
    curr_row = EOF_shift(curr_row);        
  }
  return written;
}



/* *******************************************************
   write in dest[] the first "num_word" words of "row" (unless
   the EOF is encountered).
   return the number of chars actually read 
   ******************************************************* */
int go_forw_word(int row, int num_word, char *dest, bwi_out *s, int bsize)
{
  uchar get_firstcolumn_char(int row, bwi_out *s);
  int fl_map(int,uchar,bwi_out *);
  int written,curr_word;
  uchar c;

  for(written=0,curr_word=0;curr_word<num_word; ) {
    c = get_firstcolumn_char(row,s);
    assert(c<s->alpha_size);
    dest[written] = s->inv_char_map[c];
    if (dest[written] == ' ' || dest[written]=='\n') 
      curr_word++;                           // end of word found
    if(++written == bsize)
      break;                                 // buffer is full.   
    // compute the first to last mapping
    row = fl_map(row,c,s);
    // adjust row to take the EOF symbol into account
    if(row<=s->bwt_eof_pos) row -= 1;
    if(row<0) break;   // end of file     
  }
  return written;
}


/* *********************************************************************
   skip num_url urls before the position given by row using the LF mapping
   stop if the beginning of the file is encountered
   return the row which starts with a NULL and the last skipped url.
   ********************************************************************* */
int skip_url_back(int row, int num_url, bwi_out *s)
{
  int curr_row,n,curr_word;
  uchar c,c_sb,null_remap;
  int occ_sb[256],occ_b[256];
  

  curr_row=row;

  null_remap = s->char_map[0];

  for(curr_word=0; curr_word < num_url; ) {

    get_info_sb(s,curr_row,occ_sb);  
    c = get_info_b(s,NULL_CHAR,curr_row,occ_b,WHAT_CHAR_IS);  
    c_sb = Inv_map_sb[c];

    if (s->inv_char_map[c_sb] == '\0')
      curr_word++;                           // end of word found  

    n = occ_sb[c_sb] + occ_b[c] - 1;         // # of occ before curr_row
    curr_row = s->bwt_occ[c_sb] + n;         // get next row
 
    if(curr_row==s->bwt_eof_pos)
      break;    // beginning of file

    curr_row = EOF_shift(curr_row);        
 
 }

  return curr_row; // row starts with NULL
}

/* *********************************************************************
   get the url which precedes the suffix stored at the row "row".
   the url is stored reversed and with NULL at the end.
   return the number of chars actually read (included the last NULL).
   ********************************************************************* */
int get_url_back(int row, char *url_text, bwi_out *s, int bsize)
{
  int written,curr_row,n;
  uchar c,c_sb;
  uchar null_remap;
  int occ_sb[256],occ_b[256];  

  curr_row=row;

  null_remap = s->char_map[0];

  for(written=0;; ) { // just one url

    get_info_sb(s,curr_row,occ_sb);  
    c = get_info_b(s,NULL_CHAR,curr_row,occ_b,WHAT_CHAR_IS);  
    c_sb = Inv_map_sb[c];

    url_text[written]=s->inv_char_map[c_sb];

    if (url_text[written++] == '\0')  break; // beginning of url

    if(written == bsize) break;              // buffer is full. 

    n = occ_sb[c_sb] + occ_b[c] - 1;         // # of occ before curr_row
    curr_row = s->bwt_occ[c_sb] + n;         // get next row

    if(curr_row==s->bwt_eof_pos) {
      url_text[written++] = '\0';
      break;
    }    
    curr_row = EOF_shift(curr_row);        // manage the dropping of EOF in L

  }


  return written;
}


/* *****************************************************
   compute the first-to-last map using binary search
   ***************************************************** */
int fl_map(int row, uchar ch, bwi_out *s)
{
  int i, n, rank, first, last, middle;
  int occ_sb[256],occ_b[256];
  uchar c_b,c_sb;              // char returned by get_info
  uchar ch_b;

  // rank of c in first column 
  rank = 1 + row - s->bwt_occ[ch];
  // get position in the last column using binary search
  first=0; last=s->text_size;
  // invariant: the desired position is within first and last
  while(first<last) {
    middle = (first+last)/2;
    /* -------------------------------------------------------------
       get the char in position middle. As a byproduct, occ_sb
       and occ_b are initialized
       ------------------------------------------------------------- */
    get_info_sb(s,middle,occ_sb);   // init occ_sb[]  
    c_b = get_info_b(s,NULL_CHAR,middle,occ_b,WHAT_CHAR_IS); // init occ_b[] 
    assert(c_b < Alpha_size_sb);
    c_sb = Inv_map_sb[c_b];           
    assert(c_sb < s->alpha_size);  // c_sb is the char in position middle 
    /* --------------------------------------------------------------
       count the # of occ of ch in [0,middle]
       -------------------------------------------------------------- */
    if(Bool_map_sb[ch]==0)
     n=occ_sb[ch];          // no occ of ch in this superbucket
    else {
      ch_b=0;                        // get remapped code for ch
      for(i=0;i<ch;i++)                  
	if(Bool_map_sb[i]) ch_b++;
      assert(ch_b<Alpha_size_sb);
      n = occ_sb[ch] + occ_b[ch_b];  // # of occ of ch in [0,middle]
    }
    /* --- update first or last ------------- */
    if(n>rank)
      last=middle;
    else if (n<rank)
      first=middle+1;
    else {                      // there are exactly "rank" c's in [0,middle]
      if(c_sb==ch) 
	first=last=middle;      // found!
      else 
	last=middle;            // restrict to [first,middle)
    }
  }
  // first is the desired row in the last column
  return first;
}


/* **************************************************
   return the first character of a given row.
   This routine can be improved using binary search!
   ************************************************** */
uchar get_firstcolumn_char(int row, bwi_out *s)
{
  int i;

  for(i=1;i<s->alpha_size;i++)
    if(row<s->bwt_occ[i])
      return i-1;
  return s->alpha_size-1;
}


/* **************************************************
   return the last character of a given row.
   Keep attention to the row of EOF, not dealed with.
   ************************************************** */
uchar get_lastcolumn_char(int row, bwi_out *s)
{
  int occ_sb[256],occ_b[256];
  uchar c,c_sb;

  get_info_sb(s,row,occ_sb);  
  c = get_info_b(s,NULL_CHAR,row,occ_b,WHAT_CHAR_IS);  
  c_sb = Inv_map_sb[c];
  return(s->inv_char_map[c_sb]);
}



/* **********************************************************
   This is the procedure to retrieve the locations of the pattern
   occurrences lying at row number "row". We backtrack on the text
   using the LF[] mapping until the special character, "marked" during 
   the compression process, is encountered. For it we have its explicit
   text location. So that we can derive the location of the pattern
   occurrence.
   Note that "row" here stands for the row number in the first 
   column! to get the row number in the last column we use the 
   EOF_shift as usual.
   ********************************************************** */
int get_occ_pos(bwi_out *s, int row)
{
  int occ(bwi_out *,int,uchar);
  void get_info_sb(bwi_out *,int,int *);
  uchar get_info_b(bwi_out *,uchar,int,int *,int);
  int my_fseek(FILE *f, long offset, int whence);
  int int_log2(int),bit_read(int);
  void init_bit_buffer(void);
  uchar c,c_sb;
  int i,len,occ_sb[256],occ_b[256],n,curr_row;
  
  assert(s->skip!=0);
  curr_row = row;
  len = int_log2(s->text_size); 

  /* ---- get character in L[] for the current row ----- */ 
  for(i=0; ;i++){

    if (curr_row == s->bwt_eof_pos)   // Manage  special case
      return i;

    // fetches info from the header of the superbucket
    get_info_sb(s,EOF_shift(curr_row),occ_sb);  

    // fetches occ into occ_b properly remapped and returns
    // the remapped code for occ_b of the char in the 
    // specified position
    c = get_info_b(s,NULL_CHAR,EOF_shift(curr_row),occ_b,WHAT_CHAR_IS);  
    assert(c < Alpha_size_sb);

    // Inverse Remapping: bucket --> superbucket --> remapped_text
    c_sb = Inv_map_sb[c];
    assert(c_sb < s->alpha_size);

    // Compute # occ of character c before the position given above
    n = occ_sb[c_sb] + occ_b[c] - 1;

    // If the row is marked then its starting_pos is retrieved
    if((c_sb == s->chosen_char) && ((n % s->skip) == 0)){
      n = n/s->skip;     // Rescale
      my_fseek(Infile, ((n*len)/8) + Start_prologue_occ,SEEK_SET);
      init_bit_buffer();
      if((n*len) % 8)
	bit_read( (n*len) % 8);  // skip bits not belonging to this pos  
      return(bit_read(len) + i);  // read pos-representation
    } 
    else { // Otherwise we jump back
      curr_row = s->bwt_occ[c_sb] + n;
    }
  }
}


/* **********************************************************
   This is the procedure to retrieve the locations of the url
   occurrence lying at row number "row". We backtrack on the text
   using the LF[] mapping until a NULL char is encountered. In case is
   a maked char, we have its rank and thus we use it to compute the
   rank of the input url; otherwise we increase a counter and proceed
   backward.  Note that "row" here stands for the row number in the
   first column! to get the row number in the last column we use the
   EOF_shift as usual. 
***************************************************************/
int get_url_rank(bwi_out *s, int row)
{
  void get_info_sb(bwi_out *,int,int *);
  uchar get_info_b(bwi_out *,uchar,int,int *,int);
  int my_fseek(FILE *f, long offset, int whence);
  int int_log2(int),bit_read(int);
  void init_bit_buffer(void);
  uchar c,c_sb;
  int len,occ_sb[256],occ_b[256],n,curr_row;
  int null_occ,marked_chars,rank_size_bytes;
  int Start_prologue_ranks;
  uchar null_remap;
  int count=1;

  if (s->skip==0)
    fatal_error("s->skip is zero in get_id_byrow\n");

  // Check if it is the first url
  if(row==s->bwt_eof_pos)
    return 1;

  null_occ=s->bwt_occ[1];  // # occurrences of NULL
  null_remap = s->char_map[0]; // remapped NULL
  len = int_log2(null_occ);  // bits to read for every occ

  // Compute the #marked rows
  if(null_occ % s->skip)
    marked_chars = null_occ/s->skip + 1;
  else
    marked_chars = null_occ/s->skip;

  // Compute the starting byte of marked row pos
  rank_size_bytes = marked_chars * len;
  if (rank_size_bytes % 8)
    { rank_size_bytes= rank_size_bytes/8+1; }
  else 
    { rank_size_bytes= rank_size_bytes/8; }
  Start_prologue_ranks = Start_prologue_occ + rank_size_bytes;


  // Current line is not the removed one, we checked above
  // Keep track of the removed NULL
  curr_row=EOF_shift(row);

  for( ; ; ) {

    // Retrieve the char on the last column
    get_info_sb(s,curr_row,occ_sb);  
    c = get_info_b(s,NULL_CHAR,curr_row,occ_b,WHAT_CHAR_IS);  
    c_sb = Inv_map_sb[c];

    // Determine the row on F where that char lie
    n = occ_sb[c_sb] + occ_b[c] - 1;         // # of occ before curr_row
    curr_row = s->bwt_occ[c_sb] + n;         // get next row

    // Check if beginning of file reached
    if(curr_row==s->bwt_eof_pos)
      return count;

    // Count the number of encountered NULL
    if ( c_sb == null_remap ){

      count++;

      // If the row is marked then its starting_pos is retrieved
      if ((n % s->skip) == 0){

	n = n/s->skip;     // Rescale
	my_fseek(Infile, ((n*len)/8) + Start_prologue_ranks,SEEK_SET);
	init_bit_buffer();
	if((n*len) % 8)
	  bit_read( (n*len) % 8);  // skip bits not belonging to this pos  
	
	return(bit_read(len) + count);// read pos-representation
      } 
    } 

    curr_row = EOF_shift(curr_row);        
  }

  return -1; // We should never execute this
}




/* ***********************************************************
   compute # of occ of ch among the positions 0 ... k of the bwt
   *********************************************************** */
int occ(bwi_out *s, int k, uchar ch)
{
  void get_info_sb(bwi_out *,int,int *);
  uchar get_info_b(bwi_out *,uchar,int,int *,int);
  uchar new_ch;
  int i;

  // Count for ch occurrences in the superbucket and bucket
  // containing the position k.
  int occ_sb[256],occ_b[256];  
                                  
  get_info_sb(s,k,occ_sb); // get occ  of all chars in prev superbuckets     

  if(Bool_map_sb[ch]==0)
    return(occ_sb[ch]);          // no occ of ch in this superbucket
 
  new_ch=0;                        // get remapped code for ch
  for(i=0;i<ch;i++)                  
    if(Bool_map_sb[i]) new_ch++;

  // Account for the occurrences of ch from the beginning
  // of the superbucket up to position k
  get_info_b(s,new_ch,k,occ_b,COUNT_CHAR_OCC);

  // Return the proper count of the occurrences of ch in BWT[0..k]
  return(occ_sb[ch] + occ_b[new_ch]);
 
}


/* ********************************************************
   Read informations from the header of the superbucket.
   "pos" is a position in the last column (that is the bwt). We
   we are interested in the information for the superbucket 
   containing "pos".
   Initializes the data structures:
        Inv_map_sb, Bool_map_sb, Alpha_size_sb
   Returns initialized the array occ[] containing the number 
   of occurrences of the chars in the previous superbuckets.
   ******************************************************** */
void get_info_sb(bwi_out *s, int pos, int *occ)
{
  void init_bit_buffer(void);
  int bit_read(int);
  int my_fseek(FILE *f, long offset, int whence);
  int uint_read(void); 
  int offset,i,size,sb;

  // ----- initialize the data structures for I/O-reading
  init_bit_buffer();
  sb =  pos/Bucket_size_lev1;        // superbucket containing pos
  if(pos>=s->text_size)
    fatal_error("Invalid pos! (get_info_sb)\n");
  assert(sb < Num_bucs_lev1);

  offset = Start_prologue_info_sb;

  // --------- go to the superbucket header 
  if(sb>0) {
    // skip bool map in previous superbuckets
    offset += sb*((s->alpha_size+7)/8);   
    // skip occ in previous superbuckets (except the first one)
    offset += (sb-1)*(s->alpha_size*sizeof(int));   
  }

  my_fseek(Infile,offset,SEEK_SET);          // position of sb header

  // --------- get bool_map_sb[]
  for(i=0;i<s->alpha_size;i++)         
    Bool_map_sb[i]=bit_read(1);

  // -----  compute alphabet size
  Alpha_size_sb=0;
  for(i=0;i<s->alpha_size;i++)
    if(Bool_map_sb[i]) Alpha_size_sb++;

  // ----- Invert the char-map for this superbucket
  for(i=0,size=0;i<s->alpha_size;i++) 
    if(Bool_map_sb[i]) 
      Inv_map_sb[size++]= (uchar) i;
  assert(size==Alpha_size_sb);

  // ---- for the first sb there are no previous occurrences  
  if (sb==0) {
    for(i=0; i<s->alpha_size; i++)
      occ[i] = 0;
  } else {         
    // ------- otherwise compute # occ_map in previous superbuckets
    if((s->alpha_size % 8) != 0)
      bit_read(8-s->alpha_size % 8); // restore byte alignement     
    for(i=0;i<s->alpha_size;i++)          
      occ[i] = uint_read();          
  }
}


/* **********************************************************************
   Read informations from the header of the bucket and decompresses the 
   bucket if needed. Initializes the data structures:
        Inv_map_b, Bool_map_b, Alpha_size_b
   Returns initialized the array occ[] containing the number 
   of occurrences of all the chars since the beginning of the superbucket
   Explicitely returns the character (remapped in the alphabet of 
   the superbucket) occupying the absolute position pos.
   The decompression of the bucket when ch does not occur in the 
   bucket (because of Bool_map_b[ch]=0) is not always carried out. 
   The parameter "flag" setted to COUNT_CHAR_OCC indicates that we want 
   to count the occurreces of ch; in this case when Bool_map_b[ch]==0
   the bucket is not decompressed.
   When the flag is setted to WHAT_CHAR_IS, then ch is not significant 
   and we wish to retrieve the character in position k. In this case
   the bucket is always decompressed.
   ********************************************************************** */
uchar get_info_b(bwi_out *s, uchar ch, int pos, int *occ, int flag)
{
  int int_log2(int),read7x8(void);
  void init_bit_buffer(void);
  uchar get_b_hier(int,int *);
  uchar get_b_multihuf(int,int *);
  int my_fseek(FILE *f, long offset, int whence);
  int bit_read(int);
  int buc_start_pos,len,buc,size,bits_x_char,mtf_size;
  int i;
  uchar ch_in_pos;

  buc = pos/Bucket_size_lev2;       // bucket containing pos
  assert(buc < (s->text_size + Bucket_size_lev2 - 1)/Bucket_size_lev2);

  len = (int_log2(s->text_size)+7)/8;   // variable byte-length repr.

  // read bucket starting position
  my_fseek(Infile,Start_prologue_info_b + len * buc,SEEK_SET);
  init_bit_buffer();
  buc_start_pos = bit_read(len * 8);

  // move to the beginning of the bucket
  my_fseek(Infile,buc_start_pos,SEEK_SET);
  init_bit_buffer();

  // ---- Initialize properly the occ array ----
  if((buc % (Bucket_size_lev1/Bucket_size_lev2)) == 0 ){
    for(i=0; i<Alpha_size_sb; i++)
      occ[i] = 0;
  } else {
    for(i=0;i<Alpha_size_sb;i++) // read all the occ table 
      occ[i] = read7x8(); 
  }
  
  // ---- get bool char map -----
  for(i=0;i<Alpha_size_sb;i++) // read all the bit_char_map    
    Bool_map_b[i]= bit_read(1);     

  // --- get bucket alphabet size and the code of ch in this bucket 
  Alpha_size_b=0;
  for(i=0;i<Alpha_size_sb;i++)
    if(Bool_map_b[i]) Alpha_size_b++;   // alphabet size in the bucket

  // --- if no occ of this char in the bucket then skip everything
  if((flag == COUNT_CHAR_OCC) && (Bool_map_b[ch]==0)) 
    return((uchar) 0);  //dummy return

  // ----- Invert the char-map for this bucket
  for(i=0,size=0;i<Alpha_size_sb;i++) 
    if(Bool_map_b[i]) 
      Inv_map_b[size++]= (uchar) i;
  assert(size==Alpha_size_b);

  /* --------- read initial status of the mtf list ---------- */
  mtf_size = MIN(Mtf_save,Alpha_size_b);
  bits_x_char = int_log2(Alpha_size_b);   // read initial mtf list
  for(i=0;i<mtf_size;i++) {
    Mtf[i] = bit_read(bits_x_char);
  }

  /* ----- decompress and count CH occurrences on-the-fly ------ */
  switch (Type_compression) 
    {
    case ARITH:   // ---- Arithmetic compression of the bucket -----
      fatal_error("Arithmetic coding no longer available -get_info_b-\n");
      exit(1);
    case HIER3:  // ---- three-leveled model: Fenwick's proposal -----
      ch_in_pos = get_b_hier(pos,occ);
      break;
    case UNARY: // ---- Unary compression of mtf-ranks with escape -----
      fatal_error("Unary coding no longer available -get_info_b-\n");
      exit(1);
    case MULTIH: // ---- multihuffman compression of mtf-ranks -----
      ch_in_pos = get_b_multihuf(pos,occ);
      break;
    default:
      fprintf(stderr,"Compression algorithm unknown! -get_info_b-\n");
      exit(1);
    }
  return(ch_in_pos);   // char represented in [0.. Alpha_size.sb-1]
}


#if 0
/* ************************************************************
   Arithmetic compressed bucket. Update the array occ[] summing
   up all occurrencs of the chars in its prefix preceding
   the absolute position k.  Note that ch is a bucket-remapped char.
   ************************************************************ */
uchar get_b_arith(int k, int *occ)
{
  void init_bit_buffer(void);
  void out_of_mem(char *);
  int int_log2(int);
  int unrle_only(uchar *, int, uchar *);
  int adapt_arith_decoder(uchar *,int);
  int bit_read(int);
  int i,j,h,rank,mtf_size,bpos;
  int mtf_seq_len, rle_len;
  uchar next = NULL_CHAR;
  uchar *rle, *mtf_seq; 


  mtf_size = MIN(Mtf_save,Alpha_size_b);
  bpos = k % Bucket_size_lev2;

  init_bit_buffer();  // align to the byte

  /* ------ uncompress arithmetic code -------------------- */
  rle = (uchar *) malloc(3 * Bucket_size_lev2 * sizeof(uchar));
  if(rle == NULL) out_of_mem("get_b_occ_arith");

  // Reads from file. Here +2 is a safe term for 0-bucket
  rle_len = adapt_arith_decoder(rle,Alpha_size_b+2);

  /* ------ unrle -------------------- */
  mtf_seq = (uchar *) malloc(3 * Bucket_size_lev2 * sizeof(uchar));
  if(mtf_seq == NULL) out_of_mem("get_b_occ_arith");
  mtf_seq_len = unrle_only(rle, rle_len, mtf_seq);
  assert(mtf_seq_len>bpos);

  /* ------ search or count -------------------- */
  for(j=0,h=0;j<=bpos; ) {
    rank = mtf_seq[h++];
    assert(rank<=mtf_size); 

    if(rank == 0) {
      next=Mtf[rank];     // decode mtf rank 
    }
    else if(rank<mtf_size) {
      next=Mtf[rank];     // decode mtf rank 
      for(i=rank;i>0;i--)    // update mtf list
        Mtf[i]=Mtf[i-1];
      Mtf[0]=next;            // update mtf[0]
    }
    else {                         // rank==mtf_size
      next=mtf_seq[h++];  // get char from file
      for(i=mtf_size-1;i>0;i--)    // update mtf
	Mtf[i]=Mtf[i-1];
      Mtf[0]=next;             // update mtf[0] 
    }

    occ[Inv_map_b[next]]++; // update occ counting for next
    j++;

  }

  free(mtf_seq);
  free(rle);

  return(Inv_map_b[next]);      // returning char occupying absolute position k

}
#endif


/* ************************************************************
   Hierarchically compressed bucket. Update the array occ[] summing
   up all occurrencs of the chars in its prefix preceding
   the absolute position k.  Note that ch is a bucket-remapped char.
   ************************************************************ */
extern int Mtf10log;
uchar get_b_hier(int k, int *occ)
{
  int int_log2(int);
  void init_bit_buffer(void);
  int bit_read(int);
  int decode_hierarchical(void);
  int j,i,rle,bits_x_char,rank,mtf_size,bpos;
  uchar next = NULL_CHAR; 

  mtf_size = MIN(Mtf_save,Alpha_size_b);
  bits_x_char = int_log2(Alpha_size_b); 
  bpos = k % Bucket_size_lev2;
  if(Mtf_save>10)
    Mtf10log = int_log2(mtf_size-10);  // constant used by decode_hierarchical

  /* ------ uncompress -------------------- */
  for(j=0;j<=bpos; ) {
    rank = decode_hierarchical();
    assert(rank<=mtf_size); 

    if(rank==0) {
      next=Mtf[0];       // rank zero: top of mtf list
      rle = bit_read(8);      // get # of following zeroes
      occ[Inv_map_b[next]] += MIN(1 + rle,1+bpos-j);
      j += 1+rle;
    }
    else if(rank<mtf_size) {
      next=Mtf[rank];     // decode mtf rank 
      for(i=rank;i>0;i--)    // update mtf list
        Mtf[i]=Mtf[i-1];
      Mtf[0]=next;            // update mtf[0]
      occ[Inv_map_b[next]]++; // update occ for next
      j++;                    // update j
    }
    else {                         // rank==mtf_size
      next = bit_read(bits_x_char);  // get char from file
      for(i=mtf_size-1;i>0;i--)    // update mtf
	Mtf[i]=Mtf[i-1];
      Mtf[0]=next;             // update mtf[0] 
      occ[Inv_map_b[next]]++; // update occ for next
      j++;                     // update j
    }
  }
  return(Inv_map_b[next]);   // return character occupying absolute position k
}

#if 0
/* **********************************************************************
   Unary compressed bucket. Update the array occ[] summing
   up all occurrencs of the chars in its prefix preceding
   the absolute position k.  Note that ch is a bucket-remapped char.
   ********************************************************************** */
uchar get_b_unary(int k, int *occ)
{
  static __inline__ int decode_unary(void);
  int bit_read(int),int_log2(int);
  int i,j,rle,rank,bits_x_char,bpos,mtf_size;
  uchar next = NULL_CHAR;

  mtf_size = MIN(Mtf_save,Alpha_size_b);
  bits_x_char = int_log2(Alpha_size_b); 
  bpos = k % Bucket_size_lev2;

  /* ------ uncompress and count -------------------- */
  for(j=0;j<=bpos; ) {
    rank = decode_unary();
    assert(rank<=mtf_size); 
    if(rank==0) {
      next=Mtf[0];            // rank zero: top of mtf list
      rle = bit_read(8);      // get # of following zeroes
      occ[Inv_map_b[next]] += MIN(1 + rle,1+bpos-j);
      j += 1+rle;
    }
    else if(rank<mtf_size) {
      next=Mtf[rank];     // decode mtf rank 
      for(i=rank;i>0;i--)    // update mtf list
        Mtf[i]=Mtf[i-1];
      Mtf[0]=next;            // update mtf[0]
      occ[Inv_map_b[next]]++;   // update occ for next
      j++;                    // update j
    }
    else {                         // rank==mtf_size
      next = bit_read(bits_x_char);  // get char from file
      for(i=mtf_size-1;i>0;i--)    // update mtf
	Mtf[i]=Mtf[i-1];
      Mtf[0]=next;             // update mtf[0] 
      occ[Inv_map_b[next]]++;  // update occ for next
      j++;                     // update j
    }
  }

  return(Inv_map_b[next]); // returning char at absolute position k

}
#endif

/* ************************************************************
   Multi-table-Huffman compressed bucket. Update the array occ[] 
   summing up all occurrencs of the chars in its prefix preceding
   the absolute position k. Note that ch is a bucket-remapped char.
   ************************************************************ */
uchar get_b_multihuf(int k, int *occ)
{
  void out_of_mem(char *);
  int int_log2(int);
  int multihuf_decompr(uchar *, int,int);
  int bit_read(int);
  void unmtf_unmap(uchar *,int);
  int i,bpos,mtf_seq_len;
  uchar *mtf_seq,*unmtf_bucket,char_returned; 

  assert(Mtf_save==256);    // only works with a full mtf list
  bpos = k % Bucket_size_lev2;

  // ----- first the case in which no cache is used ----
  if(Use_bwi_cache==0) {
    /* ------- decode bucket with bounded multiple huffman ----- */
    mtf_seq = (uchar *) malloc(Bucket_size_lev2*sizeof(uchar));
    if(mtf_seq==NULL) out_of_mem("-get_b_occ_multihuff- (mtf_seq)");
    mtf_seq_len = multihuf_decompr(mtf_seq,Alpha_size_b,bpos+1);
    assert(mtf_seq_len>bpos);
    assert(mtf_seq_len<=Bucket_size_lev2);

    /* --- The chars in the unmtf_bucket are already un-mapped --- */
    unmtf_unmap(mtf_seq,bpos+1);
    unmtf_bucket = mtf_seq;
    /* --- returning char at bwt-position k --> Inv[] not necessary --- */
    char_returned = unmtf_bucket[bpos];
    assert(char_returned < Alpha_size_sb);
    /* --- update occ[]array --- */  
    for(i=0; i<=bpos; i++)
      occ[unmtf_bucket[i]]++;   // unmtf_bucket contains chars un-mapped
    // deallocate mem and return 
    free(mtf_seq);
    return char_returned;
  }

  // ---- we are using some cache --------------
  if (Cache_of_buckets[k/Bucket_size_lev2] == NULL) {

    /* ------- decode bucket with multiple huffman codes ----- */
    mtf_seq = (uchar *) malloc(Bucket_size_lev2*sizeof(uchar));
    if(mtf_seq==NULL) out_of_mem("-get_b_occ_multihuff- (mtf_seq)");
    mtf_seq_len = multihuf_decompr(mtf_seq,Alpha_size_b,Bucket_size_lev2);
    assert(mtf_seq_len>bpos);

    /* --- The chars in the unmtf_bucket are already un-mapped --- */
    unmtf_unmap(mtf_seq,mtf_seq_len);
    unmtf_bucket = mtf_seq;

    /* ------- cache the uncompressed bucket ----- */
    Cache_of_buckets[k/Bucket_size_lev2] = unmtf_bucket;
    
    /* ------- store index of cached bucket in the LRU_queue ----- */
    LRU_queue_bucs[(LRU_index + Num_bucs_in_cache) % Num_bucs_lev2] 
      =  k/Bucket_size_lev2;
    Num_bucs_in_cache++;
  } 
  else {
    /* ------- retrieve the (uncompressed) cached bucket ----- */
    unmtf_bucket = Cache_of_buckets[k/Bucket_size_lev2];
  }

  /* --- returning char at bwt-position k --> Inv[] not necessary --- */
  char_returned = unmtf_bucket[bpos];
  assert(char_returned < Alpha_size_sb);
  /* --- update occ[]array --- */  
  for(i=0; i<=bpos; i++)
    occ[unmtf_bucket[i]]++;   // unmtf_bucket contains chars un-mapped

  /* ----- Update the Cache, the LRU_queue and free the memory ----- */
  if (Num_bucs_in_cache > Max_cached_buckets){

    // free the least recently used cached bucket
    free(Cache_of_buckets[LRU_queue_bucs[LRU_index]]);
    Cache_of_buckets[LRU_queue_bucs[LRU_index]] = NULL;
    Num_bucs_in_cache--;

    // determine the new least recently used cache bucket
    LRU_index = (LRU_index + 1) % Num_bucs_lev2;
  }

  return char_returned;  // returning char at absolute position k
}


/* ************************************************************
   Receives in input a bucket in the MTF form, having length 
   len_mtf; returns the original bucket where MTF-ranks have 
   been explicitely resolved. The characters obtained from Mtf[]
   are UNmapped according to the ones which actually occur into
   the superbucket. Therefore, the array Inv_map_b[] is necessary
   to unmap those chars from Alpha_size_b to Alpha_size_sb.
   ************************************************************ */
void unmtf_unmap(uchar *mtf_seq, int len_mtf)
{
  int i,j,rank;
  uchar next;

  /* ------ decode *inplace* mtf_seq -------------------- */
  for(j=0; j<len_mtf;j++ ) {
    rank = mtf_seq[j];
    assert(rank<Alpha_size_b); 
    next=Mtf[rank];              // decode mtf rank 
    for(i=rank;i>0;i--)          // update mtf list
      Mtf[i]=Mtf[i-1];
    Mtf[0]=next;                 // update mtf[0]
    mtf_seq[j]=Inv_map_b[next];  // apply invamp
  }
}


/* *********************************************************
   decode a unary code read from the proper stream. 
   The output is the # of zeroes we see before we see a 1 
   1 --> 0
   01 --> 1
   001 --> 2  
     etc.
   ********************************************************* */
static __inline__ int decode_unary(void)
{
  int bit_read(int);
  int i=0;

  while(bit_read(1)==0) i++;
  return i;
}


/* ***************************************************************
   open input files mapping them to the proper stream
   *************************************************************** */
void my_open_file(char *infile_name)
{
  FILE *my_fopen(const char *path, const char *mode);  

  /* ------ open input and output files ------ */
  if(infile_name==NULL){
    fprintf(stderr,"You should specify a file name -my_open_file-\n");
    exit(1);
  } else {
    Infile = my_fopen( infile_name, "rb"); // b is for binary: required in DOS
    if(Infile==NULL) {
      fprintf(stderr,"Unable to open file %s ",infile_name);
      perror("(my_open_file)");
      exit(1);
    }
  }
}
        


/* ***************************************************************
   Disable the cache system
   This function must be called when you do not want that any cache
   is used.
   *************************************************************** */
void disable_bwi_cache(void)
{
  Use_bwi_cache=0;
}

/* ***************************************************************
   Initialize the cache system
   This function must be called only once for each run and
   only when Infile is a bwi file
   NOTE: I have renamed this function (the old name was Init_cache_system())
   to make it clear that it can be used only for bwi files. GM 15-dic-00
   *************************************************************** */
void init_bwi_cache(void)
{
  void init_bit_buffer(void);
  void out_of_mem(char *);
  int uint_read(void),bit_read(int);
  int my_fseek(FILE *f, long offset, int whence);
  int i,text_size, buc_size_lev2;

  // compute number of level 2 buckets
  my_fseek(Infile,0,SEEK_SET);       // rewind
  init_bit_buffer();                 // initialize read-buffer
  my_fseek(Infile,1,SEEK_SET);
  text_size = uint_read();
  my_fseek(Infile,11,SEEK_SET);
  buc_size_lev2 = bit_read(8)<<10;
  Num_bucs_lev2 = (text_size + buc_size_lev2 - 1)/buc_size_lev2;

  // compute cache size
  Max_cached_buckets = (int) (Cache_percentage * Num_bucs_lev2);
  assert(Max_cached_buckets <= Num_bucs_lev2);
  if(Max_cached_buckets==0) {
    Use_bwi_cache=0;
    return;
  }
  
  // Allocate and initialize the data structures
  Use_bwi_cache = 1;
  Num_bucs_in_cache = 0;
  Cache_of_buckets = (uchar **) malloc(sizeof(uchar *) * Num_bucs_lev2);
  if(Cache_of_buckets == NULL)
    out_of_mem("init_cache_system");
  for(i=0; i<Num_bucs_lev2; i++)
    Cache_of_buckets[i] = NULL;

  // ---- Initialize the LRU_queue 
  LRU_queue_bucs = (int *) malloc(sizeof(int) * Num_bucs_lev2);
  if(LRU_queue_bucs == NULL)
    out_of_mem("Unable to alloc LRU_queue -Init_cache_system-");
  LRU_index = 0;
}


/* *******************************************************************
   this function can be called after the cache has been used
   in order to print a reporton the cache usage
   ******************************************************************* */
void report_bwi_cache_usage(void)
{
  int fsize;
  int csize = (Num_bucs_in_cache * Bucket_size_lev2);

  if(Use_bwi_cache) {
    fseek(Infile,0L,SEEK_END);
    fsize = ftell(Infile);
    fprintf(stderr, "Cache size = %d\n",csize); 
    fprintf(stderr, "Compressed file size = %d \n",fsize);
    fprintf(stderr, "Total occupied space = %d bytes\n", csize + fsize);
    fprintf(stderr, "Max # of cachable buckets = %d.\n", Max_cached_buckets);
    fprintf(stderr, "Max # of buckets in cache = %d.\n", Num_bucs_in_cache);
    fprintf(stderr, "# of buckets in file = %d.\n\n", Num_bucs_lev2);
  }
}


/* *******************************************************************
   this function reverses the passed string, ending with \0
   ******************************************************************* */
char *reverse_string(char *source)
{
  int i;
  int len = strlen(source);
  char *tmp = (char *) malloc(len+1);

  for(i=0; i<len; i++)
    tmp[i] = source[len-i-1];

  tmp[i]='\0';

  return tmp;
}


/* *******************************************************************
   These procedures are called by bwhuffw
   ******************************************************************* */




