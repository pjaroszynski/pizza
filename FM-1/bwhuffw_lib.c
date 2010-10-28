/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   pattern search in files compressed with huffword+bwi
   bwhuffw.c   
   P. Ferragina & G. Manzini, 26 Sept 2000
   (ultima modifica 31 Gen 2001 [Stefano])

   based on bwsearch.c/bwxml.c
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "common.h"
#include <ctype.h>
#include <limits.h>

// --------------------------------------------------------------------
// macro for going from the first to the last column of the bwt matrix
// -------------------------------------------------------------------- 
#define EOF_shift(n) (n < s->bwt_eof_pos) ? n+1 :  n

/* --------- external variables used for the cache system ---------- */
extern double Cache_percentage;  // size of cache (% wrt uncompressed size)

/* --- external vars used by get_word_rank() and get_huffword_rank() --- */ 
extern uchar Inv_map_sb[256];  // inverse map for the current superbucket
extern int Alpha_size_sb;      // actual size of alphabet in superbucket

/* --------- "local" global variables ----------------- */
extern char *Dict_file_name;   // Name of file containing the huffword dictionary
extern char *Header_file_name; // Name of file containing the huffword header
extern char *Body_file_name;   // Name of the file containing huffword dictionary
extern int Word_search;        // type of word search substring, prexif, etc
extern int Report_occ;         // if !=0 report the postions of the occurrences   
extern int Locate_occ;         // if !=0 compute the postions of the occurrences   
extern int Yes_no_search;      // if !=0 only compute the # of matching dict. word
extern int Chk_uncompr;        // check the result of a query by looking at the 
                        // uncompressed files 

// data structure representing a single huffman code
// the code consists of the "len" LSB of codeword
typedef struct {
  uint32 code;          // codeword
  uchar len;            // codeword len (0<len <= 4) 
} huffcode;

// data strucutre containing the Header infos
extern int FirstCW[5],Offset[5], NumCW;



/* *********************************************************************
   This is the main procedure for the huffword search.
   We start by calling the procedure dict_search() which returns
   the list of huff codes corresponding to the words containing 
   the patterns (where the meaning of containing depends onthe global 
   variable Word_search which specify a substring/prefix/suffix/fullword
   search.
   Then for each huffman codes we search its occurrences in the body
   of the compressed file.
   ********************************************************************* */
void huffword_search(char *pattern, int pattern_from_file)
{
  huffcode *search_dict(char *pattern, int pat_len, int *num_codes);
  int get_huffword_rank(bwi_out *s, int row);
  huffcode *code_list, wcode;
  void chk_occ_list(huffcode wcode, int *occ_list, int num_occ);
  int pat_len, num_codes, byte;
  int i, j, pos, tot_numocc;
  int start_pos_occ,end_pos_occ,num_occ, *occ_list=NULL;
  double start,end, start0, end0, count_time, locate_time;
  bwi_out s_main;
  bwi_out *s = &s_main; 
  uchar hpattern[4];          // huffman code of a matching dictionary word

  /* ----- pattern is read from a file? ----- */
  if(pattern_from_file)
    pattern = read_pattern_from_file((char *) pattern, &pat_len);
  else 
    pat_len = strlen(pattern);

  /* ----- get list of huffman codes -------- */
  code_list = search_dict(pattern, pat_len, &num_codes);
 
  printf("Matching dictionary words: %d\n",num_codes);
  printf("Time for dictionary searching: %.4f secs.\n",getTime());
  if(Yes_no_search || (num_codes==0)) {
    printf("\n"); return;
  }

  /* ------------ bwsearch in the compressed body ------ */
  my_open_file(Body_file_name);   // open body
  init_bwi_cache();               // init cache
  locate_time = 0;                // init timer

  // ----- read prologue ------ 
  start0 = getTime();
  read_basic_prologue(s); 
  if(Locate_occ && s->skip==0) {
    fprintf(stderr,"The file does not contain information for ");
    fprintf(stderr,"locating the occurrences!\n");
    exit(1);
  }

  // ---- loop on all huffman codes -----------------------
  tot_numocc = 0;                 // # of pattern found 
  for(i=0;i<num_codes ;i++) {
    wcode = code_list[i];
    // --- copy pattern -------
    for(j=0;j<wcode.len;j++) {
      byte = (wcode.code>>(j*8))  & 0xff;      // get 8 bits 
      hpattern[wcode.len-j-1] = (uchar) byte;  // write them to hpattern
    }
    // ----- count occurrences --------
    num_occ = bwsearch(s, hpattern, wcode.len, &start_pos_occ, &end_pos_occ);
    tot_numocc += num_occ;
    if(Verbose>1) 
      fprintf(stderr,"Pattern %08x [%d] (length %d) occurs %d times\n",
                    wcode.code,wcode.code,wcode.len,num_occ);
	
    // ----- locate occurrences ----------
    if(Locate_occ && (num_occ > 0) ) {
      assert(s->skip!=0);
      if(Chk_uncompr) {
	// --- alloc array for positions
	occ_list = (int *) malloc(num_occ * sizeof(int));
	if(occ_list==NULL) out_of_mem("huffword_search (occ_list)");
      }
      start = getTime();
      for(j = start_pos_occ; j <= end_pos_occ; j++) {
	pos = get_huffword_rank(s, j);
	if(Chk_uncompr) occ_list[j-start_pos_occ]=pos;
	if(Report_occ)
	  fprintf(stderr,"pattern %x at position %d\n",wcode.code,pos);
      }
      end = getTime();
      locate_time += end-start;
      // test the correctness of occ_list
      if(Chk_uncompr) {
	chk_occ_list(wcode,occ_list,num_occ);
	free(occ_list);
      }
    } // if(Locate_occ >0 && ... )
  }   // for
  end0 = getTime();
  count_time = (end0-start0) -locate_time;

  printf("Total number of occurrences: %d\n",tot_numocc);
  printf("Time for counting: %.4f secs.  ",count_time);
  printf("Ave. time %.4f secs\n",count_time/num_codes);
  if(Locate_occ) {
    printf("Time for locating: %.4f secs.  ",locate_time);
    printf("Ave. time %.4f secs\n",locate_time/tot_numocc);
    printf("Tot. time (search+locate) %.4f secs.  ",end0-start0);
    printf("Ave. time per search %.4f secs\n",(end0-start0)/num_codes);
  }
  printf("Overall running time %.4f secs.\n\n", getTime());
  /* ---------  Report Cache Usage Information ------------- */
  if (Verbose) 
    report_bwi_cache_usage();
  /* ---------- close body and return -------- */
  my_fclose(Infile);
}


/* *****************************************************************
   this procedure searches all occurrences of pattern in the 
   dictionary stored in Dict_file_name. 
   It returns an array of huffman codes corresponding to the 
   dictionary words enclosing the given pattern string.
   The number of returned codes is stored in *num_codes
   ***************************************************************** */
huffcode *search_dict(char *pattern, int pat_len, int *num_codes)
{
  int bwsearch_word(bwi_out *, uchar *, int, int *, int *, int);
  int get_word_rank(bwi_out *, int);
  int int_cmp(const void *, const void *);
  void extract_hbz_header(char *path);
  void print_wordlist(int, int *, huffcode *);
  huffcode ComputeCW(int, int *, int*);
  int start_pos_occ,end_pos_occ,num_occ,rank_curr_occ,i,j;
  int *list_occ;
  huffcode *wcodes;
  double start,end;
  bwi_out s_main;
  bwi_out *s = &s_main; 

  my_open_file(Dict_file_name);   // set Infile = compressed dictionary
  disable_bwi_cache();            // disable bwi cache
  /* ---------  get all matching words occurrences ------------- */
  start = getTime();
  read_basic_prologue(s);
  num_occ = bwsearch_word(s, pattern, pat_len, &start_pos_occ, 
			  &end_pos_occ,Word_search);
  if(Verbose>1) {
    end = getTime();
    fprintf(stderr,"Found %d matching dictionary words in %.2f secs!\n",
            num_occ,(end-start));
  }
  if(Yes_no_search) {
    *num_codes = num_occ;   // we only care about the # of dict words
    my_fclose(Infile);      // close compressed dictionary
    return NULL;            // we do not care about their codes
  }
  if(num_occ==0) {
    *num_codes = 0;         // no dict words found
    my_fclose(Infile);      // close compressed dictionary
    return NULL;            // no codes to return
  }

  /* ---------  accumulate all dictionary ranks in an array ------------- */
  list_occ = (int *) malloc(sizeof(int) * num_occ);
  if(list_occ == NULL)
    out_of_mem("search_dict");
  for(i = 0; i < num_occ; i++)
    list_occ[i] = get_word_rank(s, start_pos_occ+i);
  my_fclose(Infile);              // close compressed dictionary

  /* ---------  sort the ranks ------------- */
  qsort(list_occ, num_occ, sizeof(int), int_cmp);

  /* --------- Clean the list of occurrences ----------- */
  if(Word_search==WSUBSTRING) {        // check for multiple in-word occs
    for(j=0, i=1; i < num_occ; i++) {
      if (list_occ[i] != list_occ[j])  // if list_occ[i] is not a duplicate ...
	list_occ[++j] = list_occ[i];   // ... save it
    } 
    assert(j < num_occ);
    num_occ = j+1;                     // number of cleaned occurrences
  }
  if(Verbose>1) {
    end = getTime();
    fprintf(stderr,"Found %d dictionary word-ranks in %.2f secs!\n",
                    num_occ,end-start);
  }

  /* ------------------------------------------------------------  
     Determine the codewords of the matching dictionary words 
     ------------------------------------------------------------ */  
  extract_hbz_header(Header_file_name);  // read hbz header
  wcodes = (huffcode *) malloc(sizeof(huffcode) * num_occ);
  if(wcodes == NULL)
    out_of_mem("search_dict (wcodes)");

  // For each rank we compute codeword and length and store them in wcodes[] 
  for(i=0; i<num_occ; i++) {
    rank_curr_occ =  list_occ[i];  
    wcodes[i] = ComputeCW(rank_curr_occ, Offset, FirstCW);
    assert((wcodes[i].len > 0) && (wcodes[i].len < 5));
  }

  /* --- print word list from the uncompressed dict (for debugging) --- */
  if(Chk_uncompr) 
    print_wordlist(num_occ, list_occ, wcodes);

  // ---------- return -------------
  *num_codes = num_occ;
  return wcodes;
}


/* *******************************************************************
   This is similar to the bwsearch procedure in search_main.c  
   The difference resides in the fact that we are operating on a dictionary
   and assume to have various types of pattern searches: word prefix, 
   word suffix, word enclosure and word. The first and last part of the 
   procedure modify the returned pattern to reflect these kind of searches 
   since dictionary words are separated by NULL.
   ******************************************************************* */
int bwsearch_word(bwi_out *s, uchar *p, int len, int *sp, int *ep, int wsearch)
{
  int i;
  uchar c;

  /* ---- remap pattern ------ */
  assert(len>0);
  for(i=0;i<len;i++) {
    if(s->bool_char_map[p[i]]==0) return 0;  /* char not in file */
    p[i]=s->char_map[p[i]];                  /* remap char */
  }
  
  /* ---- get initial sp and ep values ----- */ 
  if ((wsearch == WFULL) || (wsearch == WSUFFIX))   // suffix or word
    { i = len; c = s->char_map[0]; }
  else 
    { i=len-1; c=p[i];}

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

  if ((wsearch == WPREFIX) || (wsearch == WFULL))  // prefix or word
    { 
      c = s->char_map[0];
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


/* **********************************************************
   This is the procedure to retrieve the rank of the word which
   "contains" the position where the present row start. 
   Recall that we have marked (a fraction of) the rows which 
   start with the NULL char. The main loop consists in 
   checking if the current row is marked (that is, starts with NULL
   and is 0 modulo s->skip). It it is not, we appy the LF mapping
   to get the row corresponding to the previous char.
   ********************************************************** */
int get_word_rank(bwi_out *s, int row)
{
  void get_info_sb(bwi_out *,int,int *);
  uchar get_info_b(bwi_out *,uchar,int,int *,int);
  uchar c,c_sb;
  int len,occ_sb[256],occ_b[256],n,curr_row;
  int null_occ, counted_null, delta, null_start_row;
  
  assert(s->skip!=0);
  // ------ compute number of bits used for each rank 
  if(s->chosen_char==s->alpha_size-1)
    null_occ = s->text_size-s->bwt_occ[s->chosen_char];
  else
    null_occ = s->bwt_occ[s->chosen_char+1] - s->bwt_occ[s->chosen_char];
  assert(null_occ>0 && null_occ < s->text_size);
  len = int_log2(null_occ);

  // ----- compute start of marked rows and current row
  null_start_row = s->bwt_occ[s->chosen_char];
  curr_row = row;
  counted_null = 0;

  /* ---- search row starting with NULL and marked  ----- */ 
  while(1){
    delta = curr_row - null_start_row;
    if(delta>=0 && delta<null_occ) {
      // ---- curr_row starts with a NULL ----
      if((delta % s->skip)==0) {
	// the row is marked! retrieve its position
	n = delta/s->skip;     // Rescale
	my_fseek(Infile, ((n*len)/8) + Start_prologue_occ,SEEK_SET);
	init_bit_buffer();
	if((n*len) % 8)
	  bit_read( (n*len) % 8);  // skip bits not belonging to this pos  
	// return the rank + number of skipped NULL
	return(bit_read(len) + counted_null); 
      }
      else
	counted_null++; // not marked 
    }

    // ----- special case: beginning of file --------
    if (curr_row == s->bwt_eof_pos) // we know that s->text[0]=NULL hence
      return (counted_null-1);      // counted_null has been just incremntd 
 
    /* ------------------------------------------------------------- 
       now we search the row corresponding to the preceeding char
       ------------------------------------------------------------- */
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

    // compute the corresponding row
    curr_row = s->bwt_occ[c_sb] + n;

  }
  fatal_error("This code should not be reached (get_word_rank)\n");
  exit(1);
}

/* **********************************************************
   This is the procedure to retrieve the rank of the word in the text
   where the present row start. This procedure derives from
   get_word_rank() so look at the comment there. The only difference
   is that instead of considering the rows starting with NULL
   we consider the rows starting with a char>=128 (>=bit8_set 
   after the recoding). Recall that these are the chars marking
   the beginning of a huffword
   ********************************************************** */
int get_huffword_rank(bwi_out *s, int row)
{
  void get_info_sb(bwi_out *,int,int *);
  uchar get_info_b(bwi_out *,uchar,int,int *,int);
  uchar c,c_sb;
  int i,len,occ_sb[256],occ_b[256],n,curr_row;
  int counted_words,bit8_set,words_occ,delta,words_start_row;
  
  assert(s->skip!=0);
  /* ---- compute codeword of the first char>=128 ------- */
  for(i=0,bit8_set=0;i<128;i++) 
    if(s->bool_char_map[i]) bit8_set++;
  if(bit8_set==s->alpha_size) 
    fatal_error("Invalid body, no chars with the 8th bit set!\n");
  else if(bit8_set>s->alpha_size)
    fatal_error("Invalid s->bool_char_map[]!\n");
  // ------ compute number of bits used for each rank 
  words_occ = s->text_size - s->bwt_occ[bit8_set];
  assert(words_occ>0 && words_occ < s->text_size);
  len = int_log2(words_occ); 

  // ----- compute start of marked rows and current row
  words_start_row = s->bwt_occ[bit8_set];
  curr_row = row;
  counted_words = 0;

  /* ---- search row starting with NULL and marked  ----- */ 
  while(1){
    delta = curr_row - words_start_row;
    if(delta>=0 && delta<words_occ) {
      // ---- curr_row starts with a char >=128 ----
      if((delta % s->skip)==0) {
	// the row is marked! retrieve its position
	n = delta/s->skip;     // Rescale
	my_fseek(Infile, ((n*len)/8) + Start_prologue_occ,SEEK_SET);
	init_bit_buffer();
	if((n*len) % 8)
	  bit_read( (n*len) % 8);  // skip bits not belonging to this pos  
	// return the rank + number of skipped words
	return(bit_read(len) + counted_words); 
      }
      else
	counted_words++;     // not marked 
    }

    // ----- special case: beginning of file --------
    if (curr_row == s->bwt_eof_pos) // we know that s->text[0]>=128 hence
      return (counted_words-1);     // counted_words has been just incremntd 

    /* ------------------------------------------------------------- 
       now we search the row corresponding to the preceeding char
       ------------------------------------------------------------- */
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

    // compute the corresponding row
    curr_row = s->bwt_occ[c_sb] + n;

  }
  fatal_error("This code should not be reached (get_word_huffword)\n");
  exit(1);
}


/* ---- comparison function to sort tags in order of apperance in T --- */
int int_cmp(const void *a, const void *b)
{
  return (*((int *)a) -  *((int *)b)); 
}

/* **********************************************************
   This procedure extracts data from the header of an hbzipped file
   and return them to the caller. Fields FirstCW and Offset are
   arrays containing the info used for constructing the Codewords.
   *********************************************************** */
void extract_hbz_header(char *path) 
{
  FILE *Header_File;	
  int i,j;

  if ((Header_File = fopen(path,"rb")) == NULL){
    fprintf(stderr,"Error opening file %s ", path);
    perror("(extract_hbz_header)");
    exit(1);
  }

  fseek(Header_File,4L,SEEK_SET);  // Skip the 4 reserved bytes

  // First codeword
  for (i=1; i<=4; i++)
    for (FirstCW[i]=0,j=0; j<4; j++)
      FirstCW[i] = (FirstCW[i] << 7) | fgetc(Header_File);

  // Offset
  for (i=1; i<=4; i++)
    for (Offset[i]=0,j=0; j<4; j++)
      Offset[i] = (Offset[i] << 8) | fgetc(Header_File);
  
  // Dictionary size
  for (NumCW=0,i=0; i<4; i++)
    NumCW = (NumCW << 8) | fgetc(Header_File);
  
  fclose(Header_File);
}


/* **********************************************************
   ComputeCW gets in input the rank of a dictionary term and the arrays
   FirstCW and Offset, and returns the codeword for that term into
   a structure of type "huffcode".
   *********************************************************** */
huffcode ComputeCW(int rank, int *ofs, int *FirstCW) {

  int i,CW,CW_temp,LenCW, offset;
  huffcode Out;
              	
  Out.code = 0;
  Out.len = 0;

  // Determines the CW length
  for(i=1; (rank < ofs[i]) && (i < 5); i++) ;
  LenCW = i;
  if (LenCW == 5){
    fprintf(stderr,"Error in CWlen computation\n");
    exit(0);  }

  // Determines the CW without tagging (i.e. 128 fan-out)
  CW = FirstCW[LenCW];   
  offset =  rank - ofs[LenCW];
  CW = CW + offset;

  // Determines the CW with tagging (i.e. 8th bit tagged of a byte)
  CW_temp=0;
  for(i=LenCW-1; i>=0; i--)
    CW_temp = (CW_temp << 8) | ((CW >> (7*i)) & 0x7F); 
  CW = CW_temp | (0x1 << (8*LenCW - 1));

  Out.code = (uint32) CW;
  Out.len = (uchar) LenCW;
		 				
  return Out;
}


/* ****************************************************************
   this function takes as input a list of word ranks and look for the 
   corresponding words in the uncompressed dictionary.
   It then prints words, ranks, and huffman codes to stderr
   The ranks must be in increasing order without duplicates! 
   **************************************************************** */ 
void print_wordlist(int num_occ, int *rank_list, huffcode *wcodes)
{
  uchar *dict_start;
  char *name_dict_uncompr;
  FILE *f;
  int i,j,len;

  /* ------- create compressed word list for extra debugging ----- */
  f = fopen("wordlist.body","wb");
  if(f==NULL) {
      fprintf(stderr,"Unable to open the word-listfile\n");
      perror("(print_wordlist)");
  }
  else {
    for(i=0; i<num_occ; i++)
      for(j=wcodes[i].len-1; j>=0; j--)
	putc((uchar)(wcodes[i].code >> (8 * j)) & 0x000000ff,f);
    fclose(f);
  }

  fprintf(stderr,"------ Retrieved words and their infos ------\n");
  // try to open the uncompressed dictionary for retrieving the words
  name_dict_uncompr = (char *) malloc(strlen(Dict_file_name)-3);
  if(name_dict_uncompr==NULL) out_of_mem("print_wordlist");
  strncpy(name_dict_uncompr,Dict_file_name,strlen(Dict_file_name)-4);
  name_dict_uncompr[strlen(Dict_file_name)-4] = (char) 0;
    if ((f = fopen(name_dict_uncompr,"rb")) == NULL){
    fprintf(stderr,"Unable to open file %s\n", name_dict_uncompr);
    perror("(print_wordlist)");
    return;
  }

  // mmap dictionary file 
  fseek(f,0L,SEEK_END);  // compute file length
  len = ftell(f);
  fseek(f,0L,SEEK_SET);  //rewind
  dict_start = (uchar *) mmap(0,len,PROT_READ,MAP_SHARED,fileno(f),0);
  
  // Word retrieval
  for(i=0,j=0;i<num_occ;i++){    
    for( ;j < rank_list[i];j++) {
      // skip one word. note: the ranks must be in increasing order 
      for(dict_start++; *dict_start!=0; dict_start++)
         ;
    }
    fprintf(stderr,"Word %16s: rank %8d ",dict_start+1,rank_list[i]);
    fprintf(stderr,"code %08x [%10u] ",wcodes[i].code, wcodes[i].code);
    fprintf(stderr,"code-len %1d\n",wcodes[i].len);
  }
  fprintf(stderr,"---------------------------------------------\n");
  munmap(dict_start,len);
  fclose(f);
  free(name_dict_uncompr);
}


/* *************************************************************
   this procedures checks (by looking in the uncompressed body)
   that the occurrences found by huffword_search() are correct.
   wcode is the huffword code to be checked and occ_list[]
   contains the word-positions found by huffword_search()
   ************************************************************* */
void chk_occ_list(huffcode wcode, int *occ_list, int num_occ)
{
  int i,words,c,next_occ;
  uint32 word_code;
  char *name_body_uncompr;
  FILE *f;

  // create name of uncompressed body
  name_body_uncompr = (char *) malloc(strlen(Body_file_name)-3);
  if(name_body_uncompr==NULL) out_of_mem("chk_occ_list");
  strncpy(name_body_uncompr,Body_file_name,strlen(Body_file_name)-4);
  name_body_uncompr[strlen(Body_file_name)-4] = (char) 0;
  // open uncompressed body
  if ((f = fopen(name_body_uncompr,"rb")) == NULL){
    fprintf(stderr,"Unable to open file %s\n", name_body_uncompr);
    perror("(chk_occ_list)");
    return;
  }

  // sort positions in order of occurrence
  qsort(occ_list, num_occ, sizeof(int), int_cmp);


  // start reding from file
  c =  getc(f);
  if(c<128 || c==EOF) {
    fprintf(stderr,"Invalid body (chk_occ_list)\n");
    exit(1);
  }
  word_code= (uint32) c;
  next_occ=words=0;  
  for(i=0;  ;i++) {
    c = getc(f); 
    // --- check if we are at the end of a word ----
    if(!(c==EOF || c>127))
      word_code = (word_code<<8) + ((uint32) c);  // no end of word
    else {
      // ------ check if the current word is equal to wcode -----
      if(word_code==wcode.code) {
	if(next_occ>=num_occ) {
	  fprintf(stderr,"!!!!! Word %x in position %d ",word_code, i);
	  fprintf(stderr,"appeared with occ_list empty\n");
	  exit(1);
	} 
	else if(words<occ_list[next_occ]) {
	  fprintf(stderr,"!!!!! Word %x in position %d ",word_code, i);
	  fprintf(stderr,"(%dth word) not in occ_list\n",words);
	  exit(1);
	}
	else if(words> occ_list[next_occ]) {
	  fprintf(stderr,"!!!!! Word %x should be in position ",word_code);
	  fprintf(stderr,"%d of body, but is not there\n",occ_list[next_occ]);
	  exit(1);
	}
	else 
	  next_occ++;   // the current word was in occ_list
      }
      // ------------- start of a new word -----------
      if(c==EOF)
	break;
      else {
      word_code= (uint32) c;
      words++;
      }
    }   // else
  }     // for
  if(next_occ<num_occ) {
    fprintf(stderr,"!!!!! Only %d occurrences (out of %d) ",next_occ,num_occ);
    fprintf(stderr,"of word %x are in the body\n", wcode.code);
    exit(1);
  }
  else {
    fprintf(stderr,"%d occurrences of %x ", num_occ,wcode.code);
    fprintf(stderr,"has been sucessfully checked\n");
  }
  fclose(f);
}
      

    
