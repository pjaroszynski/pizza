/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   random access inside .bwi files
   bwrand.c   
   P. Ferragina & G. Manzini, 25 Sept. 2000

   This program has been derived from bwsearch.c 
   The input is a row number r of the cyclic rows matrix and
   two integers n,m.
   The output are the last m chars of row r followed 
   by the first (n+m) chars of row r. In other words we find
   the text position correspondig to row r and we print 2m+n
   chars sourronding it. The idea is that this program is used to 
   display m chars preceedings and followings a pattern of length n
   which has been previously located using bwsearch. 
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "common.h"

#define BUFFER_SIZE 100000
#define EOF_shift(n) (n < s->bwt_eof_pos) ? n+1 :  n


/* --------- external variable used for the cache system ---------- */
double Cache_percentage;  // size of cache (% wrt uncompressed size)

int Report_position;    // if !=0 report the position of the input row 
int Url_id;           // rank of the url to be reported (counted from 1)



/* ***********************************************************
   read the text sourronding a given pattern. It is assumed that 
   the pattern starts at "row" (we do not know its position in 
   the input text), and we get the clen chars/words preceeding and 
   following it (according to type_count).
   ********************************************************** */
void get_stext(int row, int plen, int clen, int word_count)
{
  void output_char(char);
  int go_back(int, int, char *, bwi_out *s);
  int go_forw(int, int, char *, bwi_out *s);
  int go_back_word(int, int, char *, bwi_out *s, int bsize);
  int go_forw_word(int, int, char *, bwi_out *s, int bsize);
  char *textb;
  char *textf;
  int i, back, forw;
  bwi_out sp, *s;

  /* ---- get info: the bwi file has been already opened ---- */
  s = &sp;
  read_basic_prologue(s);
  /* --- report the row position if required ---- */
  if(Report_position) {
    if(s->skip==0) printf("-1\n");        // no locate information on file
    else printf("%d\n",get_occ_pos(s, row)); // 
  }
  /* ------- allocate memory for sourronding text ----- */
  textb=(char *) malloc(BUFFER_SIZE);
  textf=(char *) malloc(BUFFER_SIZE);
  if((textb==NULL) || (textf == NULL)) out_of_mem("get_stext");
  /* --- get clen chars preceding the current position --- */
  if (word_count) {
    back = go_back_word(row, clen, textb,s,BUFFER_SIZE); 
    assert(back<=BUFFER_SIZE);
  } else { 
    back = go_back(row, clen, textb,s); 
    assert(back<=clen);
  }
  /* --- get plen+clen chars from the current position --- */
  if (word_count) { 
    forw = go_forw_word(row, clen+1, textf,s,BUFFER_SIZE); 
    assert(forw<=BUFFER_SIZE);
  }
  else { 
    forw = go_forw(row, clen+plen, textf,s); 
    assert(forw<=clen+plen);
  }
  if(forw<plen) fatal_error("Error in pattern length -get_stext-");
  /* ---- print preceeding context ---- */
  if(back<clen)
    printf("[BOF]");
  for(i=0;i<back;i++) 
    output_char(textb[back-1-i]);  
  /* --- print pattern ----- */
  printf("<b>");
  for(i=0;i<plen;i++) 
    output_char(textf[i]);
  printf("</b>");
  /* --- print following chars ---- */
  for(i=0;i<forw;i++)
    output_char(textf[plen+i]);
  if(forw<clen)
    printf("[EOF]");
  printf("\n");
  /* ----- free mem and exit ---- */
  free(textb);
  free(textf);
  return;
}



/* ***********************************************************
   read the url corresponding to a url_id.
   ********************************************************** */
void get_url_byid(int Url_id, char **url)
{

  int get_url_back(int row, char *url_text, bwi_out *s, int bsize);
  int skip_url_back(int row, int num_url, bwi_out *s);
  void read_prologue(bwi_out *s);
  char *reverse_string(char *source);
  int pos_array,len,row,newrow,null_occ;
  int tobe_skipped,num_marked,last_marked;
  bwi_out sp, *s;

  /* ------- allocate memory for sourronding text ----- */
  *url=(char *) malloc(BUFFER_SIZE);
  if(*url==NULL) out_of_mem("alloc url space");


  /* ---- get info: the bwi file has been already opened ---- */
  s = &sp;
  read_prologue(s);

  null_occ=s->bwt_occ[1];  // # occurrences of NULL
  len = int_log2(null_occ); 
  num_marked = null_occ/s->skip + 1;
  last_marked = (null_occ/s->skip) * s->skip;

  if(Url_id > null_occ)
    fatal_error("Url id too large!\n\n");

  if (Url_id % s->skip)  // floor - 1
    pos_array = Url_id/s->skip;
  else 
    pos_array = Url_id/s->skip-1;

  // Get the row of the (jump * skip)-th url (counted from 0)
  my_fseek(Infile, ((pos_array*len)/8) + Start_prologue_occ,SEEK_SET);

  init_bit_buffer();

  if((pos_array*len) % 8)
    bit_read( (pos_array*len) % 8);  // skip bits not belonging to this pos  

  row=bit_read(len);  // read pos-representation

  if (Url_id <= last_marked) { 
    if (Url_id % s->skip)
      tobe_skipped = s->skip - (Url_id % s->skip);
    else
      tobe_skipped=0;
  } else {
      tobe_skipped = null_occ - Url_id;
  }


  newrow = skip_url_back(row, tobe_skipped, s);
  get_url_back(newrow, *url, s, BUFFER_SIZE);
  *url=reverse_string(*url);
  return;

}


/* ***********************************************************
   find the id for the given row.  ###PAOLO###
   ********************************************************** */
void get_id_byurl(char *url, int *Url_id)
{

  void read_prologue(bwi_out *s);
  int bwsearch(bwi_out *s, uchar *p, int len, int *sp, int *ep);
  int get_url_rank(bwi_out *s, int row);
  int go_back(int row, int len, char *dest, bwi_out *s);
  int sp,ep;
  bwi_out *s,ss;


  
  /* ---- get info: the bwi file has been already opened ---- */
  s = &ss;
  read_prologue(s);

  if (s->skip==0)
    fatal_error("s->skip is zero in get_id_byrow\n");

  // Searches for the url
  bwsearch(s,url,strlen(url),&sp,&ep);

  // Checks if more than one occ
  if (ep-sp>0)
    printf("There are more than 1 occurrences\n\n");

  // Compute the url rank
  *Url_id = get_url_rank(s,sp);

}


/* ********************************************************
   print a char to stdout taking care of newlines ,tab etc.
   To be completed!!!
   ******************************************************** */
void output_char(char c)
{      
  switch(c) 
    {
    case '\n': printf("[\\n]"); break;
    case '\r': printf("[\\r]"); break;
    case '>':  printf("&gt;"); break;
    case '<':  printf("&lt;"); break;
    default: printf("%c", c);
  }
}







