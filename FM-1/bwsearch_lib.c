/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   pattern search in .bwi files 
   bwsearch.c   
   P. Ferragina & G. Manzini, 10 June 2000 

   The algorithm fetches in main memory (variable s) only a few
   data regarding the basic infos of the compressed file, stored into
   the header of the file prologue. The search procedure, bwsearch(),
   moves over the buckets and superbuckets according to the
   algorithmic scheme introduced in the paper "Opportunistic Data
   Structures with Applications" (Ferragina-Manzini). Only when needed
   the infos about the (super)bucket content is fetched and
   decompressed. This way we reduce the amount of used memory and
   achieve the bound O(p) for checking the pattern occurrence. We
   refer the reader to the paper above for further details on the
   listing of pattern-occurrence locations.
   The speciality of the above approach is that we keep explicitely
   the occurrences of a selected character. This way, when searching
   for the exact location of a pattern occurrence we backtrack on the
   text until such a character is found; then we count is rank among 
   its occurrences in the BWT and thus derive its "identity" that
   can be used to retrieve its exact location in T (using the previous
   list properly stored at the end of the file).
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "common.h"

#define use_my_getline 1     // set to 0 to use standard getline()

// --------------------------------------------------------------------
// macro for going from the first to the last column of the bwt matrix
// -------------------------------------------------------------------- 
#define EOF_shift(n) (n < s->bwt_eof_pos) ? n+1 :  n

static void output_char(char);

/* --------- external variable used for the cache system ---------- */
double Cache_percentage;  // size of cache (% wrt uncompressed size)

/* --------- external variables used by unbwi() ----- */ 
uchar Inv_map_sb[256];  // inverse map for the current superbucket
int Alpha_size_sb;      // actual size of alphabet in superbucket

/* ---- "local" global variables ------- */
int Report_occ;       // if !=0 report the postions of the occurrences   
int Locate_occ;       // if !=0 compute the postions of the occurrences   
int Display_occ;      // if !=0 display the text sourronding each occurrence
int Oneline_report;   // report time occ and startrow in a single line



/* ****************************************************************
   this procedures read one line at a time from file pfile_name
   and search it in the bwi file
   **************************************************************** */
void multi_search(char *pfile_name)
{
  void read_basic_prologue(bwi_out *);
  int my_getline(char *, int, FILE *);
  FILE *pfile;   // file containing the patterns
  char *pattern;
  int i, j, plen, max_aux, tot_located;
  int start_pos_occ,end_pos_occ,num_occ;
  double start,end, start0, end0, count_time, locate_time;
  bwi_out s_main;
  bwi_out *s = &s_main; 
  int max_plen=50;

  // open file containing patterns
  pfile = fopen(pfile_name,"r"); // this will be a text file
  if(pfile==NULL) {
    fprintf(stderr, "Error opening file %s ",pfile_name);
    perror("(multi_search)");
    exit(1);
  }
  // alloc memory for pattern
  pattern = (char *) malloc(max_plen);
  if(pattern==NULL) out_of_mem("multi_search");
  max_aux=max_plen;
  locate_time = 0;               // init timer
  tot_located = 0;               // # of pattern located 

  // ----- read prologue only once ------ 
  start0 = getTime();
  read_basic_prologue(s); 
  if(Locate_occ && s->skip==0) {
    fprintf(stderr,"The file does not contain information for ");
    fprintf(stderr,"locating the occurrences!\n");
    exit(1);
  }
  // ---- loop until there are pattern in file
  for(i=0; ;i++) {
#if !use_my_getline
    plen = getline(&pattern,&max_aux,pfile);
    if(plen== -1) break;  // end of file
    if(max_aux>max_plen || plen!=(int) strlen(pattern)) {
      fprintf(stderr, "Invalid pattern in row %d (multi_search)\n",i);
      exit(1);
    }
#else
    plen = my_getline(pattern,max_aux,pfile);
    if(plen== -1) break;  // end of file
    if(plen>=max_aux) {
      fprintf(stderr, "Invalid pattern in row %d (multi_search)\n",i);
      exit(1);
    }
#endif
    // get rid of newline at the end of patterm
    pattern[plen-1]=0;
    // ----- count occurrences
    num_occ = bwsearch(s, pattern, plen-1, &start_pos_occ, &end_pos_occ);
    if(Verbose>1) 
      fprintf(stderr,"Pattern \"%s\" occurs %d times\n",pattern,num_occ);
    // ----- locate occurrences
    if(Locate_occ && (num_occ > 0) ) {
      assert(s->skip!=0);
      start = getTime();
      for(j = start_pos_occ; j <= end_pos_occ; j++)
	get_occ_pos(s, j);
      end = getTime();
      locate_time += end-start;
      tot_located += num_occ;
    }
  }
  end0 = getTime();
  count_time = (end0-start0) -locate_time;
  if(i==0) {
    fprintf(stderr,"Nothing to search!");
    return;
  }
  fprintf(stderr,"%d patterns searched. %d occs located\n",i,tot_located);
  fprintf(stderr,"Tot. time (search+locate) %.2f secs.  ",end0-start0);
  fprintf(stderr,"Ave. time per search %.4f secs\n",(end0-start0)/i);
  fprintf(stderr,"Time for searching: %.2f secs.  ",count_time);
  fprintf(stderr,"Ave. time %.4f secs\n",count_time/i);
  fprintf(stderr,"Time for locating: %.2f secs.  ",locate_time);
  fprintf(stderr,"Ave. time %.4f secs\n\n",locate_time/tot_located);
  return;
}


/* *****************************************************************
   this procedure searches all occurrences of pattern in the 
   compressed file. the exact position of each occurrence is
   computed if Locate_occ!=0 and  they are actually 
   reported if Report_occ!=0.
   If pattern_from_file!=0 the pattern is interpreted as the name
   of the file containing the pattern
   ***************************************************************** */
void single_search(char *pattern, int pattern_from_file)
{
  void read_basic_prologue(bwi_out *);
  void display_stext(int row, int plen, int clen, bwi_out *s);
  int pat_len,start_pos_occ,end_pos_occ,*list_occ,num_occ,i;
  double start,end;
  bwi_out s_main;
  bwi_out *s = &s_main; 

  /* ----- pattern is read from a file? ----- */
  if(pattern_from_file)
    pattern = read_pattern_from_file((char *) pattern, &pat_len);
  else 
    pat_len = strlen(pattern);

  /* ---------  count pattern occurrences ------------- */
  start = getTime();
  read_basic_prologue(s); 
  num_occ = bwsearch(s, pattern, pat_len, &start_pos_occ, &end_pos_occ);
  end = getTime();
  /* ------ one line outout (to be read by perl scripts) ------ */ 
  if(Oneline_report) {
    printf("%.2f %d %d\n",end-start,num_occ,start_pos_occ);
    return;
  }
  fprintf(stderr,"Found %d occ's in %.2f seconds!\n",num_occ,end-start);
  if(Verbose>0) 
    fprintf(stderr,"BWT rows from %d to %d\n",start_pos_occ,end_pos_occ); 
  /* ---------  list pattern occurrences ------------- */
  if (Locate_occ && num_occ > 0 ) {
    if(s->skip==0) {
      fprintf(stderr,"The file does not contain information for ");
      fprintf(stderr,"locating the occurrences!\n");
      exit(1);
    }
    list_occ = (int *) malloc(sizeof(int) * num_occ);
    if(list_occ == NULL)
      out_of_mem("single_search (list_occ allocation failed!)");

    start = getTime();
    // Accumulate occurrences in an array
    for(i = start_pos_occ; i <= end_pos_occ; i++)
      list_occ[i-start_pos_occ] = get_occ_pos(s, i);
    end = getTime();

    if (Report_occ) 
      for(i = 0; i < num_occ; i++) 
	fprintf(stderr,"Occ. %4d is at position %8d\n",i+1,list_occ[i]);

    fprintf(stderr,"Located %d occ's in %.2f seconds!\n",
	    num_occ,end-start);
  }
  else if(Display_occ) {
    start = getTime();
    for(i = 0; i < num_occ; i++)
      display_stext(i+start_pos_occ,pat_len,Display_occ,s);
    end = getTime();
    printf("Displayed %d occ's in %.2f seconds!\n",
	    num_occ,end-start);
  }
  fprintf(stderr,"\n");
}



/* ***********************************************************
   display the text sourronding a given pattern. Is is assumed that 
   the pattern starts at "row" (we do not know its position in 
   the input text), and we get the clen chars preceeding and 
   following it.
   ********************************************************** */
void display_stext(int row, int plen, int clen, bwi_out *s)
{
  int go_back(int, int, char *, bwi_out *s);
  int go_forw(int, int, char *, bwi_out *s);
  char *text;
  int i, back, forw;

  /* ------- allocate memory for sourronding text ----- */
  text=(char *) malloc(2*clen+plen);
  if(text==NULL) out_of_mem("get_stext");
  /* --- get clen chars preceding the current position --- */
  back = go_back(row, clen, text,s);
  assert(back<=clen);
  /* --- get plen+clen chars from the current position --- */
  forw = go_forw(row, clen+plen, text+clen, s);
  assert(forw<=clen+plen);
  if(forw<plen) fatal_error("Error in pattern length -display_stext-");
  /* ---- print preceeding context ---- */
  if(back<clen)
    printf("[BOF]");
  for(i=0;i<back;i++) 
    output_char(text[back-1-i]);  
  /* --- print pattern ----- */
  // printf("<b>");
  for(i=0;i<plen;i++) 
    output_char(text[clen+i]);
  // printf("</b>");
  /* --- print following chars ---- */
  for(i=0;i<clen;i++)
    output_char(text[plen+clen+i]);
  if(forw<clen)
    printf("[EOF]");
  printf("\n");
  return;
}


/* ********************************************************
   print a char to stdout taking care of newlines ,tab etc.
   To be completed!!!
   ******************************************************** */
static void output_char(char c)
{      
  switch(c) 
    {
    case '\n': printf("[\\n]"); break;
    case '\r': printf("[\\r]"); break;
    default: printf("%c", c);
  }
}





/* ********************************************************************
   decompress a .bwi file by a repeated application of the lf mapping
   ******************************************************************* */
void unbwi(char *infile_name)
{
  void read_basic_prologue(bwi_out *);
  int occ(bwi_out *,int,uchar);
  void get_info_sb(bwi_out *,int,int *);
  uchar get_info_b(bwi_out *,uchar,int,int *,int);
  uchar c,c_sb;
  char *outfile_name;
  int i,occ_sb[256],occ_b[256],n,curr_row;
  double start,end;
  FILE *outfile;
  bwi_out s_main;
  bwi_out *s = &s_main;

  start = getTime();
  /* ---- open output file ------ */
  outfile_name = (char *) malloc(strlen(infile_name)+3);
  outfile_name = strcpy(outfile_name, infile_name);
  outfile_name = strcat(outfile_name, ".y");
  outfile = fopen(outfile_name,"wb");
  if(outfile==NULL) {
    fprintf(stderr,"Unable to open file %s ",outfile_name);
    perror("(unbwi)");
    exit(1);
  }

  /* ---- uncompress the file ------ */
  read_basic_prologue(s);
  curr_row=0; 
  i=s->text_size;
  while(1) {
    // fetches info from the header of the superbucket
    get_info_sb(s,curr_row,occ_sb);  
    // fetches occ into occ_b properly remapped and returns
    // the remapped code for occ_b of the char in the  specified position
    c = get_info_b(s,NULL_CHAR,curr_row,occ_b,WHAT_CHAR_IS);  
    assert(c < Alpha_size_sb);
    c_sb = Inv_map_sb[c];
    assert(c_sb < s->alpha_size);
    fseek(outfile,--i,SEEK_SET);           // write char
    putc(s->inv_char_map[c_sb],outfile);    
    if(Verbose>2 && ((i%65536) == 0))
      fprintf(stderr,".");
    n = occ_sb[c_sb] + occ_b[c] - 1;       // # of occ before curr_row
    curr_row = s->bwt_occ[c_sb] + n;       // get next row
    if(curr_row==s->bwt_eof_pos) break;    
    curr_row = EOF_shift(curr_row);        
  }
  if(i!=0)
    fatal_error("Decompression failed! (unbwi)");
  if(fclose(outfile)==EOF) {
    fprintf(stderr,"Unable to close file %s ",outfile_name);
    perror("(unbwi)");
    exit(1);
  }
  end = getTime();
  fprintf(stderr,"Decompression time: %.2f seconds\n\n",end-start);
}









