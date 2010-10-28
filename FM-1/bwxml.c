/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   pattern search in .bwi files
   bwxml.c   
   P. Ferragina & G. Manzini, 26 Sept 2000

   based on bwsearch.c
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "common.h"
#include <ctype.h>

/* --------- external variable used for the cache system ---------- */
extern double Cache_percentage;  // size of cache (% wrt uncompressed size)

/* --- structure containing position and bwt row of a pattern occurrence --- */
typedef struct {
  int pos;     // position of the occ in the input text;
  int row;     // position of the occ in the bwt matrix
} occ_data;

/* ---- "local" global variables ---- */
int Report_rows;      // if !=0 report the row numbers of the occurrences
int Word_search;      // 0=substring, 1=prefix, 2=suffix, 3=full word
int Count_only;      // 0=full search, 1=count only

// -----------------------------------------------
// ----------- main procedure for xml search -----
// -----------------------------------------------
int main(int argc, char *argv[])
{
  int getopt(int argc, char * const *argv, const char *options);
  void unbwi(char *);
  void xml_search(char *, int, char *);
  extern char *optarg;
  extern int optind, opterr;
  char *infile_name, *pattern, *tag;
  int c, num_opt, pattern_from_file;

  if(argc<3) {
    fprintf(stderr, "Usage:\n\t%s ",argv[0]);
    fprintf(stderr,"[-prv] [-x 1|2|3] [-c perc] [-w 0|1|2|3]");
    fprintf(stderr,"[-t tag] pattern bwifile\n");
    fprintf(stderr,"Valid options are:\n");
    fprintf(stderr,"\t-p    *pattern* gives the pattern file name\n");    
    fprintf(stderr,"\t-r    report row # and position of the occurrences\n");
    fprintf(stderr,"\t-v    produces a verbose output\n");
    fprintf(stderr,"\t-C    only count occurrences\n");
    fprintf(stderr,"\t-c    percentage of cached buckets (default 0)\n"); 
    fprintf(stderr,"\t-t    enclosing tag (default none)\n");
    fprintf(stderr,"\t-w    word search: 0=Subs 1=Pfx 2=Sfx ");
    fprintf(stderr,"3=FullWord (default 0)\n");
    fprintf(stderr,"\t-x    data access: 1=Ext, 2=Mmap, 3=In (default 1)\n\n");
    fprintf(stderr,"The command line was: ");
    for(c=0;c<argc;c++)
      fprintf(stderr,"%s ",argv[c]);
    fprintf(stderr,"\n\n");
    exit(1);
  }

  /* ----------------- read options --------------------- */
  Verbose=0;
  infile_name=NULL;
  num_opt = opterr=0;
  pattern_from_file = 0;
  Type_mem_ops = EXT_MEM;
  Cache_percentage = 0;
  Report_rows = 0;
  Word_search = 0;
  Count_only = 0;
  tag = NULL;

  while ((c=getopt(argc, argv, "rpvCc:t:w:x:")) != -1) {
    switch (c)
    {
      case 'v':
        Verbose++; break;
      case 'r':
        Report_rows=1; break;
      case 'x':
        Type_mem_ops = atoi(optarg); break;
      case 'w':
        Word_search = atoi(optarg); break;
      case 'C':
        Count_only = 1; break;
      case 'c':
        Cache_percentage = atof(optarg); break;
      case 'p':
        pattern_from_file = 1; break;
      case 't':
        tag = optarg; break;
      case '?':
        fprintf(stderr,"Unknown option: %c -main-\n", c);
        exit(1);
    }
    num_opt++;
  }

  /* ------- read pattern and filename --------- */
  if(optind!=argc-2) {
    fprintf(stderr,"You must supply a pattern and a filename! -main-\n");
    exit(1);
  }
  else {
    pattern = (char *) argv[optind++];
    infile_name = (char *) argv[optind];
  }
  if(strlen(pattern)<=0) {
    fprintf(stderr,"Invalid pattern -main-\n");
    exit(1);
  }

  if((Cache_percentage < 0) || (Cache_percentage > 1)) {
    fprintf(stderr,"Cache percentage must be in [0,1] -main-\n");
    exit(1);
  } 
  if(strlen(infile_name)<=0) {
    fprintf(stderr,"Invalid file name -main-\n");
    exit(1);
  }
  if (check_bwi_suffix(infile_name) == 0){
     fprintf(stderr,"The file name must end with .bwi -main-\n");  
     exit(1);
  }
  if((Word_search<0) || (Word_search>3)) 
    fatal_error("Invalid word search option -main-\n");
  if((Type_mem_ops<1) || (Type_mem_ops>3))
    fatal_error("Invalid data access option -main-\n");

 
  if(Verbose>1) {
    fprintf(stderr,"\n*****************************************************");
    fprintf(stderr,"\n             bwxml  Ver 1.0\n");
    fprintf(stderr,"Created on %s at %s from %s\n",__DATE__,__TIME__,__FILE__);
    fprintf(stderr,"*****************************************************\n");
  }
  if(Verbose) {
    fprintf(stderr,"Command line: ");
    for(c=0;c<argc;c++)
      fprintf(stderr,"%s ",argv[c]);
    fprintf(stderr,"\n");
  }

  /* ----- Hilights the type of memory management ------ */
  if(Verbose) 
    switch(Type_mem_ops)
      {
      case EXT_MEM: 
	fprintf(stderr,"Memory management: via fopen, fread, getc, ...\n"); 
	break;
      case EXT_MMAP: 
	fprintf(stderr,"Memory management: via mmap() system call ...\n");
	break;
      case IN_MEM: 
	fprintf(stderr,"Memory management: internal memory ...\n"); 
	break;
      default: 
	fprintf(stderr,"Error: -- type mem ops --\n"); exit(1); break;
      }

  /* ---------  open input file ------------- */
  my_open_file(infile_name);

  /* ---------  Initialize the Cache System ------------- */
  init_bwi_cache();

  xml_search(pattern,pattern_from_file,tag);

  /* ---------  Report Cache Usage Information ------------- */
  if (Verbose) 
    report_bwi_cache_usage();

  my_fclose(Infile);
  return 0;
}



/* *****************************************************************
   if tag==NULL this procedure returns all occurrences of pattern in the 
   compressed file. Otherwise only the occurrences enclosed between
   <\tag and tag> are returned.
   If pattern_from_file!=0 the pattern is interpreted as the name
   of the file containing the pattern
   ***************************************************************** */
void xml_search(char *pattern, int pattern_from_file, char *tag)
{
  void read_basic_prologue(bwi_out *);
  int occ_cmp(const void *a, const void *b);
  int get_prefix_occ(occ_data *o, int num, bwi_out *s);
  int get_suffix_occ(occ_data *o, int num, int plen, bwi_out *s);
  int get_enclosed_occ(occ_data *occ_list, int num, char *tag, bwi_out *s);
  int pat_len,start_pos_occ,end_pos_occ,num_occ,i;
  occ_data *occ_list = NULL;
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
  if(s->skip==0) {
    fprintf(stderr,"The file does not contain information for ");
    fprintf(stderr,"locating the occurrences -xml_search-!\n");
    exit(1);
  }
  num_occ = bwsearch(s, pattern, pat_len, &start_pos_occ, &end_pos_occ);

  if(num_occ>0) {
    occ_list = (occ_data *) malloc(num_occ * sizeof(occ_data));
    if(occ_list==NULL) out_of_mem("xml_search");
    /* ------- store row numbers of occurrences ------ */
    for(i = 0; i < num_occ; i++) {
      occ_list[i].row = i+start_pos_occ;
    }  
    /* ----- consider only prefixes, suffixes, full words, etc. ----- */
    switch(Word_search) 
      {
      case 0: 
	break;  // do nothing any substring is OK
      case 1: 
	num_occ = get_prefix_occ(occ_list, num_occ,s); 
	break;        // consider only substrings which are word prefixes
      case 2: 
	num_occ = get_suffix_occ(occ_list, num_occ,pat_len,s); 
	break;        // consider only substrings which are word suffixes
      case 3: 
	num_occ = get_prefix_occ(occ_list, num_occ,s); 
	num_occ = get_suffix_occ(occ_list, num_occ,pat_len,s); 
	break;        // consider only substrings which are full words	
      default:
	fatal_error("Invalid word search option -xml_search-\n");
      }

    if (!Count_only) {
      /* --- compute text positions of the remaining occurrences --- */
      for(i = 0; i < num_occ; i++) {
	occ_list[i].pos = get_occ_pos(s, occ_list[i].row);
	if(Verbose>2)
	  fprintf(stderr,"pos: %8d row: %8d\n", occ_list[i].pos,
                          occ_list[i].row);
      }
      /* ----- sort occ_list in order of occurrence in the text --- */
      if(num_occ>0) 
	qsort(occ_list, num_occ, sizeof(occ_data), occ_cmp);
      if(tag!=NULL && num_occ>0) 
	num_occ = get_enclosed_occ(occ_list, num_occ, tag, s);
    }
  }
  end = getTime();

  // ---  print the results of the query --------
  printf("%d\n%.2f\n",num_occ,end-start); // report # occ and elapsed time   
  if (Report_rows && !Count_only)
    for(i = 0; i < num_occ; i++)
      printf("%d %d\n",occ_list[i].pos,occ_list[i].row); // position and row #
  free(occ_list);
}

/* ***************************************************************
   discard, from the num occurrences which are in occ_list, those
   which are not enclosed in a pair of opening/closing tag
   return the number of survived occurrences
   *************************************************************** */
int get_enclosed_occ(occ_data *occ_list, int num_occ, char *tag, bwi_out *s)
{
  int *get_tag_list(char *, bwi_out *, int *);
  int i, t, status, occ_accepted, tags_seen, num_tags, *tag_list; 

  if(Verbose>1)
    fprintf(stderr,"%d occurrences to be filtered\n",num_occ);

  /* ---- the tag list has a very special format (see below) ---- */
  tag_list = get_tag_list(tag,s,&num_tags);
  if(num_tags==0)
    num_occ = 0;
  else {
    assert(tag_list!=NULL);
    /* ---- find which occ are enclosed between two tags --- */
    occ_accepted = tags_seen = 0;
    status = 0;
    for(i=0;i<num_occ;i++) {
      for( ;tags_seen < num_tags; ) {
	t=tag_list[tags_seen];
	if((t>>1) > occ_list[i].pos) break;
	else {
	  status = t & 1;   //get current status in the lsb
	  tags_seen++;
	}
      }
      if(status) 
	occ_list[occ_accepted++] = occ_list[i];
    }
    free(tag_list);
    num_occ = occ_accepted;  // update number of words
  }     // num_tags==0
  return num_occ;
}    


/* **********************************************************************
   compute the tag list as follows:
   1) find all the occurrences of an opening tag "<tag"
   2) find all the occurrences of a closing tag  "/tag>"
   3) create a list containing the location of both occurrences 
      shifted by one bit. in this extra bit (the status bit) we write 1 
      for an opening tag and 0 for a closing tag
   4) sort the list 
   5) scan the sorted list and remove nested occurrences 
   6) return the list 
   ********************************************************************** */
int *get_tag_list(char *tag, bwi_out *s, int *size)
{
  int tlen,i, t, level, *tag_list;
  int num_op, num_cl, start_cl,end_cl, start_op1, end_op1, start_op2, end_op2;
  int int_cmp(const void*, const void *);
  uchar *p;

  tlen = strlen(tag);
  assert(tlen>0);
  p = (char *) malloc(tlen+2);
  if(p==NULL) out_of_mem("get_tag_list (1)");
  /* ---- search opening tag without attributes---- */
  p[0]='<';
  for(i=0;i<tlen;i++) p[i+1] = tag[i];
  p[tlen+1]='>';
  num_op = bwsearch(s, p, tlen+2, &start_op1, &end_op1);
  /* ---- search opening tag with attributes---- */
  p[0]='<';
  for(i=0;i<tlen;i++) p[i+1] = tag[i];
  p[tlen+1]=' ';
  num_op += bwsearch(s, p, tlen+2, &start_op2, &end_op2);
  /* ---- search closing tag ---- */
  p[0]='/';
  for(i=0;i<tlen;i++) p[i+1] = tag[i];
  p[tlen+1]='>';
  num_cl = bwsearch(s, p, tlen+2, &start_cl, &end_cl);
  free(p);
  /* ---- do some checking ---- */
  if(Verbose>1) {
    fprintf(stderr,"%d opening tags ",num_op);
    fprintf(stderr,"(%d with attributes)\n", 1+ end_op2-start_op2); 
    fprintf(stderr,"%d closing tags\n",num_cl);
  }
  if(num_op!=num_cl) {
    fprintf(stderr,"The number of opening and closing tags does not match\n");
    fprintf(stderr,"(%d vs %d).\n",num_op,num_cl);
    exit(2);;
  }    
  if(num_op==0) { 
    *size=0; return NULL; // no tags
  }
  /* ---- we have some work to do ---- */
  tag_list = (int *) malloc(2*num_op*sizeof(int));
  if(tag_list==NULL) 
    out_of_mem("get_tag_list (2)");
  /* --- get location of opening tags  ---- */
  for(i=start_op1;i<=end_op1;i++) {
     t = get_occ_pos(s, i); 
     tag_list[i-start_op1] = 1+ (t<<1);  // staus bit ==1
  }
  for(i=start_op2;i<=end_op2;i++) {
     t = get_occ_pos(s, i); 
     tag_list[num_op-1-(i-start_op2)] = 1+ (t<<1);  // staus bit ==1
  }
  /* --- get location of closing tags  ---- */
  for(i=start_cl;i<=end_cl;i++) {
     t = get_occ_pos(s, i); 
     tag_list[num_op+i-start_cl] = (t<<1);  // staus bit ==0
  }
  /* -------- sort tag positions ------------ */
  qsort(tag_list, 2*num_op,sizeof(int),int_cmp);
  /* ---- remove nested tags ---- */
  *size = level=0;
  for(i=0;i<2*num_op;i++) {
    if(tag_list[i]&1) {    // opening tag
      if(level++==0) 
	tag_list[(*size)++]=tag_list[i]; // copy this tag
    }
    else {                // closing tag 
      if(--level<0) {
	fprintf(stderr,"Tags are not properly nested -get_tag_list-\n");
	exit(2);
      }
      if(level==0)
	tag_list[(*size)++]=tag_list[i]; // copy this tag
    }
  }
  if(*size%2) {       // the number of tags must be even 
    fprintf(stderr,"Tags are not properly nested -get_tag_list-\n");
    exit(2);
  }
  if(Verbose>1)
    fprintf(stderr,"There are %d pairs of non-nested tags\n",*size/2);
  return tag_list;
}


/* *****************************************************************
   remove from the list of occurrences those which are not 
   word prefixes. This is done by looking at the 
   character preceeding the occurrences.
   ***************************************************************** */
int get_prefix_occ(occ_data *o, int num, bwi_out *s)
{
  int go_back(int row, int len, char *dest, bwi_out *s);
  int i, n, written=0;
  char b[1];

  for(i=0;i<num;i++) {
    n=go_back(o[i].row,1,b,s);
    if( ! ((n>0)&&isalpha(b[0])) )
      o[written++] = o[i];
  }
  return written;
}

/* *****************************************************************
   remove from the list of occurrences those which are not 
   word suffixes. This is done by looking at the character 
   following each pattern occurrence. Note that the right way to search
   for suffix occurrences would be to selected them directly using 
   a modified bwsearch procedure.
   ***************************************************************** */
int get_suffix_occ(occ_data *o, int num, int plen, bwi_out *s)
{
  int go_forw(int row, int len, char *dest, bwi_out *s);
  int i, n, written=0;
  char *b;

  b= (char *) malloc(plen+1);
  for(i=0;i<num;i++) {
    n=go_forw(o[i].row,plen+1,b,s);
    assert(n==plen || n==plen+1);
    if( ! ((n>plen)&& (isalpha(b[plen]) || (b[plen] == '&'))))
      o[written++] = o[i];
  }
  return written;
}



/* --- comparison function used to sort occ in order of appeance in T --- */
int occ_cmp(const void *a, const void *b)
{
  return ((occ_data *) a)->pos - ((occ_data *) b)->pos; 
}
/* ---- comparison function used to sort tags in order of apperance in T --- */
int int_cmp(const void *a, const void *b)
{
  return *((int *) a) -  *((int *) b); 
}









