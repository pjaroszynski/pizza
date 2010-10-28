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


/* --------- external variable used for the cache system ---------- */
extern double Cache_percentage;  // size of cache (% wrt uncompressed size)

extern int Report_position;    // if !=0 report the position of the input row 
extern int Url_id;           // rank of the url to be reported (counted from 1)

// -----------------------------------------------
// ----------- Main procedure -------------
// -----------------------------------------------

int main(int argc, char *argv[])
{
  int getopt(int argc, char * const *argv, const char *options);
  void get_stext(int, int, int, int);
  void multi_search(char *);
  void single_search(char *, int);
  void get_url_byid(int Url_id, char **url);
  extern char *optarg;
  extern int optind, opterr;
  char *infile_name;
  int c, num_opt, row, plen, clen;
  int word_count;
  char *url_txt;

  if(argc<4) {
    fprintf(stderr, "Usage:\n\t%s ",argv[0]);
    fprintf(stderr,"[-v] [-x 1|2|3] [-c perc] [-w] [-u url_id]");
    fprintf(stderr,"row p_len c_len bwifile\n\n");
    fprintf(stderr,"to print the last c_len chars and the first ");
    fprintf(stderr,"p_len+c_len chars of row.\n");
    fprintf(stderr,"or to get the url_id-th URL in a file separated by NULL\n\n\n");
    fprintf(stderr,"Valid options are:\n");
    fprintf(stderr,"\t-v    produces a verbose output\n");
    fprintf(stderr,"\t-p    report text position of row\n"); 
    fprintf(stderr,"\t-c    percentage of cached buckets (default 0)\n"); 
    fprintf(stderr,"\t-w    c_len is expressed in words\n"); 
    fprintf(stderr,"\t-u    url_id counted from 1\n");
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
  Type_mem_ops = EXT_MEM;
  Cache_percentage = 0;
  Report_position=0;
  word_count = 0;
  Url_id = 0;

  while ((c=getopt(argc, argv, "pvwc:x:u:")) != -1) {
    switch (c)
    {
      case 'v':
        Verbose++; break;
      case 'w':
        word_count = 1; break;
      case 'p':
        Report_position = 1; break;
      case 'x':
        Type_mem_ops = atoi(optarg); break;
      case 'c':
        Cache_percentage = atof(optarg); break;
      case 'u':
        Url_id = atoi(optarg); break;
      case '?':
        fprintf(stderr,"Unknown option: %c -main-\n", c);
        exit(1);
    }
    num_opt++;
  }

  /* ------- read row plen clen and filename --------- */
  if(Url_id == 0)
    {
      if(optind!=argc-4) {
	fprintf(stderr,"You must supply four arguments -main-\n");
	exit(1);
      }
      row = atoi(argv[optind++]);
      plen = atoi(argv[optind++]);
      clen = atoi(argv[optind++]);

      if(row<0) 
	fatal_error("Invalid row number -main-0");
      if(plen<0)
	fatal_error("Invalid pattern_len -main-0");
      if(clen<0)
	fatal_error("Invalid context_len -main-0");
    } else 
      {row=clen=plen=0; }

  infile_name = (char *) argv[optind];

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
 
  if(Verbose>1) {
    fprintf(stderr,"\n*****************************************************");
    fprintf(stderr,"\n             bwrand  Ver 1.0\n");
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

  /* --------- do the work ------ */
  if(Url_id > 0){
    get_url_byid(Url_id, &url_txt);
    printf("\n%s\n",url_txt);
  } else
    get_stext(row, plen, clen, word_count);

  /* ---------  Report Cache Usage Information ------------- */
  if (Verbose) 
    report_bwi_cache_usage();

  my_fclose(Infile);
  return 0;
}
