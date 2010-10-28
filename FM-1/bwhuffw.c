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
char *Dict_file_name;   // Name of file containing the huffword dictionary
char *Header_file_name; // Name of file containing the huffword header
char *Body_file_name;   // Name of the file containing huffword dictionary
int Word_search;        // type of word search substring, prexif, etc
int Report_occ;         // if !=0 report the postions of the occurrences   
int Locate_occ;         // if !=0 compute the postions of the occurrences   
int Yes_no_search;      // if !=0 only compute the # of matching dict. word
int Chk_uncompr;        // check the result of a query by looking at the 
                        // uncompressed files 

// data structure representing a single huffman code
// the code consists of the "len" LSB of codeword
typedef struct {
  uint32 code;          // codeword
  uchar len;            // codeword len (0<len <= 4) 
} huffcode;

// data strucutre containing the Header infos
int FirstCW[5],Offset[5], NumCW;


// -----------------------------------------------
// ------ main procedure for huffword search -----
// -----------------------------------------------
int main(int argc, char *argv[])
{
  int getopt(int argc, char * const *argv, const char *options);
  void huffword_search(char *, int);
  extern char *optarg;
  extern int optind, opterr, optopt;
  char *pattern;
  int c, num_opt, pattern_from_file;

  if(argc<5) {
    fprintf(stderr, "Usage:\n\t%s ",argv[0]);
    fprintf(stderr,"[-lpruvy][-x 1|2|3][-c perc][-w 0|1|2|3] ");
    fprintf(stderr,"pattern head dict body\n");
    fprintf(stderr,"Valid options are:\n");
    fprintf(stderr,"\t-c    percentage of cached buckets (default 0)\n");
    fprintf(stderr,"\t-l    locate occurrences\n");
    fprintf(stderr,"\t-p    *pattern* gives the pattern file name\n");
    fprintf(stderr,"\t-r    report row # and position of the occurrences\n");
    fprintf(stderr,"\t-u    check query looking in the uncompressed files\n"); 
    fprintf(stderr,"\t-v    produces a verbose output\n");
    fprintf(stderr,"\t-y    only compute # of matching dict. words\n");
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
  num_opt = opterr=0;
  pattern_from_file = 0;
  Type_mem_ops = EXT_MEM;
  Cache_percentage = 0;
  Report_occ = 0;
  Locate_occ = 0;
  Word_search = 0;
  Chk_uncompr = 0;
  Yes_no_search = 0;

  while ((c=getopt(argc, argv, "yurpvlc:w:x:")) != -1) {
    switch (c)
    {
      case 'v':
        Verbose++; break;
      case 'y':
        Yes_no_search=1; break;
      case 'u':
        Chk_uncompr=1; break;
      case 'r':
        Report_occ=1; break;
      case 'l':
        Locate_occ=1; break;
      case 'x':
        Type_mem_ops = atoi(optarg); break;
      case 'w':
        Word_search = atoi(optarg); break;
      case 'c':
        Cache_percentage = atof(optarg); break;
      case 'p':
        pattern_from_file = 1;
        break;
      case '?':
        fprintf(stderr,"Unknown option: %c -main-\n", optopt);
        exit(1);
    }
    num_opt++;
  }

  /* ------- read pattern and filename --------- */
  if((optind!=argc-4)) {
    fprintf(stderr,"You must supply a pattern and 3 filenames! -main-\n");
    exit(1);
  }
  else {
    pattern = (char *) argv[optind++];
    Header_file_name = (char *) argv[optind++];
    Dict_file_name = (char *) argv[optind++];
    Body_file_name = (char *) argv[optind];
  }

  if(strlen(pattern)<=0) {
    fprintf(stderr,"Invalid pattern -main-\n");
    exit(1);
  }

  if((Cache_percentage < 0) || (Cache_percentage > 1)) {
    fprintf(stderr,"Cache percentage must be in [0,1] -main-\n");
    exit(1);
  } 
  if(strlen(Header_file_name)<=0) {
    fprintf(stderr,"Invalid file name -main-\n");
    exit(1);
  }
  if (check_bwi_suffix(Dict_file_name) == 0){
     fprintf(stderr,"The dictionary file name must end with .bwi -main-\n");
     exit(1);
  }
  if (check_bwi_suffix(Body_file_name) == 0){
     fprintf(stderr,"The body file name must end with .bwi -main-\n");
     exit(1);
  }
  if(Report_occ) 
    Locate_occ = 1;  // report implies locate
  if(Yes_no_search && Locate_occ)
    fatal_error("Options -y and -r or -l are incompatible (main)\n");
  if((Word_search<0) || (Word_search>3)) 
    fatal_error("Invalid word search option -main-\n");
  if((Type_mem_ops<1) || (Type_mem_ops>3))
    fatal_error("Invalid data access option -main-\n");
 
  if(Verbose>1) {
    fprintf(stderr,"\n*****************************************************");
    fprintf(stderr,"\n             bwhuffw  Ver 1.0\n");
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
  if(Verbose>1) 
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
  huffword_search(pattern,pattern_from_file);
  return 0;
}


