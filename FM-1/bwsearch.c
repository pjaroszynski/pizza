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


/* --------- external variable used for the cache system ---------- */
extern double Cache_percentage;  // size of cache (% wrt uncompressed size)

/* --------- external variables used by unbwi() ----- */ 
extern uchar Inv_map_sb[256];  // inverse map for the current superbucket
extern int Alpha_size_sb;      // actual size of alphabet in superbucket

/* ---- "local" global variables ------- */
extern int Report_occ;       // if !=0 report the postions of the occurrences   
extern int Locate_occ;       // if !=0 compute the postions of the occurrences   
extern int Display_occ;      // if !=0 display the text sourronding each occurrence
extern int Oneline_report;   // report time occ and startrow in a single line


// ---------------------------------------------------
// ----------- Main procedure for the bwsearch tool
// ---------------------------------------------------
int main(int argc, char *argv[])
{
  int getopt(int argc, char * const *argv, const char *options);
  void unbwi(char *);
  void multi_search(char *);
  void single_search(char *, int);
  extern char *optarg;
  extern int optind, opterr, optopt;
  char *infile_name, *pattern;
  int c, num_opt, pattern_from_file, multiple_search,decompress;

  if(argc<3) {
    fprintf(stderr, "Usage:\n\t%s ",argv[0]);
    fprintf(stderr,"-d [-x 1|2|3] [-c perc] bwifile \n");
    fprintf(stderr,"to decompress a file, or\n\t%s ",argv[0]);
    fprintf(stderr,"[-lmoprv][-x 1|2|3][-s len][-c perc] pattern bwifile\n");
    fprintf(stderr,"to search a pattern. ");
    fprintf(stderr,"Valid options are:\n");
    fprintf(stderr,"\t-d    decompress file\n");
    fprintf(stderr,"\t-l    locate occurrences\n");
    fprintf(stderr,"\t-m    multiple search (kills option -osr)\n");
    fprintf(stderr,"\t-o    one line report (kills options -rls)\n");
    fprintf(stderr,"\t-p    *pattern* gives the pattern file name\n");    
    fprintf(stderr,"\t-r    report the position of the occurrences\n");
    fprintf(stderr,"\t-s    display len chars sourronding each occ ");
    fprintf(stderr,"(kills options -rl)\n");
    fprintf(stderr,"\t-v    produces a verbose output\n");
    fprintf(stderr,"\t-c    percentage of cached buckets (default 0)\n"); 
    fprintf(stderr,"\t-x    data access: 1=Ext, 2=Mmap, 3=In (default 1)\n\n");
    exit(1);
  }

  /* ----------------- read options --------------------- */
  Verbose=0;
  infile_name=NULL;
  num_opt = opterr=0;
  pattern_from_file = 0;
  multiple_search = 0;
  Type_mem_ops = EXT_MEM;
  Cache_percentage = 0;
  Report_occ = 0;
  Locate_occ = 0;
  Display_occ = 0;
  Oneline_report = 0;
  decompress = 0;

  while ((c=getopt(argc, argv, "dlomrpvs:c:x:")) != -1) {
    switch (c)
    {
      case 'd':
	decompress=1;
      case 'v':
        Verbose++; break;
      case 'l':
        Locate_occ=1; break;
      case 'r':
        Report_occ=1; break;
      case 'x':
        Type_mem_ops = atoi(optarg); break;
      case 's':
        Display_occ = atoi(optarg); break;
      case 'c':
        Cache_percentage = atof(optarg); break;
      case 'p':
        pattern_from_file = 1; break;
      case 'm':
        multiple_search = 1; break;
      case 'o':
        Oneline_report = 1; break;
      case '?':
        fprintf(stderr,"Unknown option: %c -main-\n", optopt);
        exit(1);
    }
    num_opt++;
  }

  if(decompress) {
    /* ------- read filename only ----- */
    infile_name = (char *) argv[optind];
    pattern = NULL;
  }
  else {
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
  if(Report_occ) 
    Locate_occ = 1;  // report implies locate

  if(Display_occ)
    Report_occ=Locate_occ = 0; //display text kills report and locate

  if(Oneline_report) {
    Display_occ=Report_occ=Locate_occ = 0;
  }
 
  if(Verbose>1) {
    fprintf(stderr,"\n*****************************************************");
    fprintf(stderr,"\n             bwsearch  Ver 1.0\n");
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

  if(decompress) 
    unbwi(infile_name);
  else if(multiple_search)
    multi_search(pattern);
  else
    single_search(pattern,pattern_from_file);

  /* ---------  Report Cache Usage Information ------------- */
  if (Verbose) 
    report_bwi_cache_usage();

  my_fclose(Infile);
  return 0;
}

