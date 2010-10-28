/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   bwt-based compression and indexing
   
   P. Ferragina & G. Manzini, 10 June 2000
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "common.h"


// global vars used for suffix array construction
int Use_larsson_sada;   // =1 to use larsson sadakane
char *Safile_name;      // name of possible sa file



/* ***************************************************************
        Main bwi procedure
   *************************************************************** */
int main(int argc, char *argv[])
{
  int getopt(int argc, char * const *argv, const char *options);
  void open_files(char *, char *, int);
  double getTime ( void );  
  int check_bwi_suffix(char *s);
  void compress_file(void);
  void decompress_file(void);
  extern char *optarg;
  extern int optind, opterr, optopt;
  int c, decompress, num_opt;
  double start,end;
  char *infile_name, *outfile_name;


  /* ---------------- initialize with default values --------------------- */
  Bucket_size_lev1 = 16*1024;  // 16Kb each superbucket
  Bucket_size_lev2 = 1*1024;     // 1Kb each bucket
  Mtf_save = 20;
  decompress = 0;
  Type_compression = 4;
  Verbose=0;
  Marked_char_freq = 0.02;
  infile_name = outfile_name = NULL;
  Use_larsson_sada=0;
  Type_mem_ops = EXT_MEM;  //forces to operate via I/O-ops
  Is_dictionary = 0;       // no dictionary compression
  Is_huffword = 0;         // no huffword compression
  Is_URL = 0;              // no URL compression

  /* -------------- print usage message ----------------- */
  if(argc<2) {
    fprintf(stderr, "Usage:\n\t%s [-dlvDHU] [-B Bsize] [-b bsize] ",argv[0]);
    fprintf(stderr,"[-f freq] [-m mtf_top]\n\t    [-o outfile] ");
    fprintf(stderr,"[-t 2|4] [infile] \n\n");
    fprintf(stderr,"\t-d          decompress (kills options -BDHbflmt)\n");
    fprintf(stderr,"\t-v          produces a verbose output\n");
    fprintf(stderr,"\t-l          use Larsson-Sadakane suffix sorting\n");
    fprintf(stderr,"\t-D          dictionary compression\n");
    fprintf(stderr,"\t-H          huffword file compression\n");
    fprintf(stderr,"\t-U          url compression\n");
    fprintf(stderr,"\t-B bsize    size (in Kb) of level 1 buckets ");
    fprintf(stderr,              "(default %d)\n",Bucket_size_lev1/1024);
    fprintf(stderr,"\t-b bsize    size (in Kb) of level 2 buckets ");
    fprintf(stderr,              "(default %d)\n",Bucket_size_lev2/1024);
    fprintf(stderr,"\t-f freq     frequency of marked char ");
    fprintf(stderr,              "(default %.2f)\n",Marked_char_freq); 
    fprintf(stderr,"\t-m mtf      # size of mtf list ");
    fprintf(stderr,              "(default %d)\n",Mtf_save);
    fprintf(stderr,"\t-o outfile  output file\n");
    fprintf(stderr,"\t-t 2|4      compression: 2=Hier, ");
    fprintf(stderr,              "4=MHuf (default 4)\n\n");
    fprintf(stderr,"If no input file name is given default is stdin.\n");
    fprintf(stderr,"If no output file name is given default is ");
    fprintf(stderr,"{infile}.bwi for compression\n");
    fprintf(stderr,"and {infile}.y for decompression, or stdout if ");
    fprintf(stderr,"input is stdin.\n\n");
    exit(0);
  }

  /* ------------- read options from command line ----------- */
  num_opt = opterr = 0;
  while ((c=getopt(argc, argv, "B:b:f:t:m:o:dvlDHU")) != -1) {
    switch (c)
    {
      case 'o':
        outfile_name = optarg; break;
      case 'B':
        Bucket_size_lev1 = atoi(optarg)*1024; break;
      case 'b':
        Bucket_size_lev2 = atoi(optarg)*1024; break;
      case 'm':
        Mtf_save=atoi(optarg); break;
      case 't':
        Type_compression=atoi(optarg); break;
      case 'f':
        Marked_char_freq=atof(optarg); break;
      case 'd':
        decompress=1; break;
      case 'D':
        Is_dictionary=1; break;
      case 'H':
        Is_huffword=1; break;
      case 'U':
        Is_URL=1; break;
      case 'l':
        Use_larsson_sada=1; break;
      case 'v':
        Verbose++; break;
      case '?':
        fprintf(stderr,"Unknown option: %c -main-\n", optopt);
        exit(1);
    }
    num_opt++;
  }
  if(optind<argc)
     infile_name=argv[optind];

  if(Type_compression==MULTIH)
    Mtf_save = 256;        // multihuff always use full mtf list   

  /* ------------- Check arguments ----------- */  
  if ((decompress) && (infile_name!=NULL))
    if(check_bwi_suffix(infile_name)==0) {
      fprintf(stderr,"The file name must end with .bwi -main-\n");  
      exit(1);
    }
  if(Bucket_size_lev2<=0) {
     fprintf(stderr,"The size of a bucket must be greater than 0! -main-\n");  
     exit(1);
  }
  if((Bucket_size_lev1<=Bucket_size_lev2)||
     (Bucket_size_lev1%Bucket_size_lev2) ) {
    fprintf(stderr,"The size of a bucket ");
    fprintf(stderr,"must divide the size of a superbucket -main-\n");  
    exit(1);
  }
  if(Mtf_save<=0 || Mtf_save>256) {
     fprintf(stderr,"The size of MTF list must be in the range 1-256");
     fprintf(stderr," -main-\n");  
     exit(1);
  }
  if(Marked_char_freq<0 || Marked_char_freq >1) {
     fprintf(stderr,"The frequency of the marked char must be in [0,1]");
     fprintf(stderr," -main-\n");  
     exit(1);
  }

  if((!decompress) && ((Type_compression!=2) && (Type_compression!=4))) {
     fprintf(stderr,"Available compression algorithms:\n\t"); 
     fprintf(stderr,"2=3-level Hierarchical, 4=MultiTableHuffman -main-\n");  
     exit(1);
  }


  if(Verbose>1) {
    fprintf(stderr,"\n*****************************************************");
    fprintf(stderr,"\n             bwi  Ver 2.0\n");
    fprintf(stderr,"Created on %s at %s from %s\n",__DATE__,__TIME__,__FILE__);
    fprintf(stderr,"*****************************************************\n");
  }
  if(Verbose) {
    fprintf(stderr,"Command line: ");
    for(c=0;c<argc;c++)
      fprintf(stderr,"%s ",argv[c]);
    fprintf(stderr,"\n");
  }
  if(Verbose) {
    fprintf(stderr,"Superbucket size %dK, ",Bucket_size_lev1/1024);
    fprintf(stderr,"Bucket size %dK, ",Bucket_size_lev2/1024);
    fprintf(stderr,"Marked chars %.2f%%\n",Marked_char_freq*100);
  }
  if(Verbose) {
    fprintf(stderr,"Compression method: ");
    switch(Type_compression)
    {
    case HIER3: 
      fprintf(stderr ,"Hierarchical 3-level coding. "); break;
    case MULTIH: 
      fprintf(stderr,"Huffman with multiple tables. "); break;
    default: 
      fprintf(stderr,"Error: invalid compression method\n"); exit(1);
    }
    fprintf(stderr,"Size of mtf list %d.\n",Mtf_save);
  }

  start = getTime();
  open_files(infile_name,outfile_name,decompress);
  if(decompress)
    decompress_file();
  else
    compress_file();
  end = getTime();
  fprintf(stderr,"  Total elapsed time --> %.2f seconds.\n", end-start);

  if(!decompress) {
    fseek(Outfile,0,SEEK_END);
    fprintf(stderr,"  Input: %d bytes --> Output %d bytes.\n",
            Infile_size, (int) ftell(Outfile));
    fprintf(stderr,"  Overall compression --> %.2f%% (%.2f bits per char).\n",
     (100.0*ftell(Outfile))/Infile_size,(ftell(Outfile)*8.0)/Infile_size);
  }
  fprintf(stderr,"\n");
  return 0;  
}

/* ***************************************************************
   open input and output files, and the file for the suffix array
   *************************************************************** */
void open_files(char *infile_name, char *outfile_name, int decompress)
{
  FILE *my_fopen(const char *path, const char *mode);

  /* ------ open input and output files ------ */
  if(infile_name==NULL)
    Infile=stdin;
  else {                          
    Infile=my_fopen( infile_name, "rb"); // b is for binary: required by DOS
    if(Infile==NULL) {
      fprintf(stderr,"Unable to open file %s!\n",infile_name);
      exit(1);
    }
  }

  if ( outfile_name==NULL && infile_name==NULL) {
    Outfile=stdout;
    if(isatty(fileno (Outfile))) {
      fprintf(stderr,"I can't write my output to a terminal!\n");
      exit(1);
    }
  }
  else {
    if(outfile_name==NULL) {
      /* add ".bwi" for compression ".y" for decompression */
      outfile_name = (char *) malloc(strlen(infile_name)+5);
      outfile_name = strcpy(outfile_name, infile_name);
      if(decompress)
        outfile_name = strcat(outfile_name, ".y");
      else
        outfile_name = strcat(outfile_name, ".bwi");
    }
    Outfile = fopen( outfile_name, "wb"); // b is for binary: required by DOS
    if(Outfile==NULL) {
      fprintf(stderr,"Unable to open file %s!\n",outfile_name);
      exit(1);
    }
  }

  /* ---- build the file name for the suffix array (if any) ----- */
  if(!decompress) {
    if(infile_name!=NULL) {
      Safile_name = (char *) malloc(strlen(infile_name)+4);
      Safile_name = strcpy(Safile_name, infile_name);
      Safile_name = strcat(Safile_name, ".sa");
    }
    else
      Safile_name=NULL;
  }
}

























