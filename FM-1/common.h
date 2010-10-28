/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   bwt-based compression and indexing

   common.h
   P. Ferragina & G. Manzini, 10 June 2000
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef DEBUG
#define DEBUG 1   /* set DEBUG to 0 to remove assertions and extra checks */
#endif
#if !DEBUG
#define NDEBUG 1  /* do not compile assertions */
#endif
#include <assert.h>


#define SIZE_BASIC_PROLOGUE  55   /* #bytes except char prefix_count */

#define MIN(a, b) ((a)<=(b) ? (a) : (b))


/* **************************************************************
   Type and structure definitions

   (some data structures have been overestimated to reduce dynamic
    memory allocation and since they are used at construction time.)
   ************************************************************** */
typedef unsigned int uint32;
typedef unsigned char uchar;

typedef struct {
  int *occ;               // occ chars of compact alph in prev. superbuc.
  int alpha_size;         // actual size of alphabet in this superbucket 
  uchar *bool_char_map;   // boolean map of chars occurring in this superbucket
} bucket_lev1;

typedef struct {
  uchar *text;          /* input text */
  int text_size;        /* size of text */
  int alpha_size;       /* actual size of (compact) alphabet in input text */
  int *sa;              /* suffix array */
  uchar *bwt;           /* BWT of input text */
  int bwt_eof_pos;      /* position of EOF within BWT */
  uchar bool_char_map[256];  /* ASCII[i] occurs iff bool_char_map[i] = 1 */
  uchar char_map[256];  /* ascii -> compact */
  int pfx_char_occ[256];/* entry i stores # of occ of chars 0..i-1 */
  bucket_lev1 *buclist_lev1;  /* array of num_bucs_lev1 superbuckets */
  int *start_lev2;      /* starting position of buckets */
  int *loc_occ;         /* locations of occurrences of the chosen_char */
  uchar chosen_char;    /* marked char */
  int skip;             /* one every "skip" occ of chosen_char is stored */
} bwi_input;

typedef struct {
  uchar *text;         /* input text */
  int text_size;       /* size of text */
  int alpha_size;      /* actual size of alphabet in input text */
  int *lf;             /* lf-mapping */
  uchar *bwt;          /* BWT of input text */
  int bwt_eof_pos;     /* position of EOF within BWT */
  bucket_lev1 *buclist_lev1;     /* array of num_bucs buckets*/
  uchar bool_char_map[256];
  int char_map[256];
  int inv_char_map[256];
  int bwt_occ[256];     /* entry i stores # of occ of chars 0..i-1 */
  int *start_lev2;      /* starting position of each buckets in compr file */
  int skip;             /* one occ of "chosen_char" every "skip"  */
  uchar chosen_char;    /* sampled char */
} bwi_out;


/* ================ "general" global variables ============ */
int Verbose;             /* produce verbose output */
int Type_compression;    /* indicates type of compression adopted */
int Is_dictionary;       /* a dictionary has to be processed */
int Is_URL;              /* a dictionary of URLS has to be processed */
int Is_huffword;         /* compression output of 7-bit huffword */
int Bucket_size_lev1;    /* size of superbucket (multiple of 1K) */
int Bucket_size_lev2;    /* size of bucket (multiple of 1K & divides _lev1) */
int Mtf_save;            /* # of copied MTF element */
int Num_bucs_lev1;       /* overall number of superbuckets */
int Num_bucs_lev2;       /* overall number of buckets */
int Start_prologue_info_sb; // starting byte of info superbuckets
int Start_prologue_info_b;  // starting byte of info buckets
int Start_prologue_occ;     // starting byte of char occurrences
// int Retrieved_occ;          // Number of retrieved explicit pos
double Marked_char_freq;    // maximum frequency of the marked char 


/* ======= global variables for I/O ====== */ 
FILE *Infile;            /* input file */
FILE *Outfile;           /* output file */
int Infile_size;         /* size of input file */
int  Outfile_size;       /* size of the output file */

FILE *Infile_head;            /* header file of a dictionary */
FILE *Infile_dict;            /* dictionary file */


/* ======= global variables for file management ====== */ 
uchar *File_start;  // byte where mapped file starts in memory
uchar *File_end;    // byte where mapped file starts in memory
uchar *File_pos;    // byte where read/write are currently positioned
int Type_mem_ops;   // It may assume one of the three values below

// ---- word search strategy
#define WSUBSTRING 0 // arbitrary substrings
#define WPREFIX    1 // word prefix
#define WSUFFIX    2 // word suffix
#define WFULL      3 // full word

// ---- startegy for accessing file in bwsearch
#define EXT_MEM  1  // We operate via fseek, fopen, getc, putc   
#define EXT_MMAP 2  // We map the file via mmap() function
#define IN_MEM   3  // We load the whole files in internal memory

// ---- algorithms for bucket compression ----
#define ARITH     1 // arithmetic coding 
#define HIER3     2 // hierarchical 3 level
#define UNARY     3 // usary coding
#define MULTIH    4 // huffman with multiple tables

// ----  constants used for searching in bwi files (see search_main()) 
#define WHAT_CHAR_IS 1    // used to retrieve the char in a given pos
#define COUNT_CHAR_OCC 2  // used to count char-occs before a given pos
#define NULL_CHAR 0       // used to skip the count of char-occs


// *****************************************************************
// this is a macro identical (hopefully) to the function my_getc()
// *****************************************************************
int my_getc_macro_tmp;
#define my_getc_macro(c,f)   \
if(1) {                \
  switch (Type_mem_ops)\
    {                  \
    case EXT_MEM:      \
      my_getc_macro_tmp = getc(f);  \
      if(my_getc_macro_tmp==EOF) {  \
         fprintf(stderr,"Unexpected end of file -bit_read-\n"); \
         exit(1);      \
      }                \
      c= my_getc_macro_tmp;	    \
      break;           \
    case EXT_MMAP:     \
      c = *File_pos++; \
      break;           \
    case IN_MEM:       \
      c = *File_pos++; \
      break;           \
    default:           \
      fprintf(stderr,"Error in choosing memory management -- my_getc() --");\
      exit(1);         \
      break;           \
    }                  \
} else c=0  /* this is never executed */          


/* *************************************************************************
   the following is a copy of the procedure bit_read() transformed to a macro
   to avoid the overhead of passing parameters etc
   ************************************************************************* */
uint32 t_macro;
#define bit_read_macro(dest,n)                                    \
{                                                                 \
  /* --- read groups of 8 bits until size>= n --- */              \
  while(Bit_buffer_size<n) {                                      \
    my_getc_macro(t_macro,Infile);                                \
    Bit_buffer |= (t_macro << (24-Bit_buffer_size));              \
    Bit_buffer_size += 8;                                         \
  }                                                               \
  /* ---- write n top bits in u ---- */                           \
  dest = Bit_buffer >> (32-n);                                    \
  /* ---- update buffer ---- */                                   \
  Bit_buffer <<= n;                                               \
  Bit_buffer_size -= n;                                           \
}
/* ************************************************************
   this macro reads a single bit from Bit_buffer. 
   ************************************************************ */
#define single_bit_read_macro(dest)         \
{                                           \
  if(Bit_buffer_size==0) {                  \
    my_getc_macro(t_macro,Infile);          \
    dest = t_macro >>7;                     \
    Bit_buffer = t_macro << 25;             \
    Bit_buffer_size=7;                      \
  }                                         \
  else {                                    \
    Bit_buffer_size--;                      \
    dest = Bit_buffer >> 31;                \
    Bit_buffer <<= 1;                       \
  }                                         \
}
    

// -------------------------------------------------------
// finally, some prototypes common to many functions
// -------------------------------------------------------
void fatal_error(char *s);
void out_of_mem(char *s);
double getTime ( void );  
int int_log2(int);
int int_pow2(int);
int uint_read(void);
int bit_read(int);
uchar my_getc(FILE *);
int my_fseek(FILE *, long, int);
void init_bit_buffer(void);

void read_basic_prologue(bwi_out *);
uchar *read_pattern_from_file(char *, int *);
int bwsearch(bwi_out *, uchar *, int, int *, int *);
int get_occ_pos(bwi_out *, int);
int check_bwi_suffix(char *s);
int occ(bwi_out *,int,uchar);
void get_info_sb(bwi_out *,int,int *);
uchar get_info_b(bwi_out *,uchar,int,int *,int);

void init_bwi_cache(void);
void report_bwi_cache_usage(void);
void disable_bwi_cache(void);

void my_open_file(char *);
void my_fclose(FILE *);





