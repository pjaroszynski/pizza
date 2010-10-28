/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   pattern search with plain suffix array
   
   sa.c
   Ver 1.0 (2-giu-00)
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "common.h"
#include "sacopy.h"
#include <limits.h>

#define use_my_getline 1     // set to 0 to use standard getline()

FILE *Textfile, *Safile;
int Pointer_size;
int Textfile_size, Verbose, Report_occ, Locate_occ;
int Use_larsson_sada;

static int fbit_read(FILE *,int);
static void fbit_flush( FILE * );


/* ***************************************************************
        Main procedure
   *************************************************************** */
int main(int argc, char *argv[])
{
  int getopt(int argc, char * const *argv, const char *options);
  void sa_multi_search(char *pfile_name);
  void sa_single_search(char *pattern, int pattern_from_file);
  void open_files(char *);
  extern int optind, opterr, optopt;
  char *textfile_name;
  int c, num_opt, pattern_from_file, multiple_search;
  uchar *pattern;

  if(argc<3) {
    fprintf(stderr, "Usage:\n\t%s [-clmprv] pattern input_file\n\n",argv[0]);
    fprintf(stderr,"\t-c      count the occurrences only (no locate)\n");    
    fprintf(stderr,"\t-l      use Larsson-Sadakane suffix sorting\n");    
    fprintf(stderr,"\t-m      multiple pattern search\n");    
    fprintf(stderr,"\t-r      report the position of the occurrences\n");
    fprintf(stderr,"\t-p      pattern is a filename containig the pattern\n");
    fprintf(stderr,"\t-v      produces a verbose output\n\n");
    fprintf(stderr,"\tSearch {pattern} in {input_file}. ");
    fprintf(stderr,"If a suffix array for {input_file}\n\tdoes not ");
    fprintf(stderr,"exist creates one in {input_file}.sa\n\n");
    exit(0);
  }

  /* ---------------- read options --------------------- */
  Verbose=0;
  Report_occ = 0;
  Locate_occ=1;
  Use_larsson_sada=0;
  textfile_name=NULL;
  pattern_from_file = 0;
  multiple_search = 0;
  num_opt = opterr = 0;

  while ((c=getopt(argc, argv, "cmrvpl")) != -1) {
    switch (c)
    {
      case 'c':
        Locate_occ = 0; break;
      case 'r':
        Report_occ = 1; break;
      case 'l':
        Use_larsson_sada = 1; break;
      case 'm':
        multiple_search = 1; break;
      case 'v':
        Verbose++; break;
      case 'p':
        pattern_from_file = 1; break;
      case '?':
        fprintf(stderr,"Unknown option: %c -main-\n", optopt);
        exit(1);
    }
    num_opt++;
  }
  /* ------- read pattern and filename --------- */
  if(optind!=argc-2) {
    fprintf(stderr,"You must supply a patter and a filename! -main-\n");
    exit(1);
  }
  else {
    pattern = (uchar *) argv[optind++];
    textfile_name = (uchar *) argv[optind];
  }

  // open Textfile and Safile
  open_files(textfile_name);
  // do the search
  if(multiple_search)
    sa_multi_search(pattern);
  else
    sa_single_search(pattern,pattern_from_file);

  return 0;
}


/* ****************************************************************
   this procedures read one line at a time from file pfile_name
   and search it in the input file
   **************************************************************** */
void sa_multi_search(char *pfile_name)
{
  void sa_bin_search(uchar *, int, int *, int *);
  int getline(char **, int *, FILE *);
  int my_getline(char *, int, FILE *);
  int get_suffix(int);  
  FILE *pfile;   // file containing the patterns
  char *pattern;
  int i, j, plen, max_aux, tot_located;
  int sp,ep,num_occ;
  double start,end;
  int max_plen=40;

  // open file containing patterns
  pfile = fopen(pfile_name,"r"); // this will be a text file
  if(pfile==NULL) {
    fprintf(stderr, "Error opening file %s ",pfile_name);
    perror("(sa_multi_search)");
    exit(1);
  }
  // alloc memory for pattern
  pattern = (char *) malloc(max_plen);
  if(pattern==NULL) out_of_mem("sa_multi_search");
  max_aux=max_plen;
  tot_located = 0;               // # of pattern located 

  // ---- loop until there are pattern in file
  start = getTime();
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
    sa_bin_search(pattern, plen-1, &sp, &ep);
    num_occ = ep-sp;
    if(Verbose>1) 
      fprintf(stderr,"Pattern \"%s\" occurs %d times\n",pattern,num_occ);
    // ----- locate occurrences
    if(num_occ > 0 && Locate_occ) {
      for(j=sp;j<ep;j++) get_suffix(j);
      tot_located += num_occ;
    }
  }
  end = getTime();
  if(i==0) {
    fprintf(stderr,"Nothing to search!");
    return;
  }
  fprintf(stderr,"%d patterns searched. %d occs located\n",i,tot_located);
  fprintf(stderr,"Tot. time (search+locate) %.2f secs.  ",end-start);
  fprintf(stderr,"Ave. time per search %.4f secs\n\n",(end-start)/i);
  return;
}


/* *****************************************************************
   this procedure searches all occurrences of pattern in the 
   input file. the exact position of each occurrence is always
   computed but they are actually reported only if Report_occ!=0.
   If pattern_from_file!=0 the pattern is interpreted as the name
   of the file containing the pattern
   ***************************************************************** */
void sa_single_search(char *pattern, int pattern_from_file)
{
  void sa_bin_search(uchar *, int, int *, int *);
  int get_suffix(int);  
  uchar *read_pattern_from_file(char *, int *);
  int pat_len,sp,ep,*list_occ,num_occ,i;
  double start,end;

  /* ----- pattern is read from a file? ----- */
  if(pattern_from_file)
    pattern = read_pattern_from_file((char *) pattern, &pat_len);
  else 
    pat_len = strlen(pattern);

  start = getTime();
  sa_bin_search(pattern, pat_len, &sp, &ep);
  end = getTime();
  num_occ = ep-sp;
  fprintf(stderr,"Found %d occ in %.2f seconds!\n",num_occ,end-start);
  
  /* ---------  list pattern occurrences ------------- */
  if (num_occ > 0 && Locate_occ) {
    list_occ = (int *) malloc(sizeof(int) * num_occ);
    if(list_occ == NULL)
      out_of_mem("sa_single_search (list_occ allocation failed!)");

    start = getTime();
    // Accumulate occurrences in an array
    for(i = sp; i < ep; i++)
      list_occ[i-sp] = get_suffix(i);
    end = getTime();

    if (Report_occ)
      for(i = 0; i < num_occ; i++)
	fprintf(stderr,"Occ. %4d is at position %8d\n",i+1,list_occ[i]);

    fprintf(stderr,"Retrieved %d occ's in %.2f seconds!\n\n",
	    num_occ,end-start);
  }
  else
    fprintf(stderr,"\n");
}





void sa_bin_search(uchar *s, int len, int *sp, int *ep)
{
  int get_suffix(int);
  int first,last,middle,pos,i,c;

  assert(len>0);
  last=Textfile_size;
  middle=0;
  first=0;
  // invariant: pattern can be only within [first,last)
  while(first<last) {
    middle = (first+last)/2;
    pos = get_suffix(middle);
    fseek(Textfile,pos,SEEK_SET);
    // compare s[] with text starting at pos
    for(i=0;i<len;i++) {
      c=getc(Textfile);
      if(c==EOF) c = -1;  // EOF is the smallest char
      if(c<s[i]) {
	first = middle+1; 
	break;
      }
      else if(c>s[i]) {
	last = middle;
	break;
      }
    }
    if(i==len) break;    // pattern found at middle. exit while 
  }
  // binary search concluded. find all occurrences of pattern
  if(first<last) {
    int ff,ll,mm;
    // --------- get first occ of pattern ----------
    ff=first; ll=middle;
    // invariance: the first occ is between ff and ll-1
    //             and there is an occ in ll 
    while(ff<ll) {
      mm = (ff+ll)/2;
      pos = get_suffix(mm);
      fseek(Textfile,pos,SEEK_SET);
      for(i=0;i<len;i++) {
	c=getc(Textfile);
	if(c==EOF) c = -1;  // EOF is the smallest char
	if(c<s[i]) break;
	else if(c>s[i]) {
	  fprintf(stderr,"suffix out of order -sa_bin_search-\n");
	  exit(1);
	}
      }
      if(i==len) ll=mm;
      else ff=mm+1;
    }
    first=ll;
    // --------- get last occ of pattern ----- 
    ff=middle+1; ll=last;
    // invariance: the last occ is between ff and ll-1
    //             and there is an occ in ff-1 
    while(ff<ll) {
      mm = (ff+ll)/2;
      pos = get_suffix(mm);
      fseek(Textfile,pos,SEEK_SET);
      for(i=0;i<len;i++) {
	c=getc(Textfile);
	if(c==EOF) c = -1;  // EOF is the smallest char
	if(c>s[i]) break;
	else if(c<s[i]) {
	  fprintf(stderr,"suffix out of order -sa_bin_search-\n");
	  exit(1);
	}
      }
      if(i==len) ff=mm+1;
      else ll=mm;
    }
    last=ff;
  }
  *sp=first;
  *ep=last;
  return;
}


#if 0
void print_occ(int sp, int ep)
{
  int i; 
  __inline__ int get_suffix(int);

  for(i=sp;i<ep;i++) {
    printf("Occ. %4d is at position %8d\n",1+i-sp,get_suffix(i));
  }
}
#endif

// ----- return  the position of the n-th suffix
__inline__ int get_suffix(int n)
{
  int pos,rem;

  assert(n>=0 && n<Textfile_size);
  pos = n*Pointer_size;
  rem = pos& 0x7;         // last 3 bits
  fseek(Safile, pos>>3, SEEK_SET);
  init_bit_buffer();
  if(rem)  fbit_read(Safile,rem);
  return fbit_read(Safile, Pointer_size);
}


/* ***************************************************************
   open input and output files
   *************************************************************** */
void open_files(char *textfile_name)
{
  int int_log2(int), n, sa_size;
  void create_sa(FILE *, FILE *, int);
  void out_of_mem(char *); 
  double getTime(void);
  char *safile_name;
  double start,end;
  long long sa_size_ll;

  /* --------------- open input file ----------------- */
  assert(textfile_name!=NULL);
  Textfile=fopen( textfile_name, "rb");
  if(Textfile==NULL) {
    fprintf(stderr,"Unable to open file %s!\n",textfile_name);
    exit(1);
  }
  if (fseek(Textfile, 0L, SEEK_END)) {
    perror(textfile_name);
    exit(1);
  }
  Textfile_size=ftell(Textfile);
  if (Textfile_size==0) {
    fprintf(stderr, "%s: file empty\n", textfile_name);
    exit(1);
  }
  Pointer_size = int_log2(Textfile_size);
  sa_size_ll = (((long long) Pointer_size)*Textfile_size+7)/8;
  sa_size = (int) sa_size_ll;

  /* ------------ open sa file ---------------- */
  safile_name = (char *) malloc(strlen(textfile_name)+4);
  if(safile_name==NULL) out_of_mem("open_files");
  safile_name = strcpy(safile_name, textfile_name);
  safile_name = strcat(safile_name, ".sa");
  Safile = fopen(safile_name, "a+b"); 
  n=ftell(Safile);

  if (n==0) {
    if(Verbose)
      fprintf(stderr, "Writing the suffix array to file %s\n",safile_name);
    start = getTime();
    create_sa(Textfile,Safile,Textfile_size);
    end = getTime();
    if(Verbose)
      fprintf(stderr, "Time for suffix array construction and writing: %f\n",
	      end-start);
    if (fseek(Safile, 0L, SEEK_END)) {
      perror(safile_name);
      exit(1);
    }
    n=ftell(Safile);
  }

  if(n!=sa_size) {
    fprintf(stderr, "Invalid suffix array file (%d vs %d) -open_files-\n",
	    n,sa_size);
    exit(1);
  }
  if(Verbose>1) 
    fprintf(stderr,"Compression ratio: %.2f%%\n",100+
	    (100.0*sa_size)/Textfile_size);
}


void create_sa(FILE *text, FILE *sa, int n)
{
  int i, psize, *p, *paux;
  int *create_sa_ls(FILE *, int);
  int *create_sa_5n(FILE *, int);
  int int_log2(int);

  if(Use_larsson_sada) {
    paux = create_sa_ls(text,n);
    p = paux+1;
  } else 
    p = paux = create_sa_5n(text,n);

  // write sa to file
  init_bit_buffer();
  psize = int_log2(n);
  for(i=0;i<n;i++)
    fbit_write(sa,psize,p[i]);
  fbit_flush(sa);

  free(paux);
  return;
}


int *create_sa_5n(FILE *text, int n)
{
  int scmp3(unsigned char *p, unsigned char *q, int maxl);
  void out_of_mem(char *);
  int t, *p, *suffixsort5n(char *,int);
  char *s;

   // ----- allocate
   s=malloc(n+N_OVERSHOOT);
   if (!s) out_of_mem("create_sa_5x");
   rewind(text);
   t=read(fileno(text),s,n);
   if(t!=n) {
     fprintf(stderr,"error reading input file -create_sa_5n-\n");
     exit(1);
   }
   p=suffixsort5n(s,n);
#if 0
   { int i;
   for (i=0; i<n-1; ++i)
     if (scmp3(s+p[i], s+p[i+1], 
                MIN(n-p[i], n-p[i+1]))>=0) {
       fprintf(stderr, "Suffix array check failed at position %d\n", i);
       exit(1);
     }
   }
#endif
   free(s);
   return p; 
}

#define COMPACT 2

int *create_sa_ls(FILE *text, int n)
{
   int i, *x, *p, *pi;
   int k, l;
   void out_of_mem(char *);
   void suffixsort(int *x, int *p, int n, int k, int l);
#if COMPACT==2
   unsigned char q[UCHAR_MAX+1];
#endif

   // ----- allocate
   p=malloc((n+1)*sizeof *p);
   x=malloc((n+1)*sizeof *x);
   if (! p || ! x) out_of_mem("create_sa_ls");

   // remap alphabet
#if COMPACT==1
   l=UCHAR_MAX;
   k=1;
   for (rewind(text), pi=x; pi<x+n; ++pi) {
      *pi=c=getc(text);
      if (c<l)
         l=c;
      if (c>=k)
         k=c+1;
   }
#else
   for (rewind(text), pi=x; pi<x+n; ++pi)
      *pi=getc(text);
#if COMPACT==0
   l=0;
   k=UCHAR_MAX+1;
#elif COMPACT==2
   for (i=0; i<=UCHAR_MAX; ++i)
      q[i]=0;
   for (pi=x; pi<x+n; ++pi)
      q[*pi]=1;
   for (i=k=0; i<=UCHAR_MAX; ++i)
      if (q[i])
         q[i]=k++;
   for (pi=x; pi<x+n; ++pi)
      *pi=q[*pi]+1;
   l=1;
   ++k;
#endif
#endif
   // compute sa
   suffixsort(x, p, n, k, l);
   // deallocate and return 
   free(x);
   return p;
}

#if DEBUG
/* ******* lexicographic string comparison ******* */ 
int scmp3(unsigned char *p, unsigned char *q, int maxl)
{
   int i;
   i = 0;
   while (maxl>0 && *p==*q) {
      p++; q++; i++;
      maxl--;
   }
   if (maxl>0) return *p-*q;
   return q-p;
}
#endif












