/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   filter and unsort a list of words 
   
   scramble.c
   Ver 1.0 (21-jun-00)
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

// ---- global prototypes and vars
void fatal_error(char *);
int getline(char **, int *, FILE *);
void sprand(long seed);
int randomint(int n);
int Verbose;


/* ***************************************************************
        Main procedure
   *************************************************************** */
int main(int argc, char *argv[])
{
  int getopt(int argc, char * const *argv, const char *options);
  extern int optind, opterr, optopt;
  int c, k,i,j,plen,wlen,minlen,maxlen,seed,flen;
  char *infile_name, **words, *t;
  FILE *f;

  if(argc<2) {
    fprintf(stderr, "Usage:\n\t%s -v ",argv[0]);
    fprintf(stderr, "[-m minlen][-M maxlen][-s seed] input_file\n\n");
    fprintf(stderr,"\t-m minlen  minimum word-len accepted (def. 1)\n");    
    fprintf(stderr,"\t-M maxlen  maximum word-len accepted (def. 9999)\n"); 
    fprintf(stderr,"\t-s seed    seed for random nums (def. 0=no scramble)\n");
    fprintf(stderr,"\t-v         produces a verbose output\n");
    fprintf(stderr,"\n");
    exit(0);
  }

  /* ---------------- read options --------------------- */
  Verbose = 0;
  minlen=1;
  maxlen = 9999;
  seed=0;
  opterr = 0;
  while ((c=getopt(argc, argv, "vm:M:s:")) != -1) {
    switch (c)
    {
    case 'm':
      minlen = atoi(optarg); break;
    case 'M':
      maxlen = atoi(optarg); break;
    case 's':
      seed = atoi(optarg); break;
    case 'v':
      Verbose++; break;
    case '?':
      fprintf(stderr,"Unknown option: %c -main-\n", optopt);
      exit(1);
    }
  }
  infile_name=NULL;
  if(optind<argc)
     infile_name=argv[optind];
  if(infile_name==NULL)
    fatal_error("You must supply a file name (main)");

  //--------- open file
  if((f=fopen(infile_name,"r"))==NULL) {
    fprintf(stderr,"Error opening file %s ", infile_name);
    perror("(main)");
    exit(1);
  }
  fseek(f,0,SEEK_END);
  flen = ftell(f);
  if(flen<1) 
    fatal_error("Empty file (main)");
  if(Verbose)
    fprintf(stderr, "file length :%d\n",flen);
  // --- alloc word list ----
  words = (char **) malloc(flen/2 * sizeof(char *));
  if(words==NULL)
    fatal_error("Out of memory (main)");

  //---- read words ------
  rewind(f);
  for(i=j=0; ;i++) {
    words[j]=NULL;
    wlen = 0;
    plen = getline(&words[j],&wlen,f);
    if(plen== -1) break;  // end of file
    if(plen!=(int) strlen(words[j])) {
      fprintf(stderr, "Invalid pattern in row %d (multi_search)\n",i);
      exit(1);
    }
    if( (plen-1 >= minlen) && (plen-1<=maxlen) ) {
      if(Verbose>1)
	fprintf(stderr,"%s",words[j]);
      j++;
    }
  }
  if(Verbose) 
    fprintf(stderr,"%d lines read, %d selected \n",i,j);

  // ---- scramble lines ----
  if(seed) {
    sprand(seed);                // init random number generator
    for(i=j-1;i>0;i--) {
      k=randomint(i+1);    /*  0 <= k <= i e' la destinazione di p[i] */
      t=words[i]; words[i]=words[k]; words[k]=t;
    }
  }

  // ---- print result ------
  for(i=0;i<j;i++)
    printf("%s",words[i]);
  return 0;
}

// --- print error message and exit
void fatal_error(char *s)
{
  fprintf(stderr,"%s\n",s);
  exit(1);
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   portable random number generator
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

/* =============================================================
   global vars
   ============================================================= */
#define PRANDMAX 1000000000L   /* max random number+1 */
long arr[55];
int axx,bxx;
/* *********** return a random number between 0 and PRANDMAX-1 ********* */
long lprand(void)
{
  long t;

  if(axx-- == 0)
    axx=54;
  if(bxx-- == 0)
    bxx=54;
  t = arr[axx] - arr[bxx];
  if(t<0)
    t+= PRANDMAX;
  arr[axx] = t;
  return t;
}

/* ************ init the random number generator *********** */
void sprand(long seed)
{
  long lprand(void);
  int i, ii;
  long last, next;

  if(seed<0)
    seed = -seed;
  seed = seed % PRANDMAX;
  arr[0]=last=seed;
  next=1;
  for(i=1;i<55;i++) {
    ii = (21*i) % 55;
    arr[ii] = next;
    next = last - next;
    if(next<0)
      next += PRANDMAX;
    last = arr[ii];
  }
  axx = 0;
  bxx = 24;
  for(i=0;i<165;i++)
    last = lprand();
}

/* *******************************************************************
   generates a random integer between 0 and n-1
   ******************************************************************* */
int randomint(int n)
{
  long lprand(void);
  double f;
  int x;

  f = ((double) lprand())/PRANDMAX;
  x = (int) floor(f*n);
  if((x<0)||(x>=n))
    fatal_error("(randomint)");
  return x;
}

