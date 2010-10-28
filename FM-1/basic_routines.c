/* **************************************************************
   Basic_routines.c
   P. Ferragina & G. Manzini, 10 June 2000
   ************************************************************** */
#include "common.h"


int check_bwi_suffix(char *s)
{
  int len;

  if (s == NULL)
    return 0;
  len = strlen(s);
  if (len <= 4)
    return 0;
  return(!strcmp(s + len - 4, ".bwi"));

}

void out_of_mem(char *s)
{
  fprintf(stderr,"Ouf of memory in function %s!\n",s);
  exit(1);
}

void fatal_error(char *s)
{
  fprintf(stderr,"%s",s);
  exit(1);
}

double getTime ( void )
{
   double usertime,systime;
   struct rusage usage;

   getrusage ( RUSAGE_SELF, &usage );

   usertime = (double)usage.ru_utime.tv_sec +
     (double)usage.ru_utime.tv_usec / 1000000.0;

   systime = (double)usage.ru_stime.tv_sec +
     (double)usage.ru_stime.tv_usec / 1000000.0;

   return(usertime+systime);
}


int int_log2(int u)    // compute # bits to represent u
{
  int i = 1;
  int r = 1;
  
  while((i<=32) && (r<u)){
    r=2*r+1;
    i = i+1;
  }
    
  assert(i<=32);
  return i;
}

int int_pow2(int u)    // compute 2^u
{
  int i,val;

  for(i=0,val=1; i<u; i++)
    val *= 2;

  assert(i==u);
  return val;

}
       

/* *****************************************************
   Basic procedures to read and write bits using a Bit_buffer. 
   Unread/unwritten bits of Bit_buffer are the most significant ones.
   ***************************************************** */

uint32 Bit_buffer;       
int  Bit_buffer_size;  /* number of unread/unwritten bits in Bit_buffer */

// ***** Initializes the Bit_buffer to zero 
void init_bit_buffer(void)
{
  Bit_buffer= (uint32) 0;
  Bit_buffer_size=0;
}



// ***** Write n (<= 24) bits taken from vv. The content of 
// ***** Bit_buffer is flushed out until it contains <8 bits 
  
void bit_write24(int n, int vv)  
{                              // v contains bits to read starting from
                               // the least significant bits
  uint32 v = (uint32) vv;

  assert(Bit_buffer_size<8);
  assert(n>0 && n<=24);
  assert( v < 1u << n );
   
  /* ------- add n bits to Bit_buffer -------- */
  Bit_buffer_size += n;       // add first, to compute the correct shift
  Bit_buffer |= (v << (32 - Bit_buffer_size));  // compact to end of the buffer

  /* ------- flush Bit_buffer as much as possible ----- */
  while (Bit_buffer_size>=8) {            
    if( putc((Bit_buffer>>24),Outfile) == EOF) {
      fprintf(stderr,"Error writing to output file -bit_write-\n");
      exit(1);
    }
    Outfile_size++;
    Bit_buffer <<= 8;                       
    Bit_buffer_size -= 8;                 
  }                                           
} 

// ****** Write in Bit_buffer n bits taken from vv (possibly n > 24) 
void bit_write(int n, int vv)
{  
  void bit_write24(int n, int vv);  
  uint32 v = (uint32) vv;

  assert(n <= 32);
  if (n > 24){
    bit_write24(n-24, (v>>24) & 0xffL);
    bit_write24(24, v & 0xffffffL);
  } else {
    bit_write24(n,v);
  }

}

// ****** Write in Bit_buffer four bytes 
void uint_write(int uu)
{  
  void bit_write(int n, int vv);
  uint32 u = (uint32) uu;

  bit_write(8, (u>>24) & 0xffL);
  bit_write(8, (u>>16) & 0xffL);
  bit_write(8, (u>> 8) & 0xffL);
  bit_write(8,  u      & 0xffL);

}


// ***** write n using the 7x8 scheme. The most
// ***** significant bit indicates if the representation extends
// ***** to the successive bytes.

void write7x8(int nn)
{
  void bit_write(int n, int vv);
  uint32 t;
  uint32 n = (uint32) nn;

  do {
    t = n & 0x7f;       // takes last seven bits
    n = n >> 7;
    if(n>0) t |= 0x80;  // set 8th bit
    bit_write(8,t);
  } while(n>0);
}


// *****  Complete with zeroes the first byte of Bit_buffer 
// *****  This way, the content of Bit_buffer is entirely flushed out

void bit_flush( void )
{
  if(Bit_buffer_size!=0)
    bit_write(8 - (Bit_buffer_size%8) , 0);  // pad with zero !
}





// ********************************************************************
// ********* Procedures to access the memory in various ways  *********
// ********************************************************************

// ***** Return a character taken from the proper stream
// ***** This stream can be: I/O, MMap, Internal Memory 
// ***** and it is indicated in the global variable Type_mem_ops

static int my_getc_tmp;
__inline__ uchar my_getc(FILE *f)
{

  switch (Type_mem_ops)
    {
    case EXT_MEM: 
      if ((my_getc_tmp=getc(f)) == EOF){  
	fprintf(stderr,"Unexpected end of file -my_getc-\n");
	exit(1);
      }
      return my_getc_tmp;
    case EXT_MMAP:
      return *File_pos++;
    case IN_MEM:
      return *File_pos++;
    default:
      fprintf(stderr,"Error in choosing memory management! -my_getc()-\n");
      exit(1);
    }
  fprintf(stderr,"Error in the code! -my_getc()-\n");
  exit(1);
  return 0;  
}


// ***** Moves over the proper stream.
// ***** This stream can be: I/O, MMap, Internal Memory 
// ***** and it is indicated in the global variable Type_mem_ops

int my_fseek(FILE *f, long offset, int whence)
{
  int res = 0;    // code for no error

  switch (Type_mem_ops)
    {
    case EXT_MEM: 
      res = fseek(f,offset,whence);
      break;

    case EXT_MMAP:
      if ((whence == SEEK_SET) && (offset < File_end - File_start + 1))
	{ File_pos = File_start + offset;}
      else if ((whence == SEEK_END) &&  (offset < File_end - File_start + 1))
	{ File_pos = File_end - offset; }
      else res = 1;  // error
      break;

    case IN_MEM:
      if ((whence == SEEK_SET) && (offset < File_end - File_start + 1))
	{ File_pos = File_start + offset; }
      else if ((whence == SEEK_END) &&  (offset < File_end - File_start + 1))
	{ File_pos = File_end - offset; }
      else res = 1;  // error
      break;

    default:
      fprintf(stderr,"Error in choosing memory management -my_fseek-");
      exit(1);
    }
  return res;
}

// ***** fopen, fread, or mmap a file.
// ***** The effect depends on the chosen memory management
// ***** according to what is stored in Type_mem_ops
 
FILE *my_fopen(const char *path, const char *mode)
{
  FILE *res;
  int len;

  if ((res = fopen(path,mode)) == NULL){
    fprintf(stderr,"Error opening file %s ", path);
    perror("(my_fopen)");
    exit(1);
  }
  fseek(res,0L,SEEK_END);  // compute file length
  len = ftell(res);

  fseek(res,0L,SEEK_SET);  //rewind

  switch (Type_mem_ops)
    {
    case EXT_MEM:
      break;
      
    case EXT_MMAP:
      // fd = open(path,O_RDONLY);
      File_start = (uchar *) mmap(0,len,PROT_READ,MAP_SHARED,fileno(res),0);
      File_pos = File_start;
      File_end = File_start + len - 1;
      break;

    case IN_MEM:
      if ((File_start = (uchar *) malloc(len)) == NULL){
	fprintf(stderr,"Error in allocating memory -my_fopen-");
	exit(1);
      }
      if (len != (int) fread(File_start,sizeof(uchar),len,res)) {
	// (int) read(fileno(res),File_start,len)) {
	fprintf(stderr,"Error in reading form file -my_fopen-");
	exit(1);
      }
      File_pos = File_start;
      File_end = File_start + len - 1;
      break;

    default:
      fprintf(stderr,"Error in choosing memory management -my_fopen-");
      break;
    }

  return res;

}

// ***** fclose or munmap a file.
// ***** The effect depends on the chosen memory management
// ***** according to what is stored in Type_mem_ops
 
void my_fclose(FILE *f)
{

  switch (Type_mem_ops)
    {
    case EXT_MEM:
      break;
      
    case EXT_MMAP:
      munmap(File_start, (int)(File_end - File_start + 1));
      break;

    case IN_MEM:
      free(File_start);
      break;

    default:
      fprintf(stderr,"Error in choosing memory management -- my_fclose() --");
      break;
    }

  if (fclose(f) == EOF){
    fprintf(stderr,"Error in closing the file -- my_fclose --");
    exit(1);
  }
}




// ***** Load Bit_buffer 8 bits at time, until it contains at 
// ***** least n (<=24) bits. Then return n bits from Bit_buffer.
// ***** This procedure exploites the function my_getc which
// ***** accesses the proper stream to which the file is mapped.

__inline__ int bit_read24(int n)
{
  uchar my_getc(FILE *f);
  uint32 t,u;

  assert(Bit_buffer_size<8);
  assert(n>0 && n<=24);

  /* --- read groups of 8 bits until size>= n --- */
  while(Bit_buffer_size<n) {
    t = (uint32) my_getc(Infile);
    Bit_buffer |= (t << (24-Bit_buffer_size));
    Bit_buffer_size += 8;
  }
  /* ---- write n top bits in u ---- */
  u = Bit_buffer >> (32-n);
  /* ---- update buffer ---- */
  Bit_buffer <<= n;
  Bit_buffer_size -= n;
  return((int)u);
}


// ***** Return 32 bits taken from Bit_buffer
int uint_read(void)
{
  int bit_read(int n);  
  uint32 u;

  u =  bit_read(8)<<24;
  u |= bit_read(8)<<16;
  u |= bit_read(8)<<8;
  u |= bit_read(8);
  return((int)u);
}

// ****** Read n bits from Bit_buffer 
__inline__ int bit_read(int n)
{  
  int bit_read24(int n);  
  uint32 u = 0;

  assert(n <= 32);
  if (n > 24){
    u =  bit_read24(n-24)<<24;
    u |= bit_read24(24);
    return((int)u);
  } else {
    return(bit_read24(n));
  }

}

// ***** Return an integer coded according to a var-length
// ***** representation which is byte-aligned. The most
// ***** significant bit indicates if the representation extends
// ***** to the successive bytes.

int read7x8(void)
{
  int bit_read(int);
  int i;
  uint32 ris,flag,n;


  ris=0;
  flag=1;
  for(i=0;flag!=0;i+=7) {
    n= (uint32) bit_read(8);   // read eigth bit
    flag= n & 0x80;  // continuation bit on !?  
    n &= 0x7f;      // set the 8th bit to off  
    assert(i<=28);
    ris |= n << i;     // add 7 bits (from most to less signif.)
  }
  return ris;
}


      
/* ************************************************************
   read a pettren from file "pattern_file_name". 
   return a pointer to the pattern and store the lenght in *pat_len
   *********************************************************** */
uchar *read_pattern_from_file(char *pattern_file_name, int *pat_len)
{
  FILE *pattern_file;
  uchar *pattern;
  int i;

  // -------- open pattern file
  assert(pattern_file_name!=NULL);
  pattern_file=fopen(pattern_file_name, "rb");
  if(pattern_file == NULL) {
    fprintf(stderr,"Unable to open file %s!\n",pattern_file_name);
    perror("(read_pattern_from_file)");
    exit(1);
  }
  //------- read pattern
  fseek(pattern_file,0L,SEEK_END);
  *pat_len = ftell(pattern_file);
  assert(*pat_len > 0);
  pattern = (uchar *) malloc(*pat_len);
  if(pattern==NULL) 
    out_of_mem("read_pattern_from_file");
  fseek(pattern_file,0L,SEEK_SET);
  for(i=0; i< (int)(*pat_len); i++)
    pattern[i] = getc(pattern_file);
  return pattern;
}
      

// *****  Complete with zeroes the first byte of Bit_buffer 
//        This way, the content of Bit_buffer is entirely flushed out *****

void fbit_flush(FILE *f)
{
  void fbit_write24(FILE *f, int n, int vv);
  
  if(Bit_buffer_size!=0)
    fbit_write24(f, 8 - (Bit_buffer_size%8) , 0);  // pad with zero !
}


// ********* Write n (<= 24) bits taken from v. The content of 
//           Bit_buffer is flushed out until it contains <8 bits ********* 
  
void fbit_write24(FILE *f, int n, int vv)  
{                              // v contains bits to read starting from
                               // the least significant bits
  uint32 v = (uint32) vv;

  assert(Bit_buffer_size<8);
  assert(n>0 && n<=24);
  assert( v < 1u <<(n+1) );
   
  /* ------- add n bits to Bit_buffer -------- */
  Bit_buffer_size += n;       // add first, to compute the correct shift
  Bit_buffer |= (v << (32 - Bit_buffer_size));  // compact to end of the buffer

  /* ------- flush Bit_buffer as much as possible ----- */
  while (Bit_buffer_size>=8) {            
    if( putc((Bit_buffer>>24),f) == EOF) {
      fprintf(stderr,"Error writing to output file -fbit_write-\n");
      exit(1);
    }
    Bit_buffer <<= 8;                       
    Bit_buffer_size -= 8;                 
  }                                           
} 



// *** Load Bit_buffer 8 bits at time, until it contains at 
//     least n (<=24) bits. Then return n bits from Bit_buffer. ***

int fbit_read24(FILE *f, int n)
{
  int t0;
  uint32 t,u;

  assert(Bit_buffer_size<8);
  assert(n>0 && n<=24);

  /* --- read groups of 8 bits until size>= n --- */
  while(Bit_buffer_size<n) {
    if ((t0=getc(f)) == EOF){  
      fprintf(stderr,"Unexpected end of file -bit_read-\n");
      exit(1);
    }
    t = (uint32) t0;
    Bit_buffer |= (t << (24-Bit_buffer_size));
    Bit_buffer_size += 8;
  }
  /* ---- write n top bits in u ---- */
  u = Bit_buffer >> (32-n);
  /* ---- update buffer ---- */
  Bit_buffer <<= n;
  Bit_buffer_size -= n;
  return((int)u);
}


// ****** Write in file f the n bits taken from vv (possibly n > 24) 
void fbit_write(FILE *f, int n, int vv)
{  
  void fbit_write24(FILE *f,int n, int vv);  
  uint32 v = (uint32) vv;

  assert(n <= 32);
  if (n > 24){
    fbit_write24(f,n-24, (v>>24) & 0xffL);
    fbit_write24(f,24, v & 0xffffffL);
  } else {
    fbit_write24(f,n,v);
  }
}

// ****** Read n bits from file f 
int fbit_read(FILE *f,int n)
{  
  int fbit_read24(FILE *f,int n);  
  uint32 u = 0;

  assert(n <= 32);
  if (n > 24){
    u =  fbit_read24(f,n-24)<<24;
    u |= fbit_read24(f,24);
    return((int)u);
  } else {
    return(fbit_read24(f,n));
  }
}


/* ******************************************************************
   this function is a simple replacement of the getline() GNU procedure
   which is not available under windows.
   The input is a buffer "line" of size "max_len" and a file *f.
   The routins read from *f until an EOL is found and returns the
   chars read in line[] including the EOL and a null char. 
   max_len contains the number of chars read (excluding the terminating null)
   if EOF or an error occurs the function returns -1.
   ***************************************************************** */
int my_getline(char *line, int max_len, FILE *f)
{
  int c,i;

  for(i=0;i<max_len-1;i++) {
    c  = getc(f);
    if(c==EOF) return -1;     // error or end of file
    line[i]=(char) c;         // store char
    if(c=='\n') {             // end of line?
      line[++i] = '\0';
      return i;               // number of chars written excluding \0
    } 
  }  
  return max_len;             // line too long
}

















