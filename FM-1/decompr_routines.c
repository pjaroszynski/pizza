/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   bwt-based compression and indexing
   Routines for undoing run-length encoding and uncompress buckets
   decompr_routines.c   -----  P. ferragina & G. Manzini, 10 June 2000

       uncompress_bucket_unary();
       uncompress_bucket_hierarchy();
       uncompress_bucket_arith();
       unrle_only();
       decode_unary();
       decode_hierarchical();
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "common.h"


#if 0
/* ************************************************************
   uncompress the bucket which starts at the current position
   of infile. The bucket (without compression and rle) is "len" 
   bytes long and should be written in array dest[]. The
   bucket is read from the proper stream using bit_read().
   ************************************************************ */
void uncompress_bucket_unary(uchar *dest, int len, int alpha_size)
{
  int int_log2(int);
  void init_bit_buffer(void);
  int bit_read(int);
  static __inline__ int decode_unary(void);
  int k,j,i,rle,bits_x_char,rank,local_alpha_size,mtf_size;
  uchar c,mtf[256], inv_local_map[256]; 

  /* ---------- read local boolean map and compute inverse map ------ */
  local_alpha_size=0;
  for(k=0;k<alpha_size;k++)     
    if(bit_read(1)) 
      inv_local_map[local_alpha_size++] = k;

  /* --------- read initial status of the mtf list ---------- */
  mtf_size = MIN(Mtf_save,local_alpha_size);
  bits_x_char = int_log2(local_alpha_size);   // read initial mtf list
  for(k=0;k<mtf_size;k++) {
    mtf[k] = bit_read(bits_x_char);
  }
  
  /* ------ uncompress -------------------- */
  for(j=0;j<len; ) {
    rank = decode_unary();
    assert(rank<=mtf_size); 
    if(rank==0) {
      dest[j++]=mtf[0];       // rank zero: top of mtf list
      rle = bit_read(8);      // get # of following zeroes
      for(i=0;i<rle;i++)
	dest[j++]=mtf[0];
    }
    else if(rank<mtf_size) {
      dest[j]=mtf[rank];     // decode mtf rank 
      for(i=rank;i>0;i--)    // update mtf list
        mtf[i]=mtf[i-1];
      mtf[0]=dest[j++];      // update mtf[0] and j
    }
    else {                   // rank==mtf_size
      dest[j]=bit_read(bits_x_char);  // get char from file
      for(i=mtf_size-1;i>0;i--)       // update mtf
	mtf[i]=mtf[i-1];
      mtf[0]=dest[j++];     // update mtf[0] and j 
    }
  }

  /* ----- remap the bucket according to the inverse local map --- */
  for(j=0;j<len;j++) {             
    c=dest[j];                    
    assert(c<local_alpha_size);
    dest[j]=inv_local_map[c];   // expressed in the alphabet of superbucket
    assert(dest[j]<alpha_size);     
  }
}
#endif

/* ************************************************************
   uncompress the bucket which starts at the current position
   of infile. The bucket (without compression and rle) is "len" bytes 
   long and should be written in array dest[]. The bucket is
   read from the proper stream using bit_read().
   ************************************************************ */
int Mtf10log;
void uncompress_bucket_hierarchy(uchar *dest, int len, int alpha_size)
{
  int int_log2(int);
  void init_bit_buffer(void);
  int bit_read(int);
  int decode_hierarchical(void);
  int k,j,i,rle,bits_x_char,rank,local_alpha_size,mtf_size;
  uchar c,mtf[256], inv_local_map[256]; 

  /* ---------- read local boolean map and compute inverse map ------ */
  local_alpha_size=0;
  for(k=0;k<alpha_size;k++)     
    if(bit_read(1)) 
      inv_local_map[local_alpha_size++] = k;

  /* --------- read initial status of the mtf list ---------- */
  mtf_size = MIN(Mtf_save,local_alpha_size);
  bits_x_char = int_log2(local_alpha_size);   // read initial mtf list
  for(k=0;k<mtf_size;k++) {
    mtf[k] = bit_read(bits_x_char);
  }
  if(mtf_size>10)
    Mtf10log = int_log2(mtf_size-10); // constant used by decode_hierarchical
  
  /* ------ uncompress -------------------- */
  for(j=0;j<len; ) {
    rank = decode_hierarchical();
    assert(rank<=mtf_size); 

    if(rank==0) {
      dest[j++]=mtf[0];       // rank zero: top of mtf list
      rle = bit_read(8);      // get # of following zeroes
      for(i=0;i<rle;i++)
	dest[j++]=mtf[0];
    }
    else if(rank<mtf_size) { 
      dest[j]=mtf[rank];     // decode mtf rank 
      for(i=rank;i>0;i--)    // update mtf list
        mtf[i]=mtf[i-1];
      mtf[0]=dest[j++];      // update mtf[0] and j
    }
    else {                   // rank==mtf_size
      dest[j]=bit_read(bits_x_char);  // get char from file
      for(i=mtf_size-1;i>0;i--)       // update mtf
	mtf[i]=mtf[i-1];
      mtf[0]=dest[j++];     // update mtf[0] and j 
    }
  }

  /* ----- remap the bucket according to the inverse local map --- */
  for(j=0;j<len;j++) {             
    c=dest[j];                    
    assert(c<local_alpha_size);
    dest[j]=inv_local_map[c];   // expressed in the alphabet of superbucket
    assert(dest[j]<alpha_size);     
  }
}

#if 0
/* ************************************************************
   Uncompression procedure for a bucket compressed
   via an adpative 0-order arithmetic coder. "len" is the
   bucket length without compression and rle-encoding. The bucket
   is read from the proper stream.
   ************************************************************ */
void uncompress_bucket_arith(uchar *dest, int len, int alpha_size)
{
  void init_bit_buffer(void);
  void out_of_mem(char *);
  int int_log2(int);
  int unrle_only(uchar *, int, uchar *);
  int adapt_arith_decoder(uchar *,int);
  int bit_read(int);
  int h,k,j,i,bits_x_char,rank,local_alpha_size,mtf_size;
  int mtf_seq_len, rle_len;
  uchar c,mtf[256], inv_local_map[256], *rle, *mtf_seq; 

  /* ---------- read local boolean map and compute inverse map ------ */
  local_alpha_size=0;
  for(k=0;k<alpha_size;k++)     
    if(bit_read(1)) 
      inv_local_map[local_alpha_size++] = k;

  /* --------- read initial status of the mtf list ---------- */
  mtf_size = MIN(Mtf_save,local_alpha_size);
  bits_x_char = int_log2(local_alpha_size);   // read initial mtf list
  for(k=0;k<mtf_size;k++) {
    mtf[k] = bit_read(bits_x_char);
  }

  // Here we do not need to explicitely align to the byte, as in 
  // the compression step since the dummy-bits present in the 
  // buffer are automatically loosen by resetting the buffer.
  init_bit_buffer();

  /* ------ uncompress arithmetic code -------------------- */
  rle = (uchar *) malloc(2 * len * sizeof(uchar));
  if(rle == NULL) out_of_mem("uncompress_data");

  // Reads from file. Here +2 is a safe term for 0-bucket
  rle_len = adapt_arith_decoder(rle,local_alpha_size+2);

  /* ------ unrle -------------------- */

  mtf_seq = (uchar *) malloc(2 * len * sizeof(uchar));
  if(mtf_seq == NULL) out_of_mem("uncompress_data");

  mtf_seq_len = unrle_only(rle, rle_len, mtf_seq);

  /* ------ unmtf -------------------- */
  for(j=0,h=0;j<mtf_seq_len; ) {
    rank = mtf_seq[j++];

    if(rank<Mtf_save) {
      dest[h]=mtf[rank];     // decode mtf rank 
      for(i=rank;i>0;i--)    // update mtf list
        mtf[i]=mtf[i-1];
      mtf[0]=dest[h++];      // update mtf[0] and h
    }
    else {                   // rank==mtf_size

      dest[h]=mtf_seq[j++]; // get char stored explicitely
      for(i=mtf_size-1;i>0;i--)       // update mtf
	mtf[i]=mtf[i-1];
      mtf[0]=dest[h++];     // update mtf[0] and h 
    }
  }

  /* ----- remap the bucket according to the inverse local map --- */
  for(j=0;j<len;j++) {             
    c=dest[j];                    
    assert(c<local_alpha_size);
    dest[j]=inv_local_map[c];   // expressed in the alphabet of superbucket
    assert(dest[j]<alpha_size);     
  }

  free(mtf_seq);
  free(rle);
}


/* ********************************************************************
      Expand Run-length Encoding: 0^m --> 0 foll'd by Gamma-encoding of m 
   ******************************************************************* */
int unrle_only(uchar *in, int len, uchar *out)
{
  int int_pow2(int);
  int i,j,z,len_repr,len_out;

  len_out = 0;
  i = 0;
  while (i<len) 
    {
      if(in[i]==0){
	i++;
	for(len_repr=1; in[i] == 0;i++)  // it's represented len
	  len_repr++;

	assert(in[i]==1);
	z = int_pow2(len_repr-1);   // value of 0-run
	i++;     // skip the 1_escape bit

	for(j=0; j < len_repr - 1; j++) // decode the binary repr.
	  z += in[i++]<<j;

	for(j=0; j<z;j++)  // 0-run
	  out[len_out++] = (uchar) 0;
      } else {
	out[len_out++] = (uchar) in[i++];
      }
    }

  return(len_out);
}
#endif

// ----------------------------------------------------------
// these external variables are required for the use of 
// the macro bit_read_macro() and single_bit_read()
// -----------------------------------------------------------
extern uint32 Bit_buffer;       
extern int  Bit_buffer_size;  

#if 0
/* *********************************************************
   decode a unary code read from Infile. the output is the
   # of zeroes we see before we see a 1 
   1 --> 0
   01 --> 1
   001 --> 2  
     etc.
   ********************************************************* */
static __inline__ int decode_unary(void)
{
  int t,i=0;

  do {
    single_bit_read_macro(t);
    if(t!=0) break;
    i++;
  } while(1);
  return i;
}
#endif

/* *********************************************************
   Decode the 3-level Fenwick's code
   ********************************************************* */
int decode_hierarchical(void)
{
  int rank,tmp;
  
  bit_read_macro(rank,2);  // rank in [0..2]

  if (rank == 3) {   // escape code for first level
    bit_read_macro(tmp,3);
    rank +=  tmp;   // rank in [3..9]
  }

  if (rank == 10) { // escape code for second level
    bit_read_macro(tmp,Mtf10log);
    rank += tmp; // full code
  }

  return rank;
}




