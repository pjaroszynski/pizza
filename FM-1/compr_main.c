/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   bwt-based compression and indexing
   Main compression procedures
   compr_main.c   -----  P. Ferragina & G. Manzini 2000
   This files contains most of the routines used during compression,
   the only excpetions are the routines for multiple huffman coding
   which are in multihuf.c 

       compress_file();
       compute_info_buckets();
       compress_superbucket();
       compress_bucket();
       compute_locations();
       compute_locations_dict();
       compute_ranks_dict();
       compute_locations_huffword();
       mtf_string();
       write_prologue();
       write_susp_infos();
       read_text();
       remap_alphabet();
       build_sa();
       compute_bwt();
       rle_hierarchical();
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "common.h"
#include "sacopy.h"

#define EOF_shift(n) (n < s->bwt_eof_pos) ? n+1 :  n

/* =============================================================
   strucuture used by compute_locations_dict() and 
   compute_locations_huffword() to store the text position and
   bwt position of a char
   ============================================================= */  
typedef struct { 
  int text_pos; 
  int bwt_pos;
} mc_pos;


/* ************************************************************
   *                                                          *
   * main compression routine                                 *
   *                                                          *
   ********************************************************** */
void compress_file(void)
{
  void read_text(FILE *, bwi_input *s);
  void remap_alphabet(bwi_input *s);
  void build_sa(bwi_input *s);
  void compute_bwt(bwi_input *s);
  void compute_info_superbuckets(bwi_input *s);
  void compute_info_buckets(bwi_input *s);
  void write_prologue(bwi_input *s);
  void compress_superbucket(bwi_input *s, int);
  int compute_locations(bwi_input *s);
  int compute_locations_dict(bwi_input *s, int*);
  int compute_ranks_dict(bwi_input *s, int*);
  int compute_locations_huffword(bwi_input *s, int *);
  void bit_flush( void );
  void bit_write(int,int);  
  void init_bit_buffer(void);
  void write_susp_infos(bwi_input *s);
  bwi_input s;
  int i,len, retr_occ, retr_occ2, loc_occ_range;
  int Start_prologue_ranks;

  /* --------- Load the text file from disk ------- */  
  if(Verbose) fprintf(stderr,"Reading input file... ");  
  read_text(Infile, &s);
  if(Verbose) fprintf(stderr,"done! (%f seconds)\n",getTime());
  
  /* --------- Compact alphabet ------- */  
  if(Verbose>1) fprintf(stderr,"Remapping alphabet... ");  
  remap_alphabet(&s);
  if(Verbose>1) fprintf(stderr,"done! (%f seconds). ",getTime());
  if(Verbose>1) fprintf(stderr,"Compact alphabet size = %d\n",s.alpha_size);

  /* --------- Build suffix array ------- */  
  if(Verbose) fprintf(stderr,"Building suffix array");  
  build_sa(&s);
  if(Verbose) fprintf(stderr,"done! (%f seconds)\n",getTime());

  /* --------- Compute BWT ------- */  
  if(Verbose>1) fprintf(stderr,"Computing BWT... ");
  compute_bwt(&s);
  if(Verbose>1) fprintf(stderr,"done! (%f seconds)\n",getTime());

  /* ------- mark chars and compute locations ----- */ 
  if (Is_dictionary)
    retr_occ = compute_locations_dict(&s,&loc_occ_range);    // dictionary
  else if (Is_huffword)
    retr_occ = compute_locations_huffword(&s,&loc_occ_range);// huffword 
  else if (Is_URL)
    retr_occ = compute_ranks_dict(&s,&loc_occ_range);        // URL
  else
    retr_occ = compute_locations(&s);                        // standard


  /* --------- Compute various infos for each superbucket ------- */  
  if(Verbose>1) fprintf(stderr,"Computing infos superbukets... ");
  compute_info_superbuckets(&s);
  if(Verbose>1) fprintf(stderr,"done! (%f seconds)\n", getTime());

  /* --------- Compute various infos for each bucket ------- */  
  if(Verbose>1) fprintf(stderr,"Computing infos buckets... ");
  compute_info_buckets(&s);
  if(Verbose>1) fprintf(stderr,"done! (%f seconds)\n", getTime());

  /* --------- Writing the compressed file ------- */
  Infile_size = s.text_size; 
  Outfile_size=0;

  write_prologue(&s);
  if(Verbose) fprintf(stderr,"Prologue --> %d bytes!\n",Outfile_size);

  for(i=0;i<Num_bucs_lev1;i++)
    compress_superbucket(&s,i);

  /* ---- keep starting positions of occ-explicit list ---- */
  Start_prologue_occ = Outfile_size;

  /* -- write the starting position of buckets -- */
  write_susp_infos(&s);

  if (fseek(Outfile,Start_prologue_occ,SEEK_SET)) {
    fprintf(stderr, "Seek error on output file -compress_file-\n");
    exit(1);
  }


  /* -- write the position of the marked chars ---- */
  init_bit_buffer();
  if(Is_dictionary || Is_huffword || Is_URL)
    len = int_log2(loc_occ_range);     // bits required for each rank
  else  
    len = int_log2(s.text_size);       // bits required for each pos 

  for(i=0; i < retr_occ; i++)
    bit_write(len,s.loc_occ[i]);

  bit_flush();

  Start_prologue_ranks = (int)ftell(Outfile);

  if(Verbose)  
    fprintf(stderr,"List of %d marked ranks --> %d bytes!\n",
	    retr_occ,Start_prologue_ranks-Start_prologue_occ);

  /* -- in the case of URL we also store the DICT info -- */
  /* It should be put together with the computation above --*/
  /* Thus removing these differences in the code --*/
  /* Hence Start_prologue_occ indicates the starting position of RANKS. */
  /* After retr_occ RANKS start the LOCATIONS, which are again retr_occ */
  /* in number. The value of retr_occ can be computed at decompression time */
  /* by using the same formula adopted in compute_ranks_dict() */
  if (Is_URL) {
    retr_occ2 = compute_locations_dict(&s,&loc_occ_range);  // DICT
    
    if (retr_occ != retr_occ2)
      out_of_mem("Unequal number of sampled NULLs\n");

    for(i=0; i < retr_occ; i++)
      bit_write(len,s.loc_occ[i]);
    
    bit_flush();
  
    
  if(Verbose) 
    fprintf(stderr,"List of %d marked locations --> %d bytes!\n",
	    retr_occ2,(int)ftell(Outfile) - Start_prologue_ranks);
  }
}

/* *********************************************************
   Init s->start_lev2 which is the starting position of each bucket
   in the compressed file
   ********************************************************* */             
void compute_info_buckets(bwi_input *s)
{
  void out_of_mem(char *);

  // ---- compute number of buckets -----
  assert((Bucket_size_lev1 % Bucket_size_lev2)==0);
  Num_bucs_lev2 = (s->text_size + Bucket_size_lev2 - 1) / Bucket_size_lev2;

  // ----- alloc array for buckets starting positions ---------- 
  s->start_lev2 =  (int *) malloc((Num_bucs_lev2)* sizeof(int));
  if(s->start_lev2==NULL) out_of_mem("compute_info_buckets"); 
}


/* *********************************************************
   Init s->buclist_lev1 (list of superbuckets)
   For each superbuckets init the fields:  
     occ[], alpha_size, bool_char_map[]
   ********************************************************* */             
void compute_info_superbuckets(bwi_input *s)
{
  void out_of_mem(char *);
  int i,b,temp,occ,k;
  bucket_lev1 *sb;
  
  // ---- compute number of superbuckets  -----
  Num_bucs_lev1 = (s->text_size + Bucket_size_lev1 - 1) / Bucket_size_lev1;

  // ----- alloc superbuckets ---------- 
  s->buclist_lev1 = (bucket_lev1 *)malloc(Num_bucs_lev1 * sizeof(bucket_lev1));
  if(s->buclist_lev1==NULL) out_of_mem("compute_info_superbuckets"); 

  // ------ alloc aux array for each superbucket
  for(i=0; i< Num_bucs_lev1; i++){
    sb = &s->buclist_lev1[i];   // sb points to current superbucket

    // Allocate space for data structures
    sb->occ = (int *) malloc((s->alpha_size)* sizeof(int));
    if(sb->occ==NULL) out_of_mem("compute_info_superbuckets"); 
    sb->bool_char_map = (uchar *)malloc((s->alpha_size)*sizeof(uchar));
    if(sb->bool_char_map == NULL) out_of_mem("compute_info_superbuckets"); 

    // Initialize data structures
    sb->alpha_size=0;    
    for(k=0; k<s->alpha_size; k++){
      sb->bool_char_map[k] = 0;
      sb->occ[k] = 0;
    }
  }  

  // ---- scan bwt and gather information ----------
  for(i=0;i<s->text_size;i++) {
    sb =  &s->buclist_lev1[i/Bucket_size_lev1];
    k = s->bwt[i];                           // current char
    assert(k<s->alpha_size);
    if(sb->bool_char_map[k]==0) { // build char_map of current sb
      sb->bool_char_map[k]=1; 
      sb->alpha_size++;          // compute alphabet size
    }
    sb->occ[k]++;
  }

  //----- compute occ in previous buckets ----
  for(k=0;k<s->alpha_size;k++) {
    occ=0;
    for(b=0;b<Num_bucs_lev1;b++) {   //prefix sum on OCC
      temp=s->buclist_lev1[b].occ[k];
      s->buclist_lev1[b].occ[k]=occ;
      occ += temp;
    }
  }
}


/* **************************************************************
   compress a superbucket

   NOTE the content of s->bwt is changed (remapped) !!!! 

   The number of occurrences of each char is represented using
   a special coding byte-aligned: if 7 bits are not sufficient
   then we use a further byte and indicate this by setting the 8th bit 
   of the first byte to 1. Other bytes can be used if necessary by 
   employing the same trick.
   ************************************************************** */
void compress_superbucket(bwi_input *s, int num)
{
  void write7x8(int);
  void compress_bucket(uchar *, int, int);

  bucket_lev1 sb;  
  uchar *in, c, char_map[256];
  int k, temp;
  int sb_start,sb_end,start,len,bocc[256];
  int b2,i,j;

  assert(num<Num_bucs_lev1);
  sb = s->buclist_lev1[num];       // current superbucket
  sb_start = num * Bucket_size_lev1;   // starting position of superbucket 
  sb_end = MIN(sb_start+Bucket_size_lev1,s->text_size);    
  b2 = sb_start/Bucket_size_lev2;  // initial level 2 bucket 

  temp=0;                             // build char map for superbucket 
  for(k=0;k<s->alpha_size;k++) 
    if(sb.bool_char_map[k]) 
      char_map[k]=temp++;
  assert(temp==sb.alpha_size);

  for(i=0;i< sb.alpha_size;i++)              // init occ
    bocc[i]=0;

  for(start=sb_start;start<sb_end;start+=Bucket_size_lev2,b2++) {
    s->start_lev2[b2]=Outfile_size;           // start of bucket in compr file
    len = MIN(Bucket_size_lev2,sb_end-start); // length of bucket
    in = s->bwt + start;                      // start of bucket

    if(start!=sb_start)     // if not the first bucket write occ to file
      for(k=0;k<sb.alpha_size;k++)
         write7x8(bocc[k]);

    for(j=0;j<len;j++) {               // update occ[] and remap
      assert(in[j]<s->alpha_size);
      c=char_map[in[j]];               // compute remapped char
      assert(c<sb.alpha_size);          
      in[j]=c;                         // remap
      bocc[c]++;                       // bocc is used in the next bucket
    }    

    compress_bucket(in,len,sb.alpha_size);
  }
}

/* **********************************************************************
   compress and write to file a bucket of length "len" starting at in[0].
   the compression is done as follows:
   first the charatcters are remapped (we expect only a few distinct chars
   in a single bucket)  then we use mtf (with a list of size Mtf_save)
   then we rle and compress using a unary code. 
   ********************************************************************** */ 
void compress_bucket(uchar *in, int len, int alpha_size)
{
  int int_log2(int);
  void init_bit_buffer(void);
  void bit_write(int,int);
  void bit_flush( void );
  void out_of_mem(char *);
  int mtf_string(uchar *, uchar *, uchar *, int);
  void rle_hierarchical(uchar *, int, int);  
  void multihuf_compr(uchar *, int, int);
  int k,j,bits_x_char,local_alpha_size,mtf_len;
  uchar c,mtf[256],local_bool_map[256], local_map[256]; 
  uchar *mtf_seq;

  /* ---------- init ------------ */
  init_bit_buffer();

  /* ---------- compute and write local boolean map ------ */
  for(k=0;k<alpha_size;k++)     
    local_bool_map[k]=local_map[k]=0;
  local_alpha_size=0;

  for(j=0;j<len;j++) {             // compute local boolean map
    c=in[j];                       // remapped char
    assert(c<alpha_size);                              
    local_bool_map[c]=1;     
  }

  for(k=0;k<alpha_size;k++)      // compute local map
    if(local_bool_map[k])
      local_map[k]=local_alpha_size++;  

  for(j=0;j<len;j++)             // remap bucket
    in[j]=local_map[in[j]];       
  
  for(k=0;k<alpha_size;k++)     // write bool char map to file 
    if(local_bool_map[k]) bit_write(1,1); 
    else bit_write(1,0);

  /* ----------- compute and write mtf picture ------------- */
  mtf_seq = (uchar *) malloc(2*len*sizeof(uchar)); // mtf temporary buffer
  if(mtf_seq==NULL) out_of_mem("compress_bucket (mtf_seq)");
  mtf_len = mtf_string(in,mtf_seq,mtf,len);  // mtf_seq=mtf(in), init mtf-list
  bits_x_char = int_log2(local_alpha_size);   // write mtf to file
  for(k=0;k<MIN(Mtf_save,local_alpha_size);k++) {
    bit_write(bits_x_char,mtf[k]);  
  }


  // -- Applies the proper compression routine --
  switch (Type_compression) 
    {
    case ARITH:  // ---- Arithmetic compression of the bucket -----
      fatal_error("Arithmetic coding no longer available -compress_bucket-\n");
      exit(1);
    case HIER3:  // ---- three-leveled model: Fenwick's proposal -----
      rle_hierarchical(mtf_seq, mtf_len,local_alpha_size);
      break;
    case UNARY:  // ---- Unary compression of mtf-ranks with escape -----
      fatal_error("Unary coding no longer available -compress_bucket-\n");
      exit(1);
    case MULTIH: // ---- RLE + MultiHuffman compression of the bucket -----
      multihuf_compr(mtf_seq,mtf_len,local_alpha_size);  
      break;
    default:
      fprintf(stderr,"\n Compression algorithm unknown! ");
      fprintf(stderr,"-compress_superbucket-\n");
      exit(1);
    }
  bit_flush();         // Byte-align the next compressed bucket
  free(mtf_seq);
}



/* ********************************************************************
   compute mtf for a string
   input
     int   len         size of input string
     uchar *in         input string of size len
     uchar *out        empty string of size 2*len
     uchar *mtf_start  empty string of size Mtf_save
   output
     write the compressed string in *out 
     return the size of the compressed string
     return in mtf_start the starting picture of the mtf list
   note
     this is a plain mtf encoder with a mtf list of size Mtf_save. If
     a char C is not in the mtf list it is encoded with the pair 
     (Mtf_save, C). The only unusual point is that the initial status
     of the mtf list is determined during the encoding: obviously we
     choose the initial status wich produces as few as possible escapes.
     Note that if fewer than Mtf_save distinct chars appear in the input
     the mtf list is not completely determined, but this is ok since
     it does not affect decompression. 
   ****************************************************************** */
int mtf_string(uchar *in, uchar *out, uchar *mtf_start, int len)
{
  int i,j,o,m,h;
  uchar c, mtf[256];

  m=0;   // # of chars in the mtf list
  o=0;   // # of char in the output string

  for(i=0;i<len;i++) {
    c=in[i];
    // --- search c in mtf ----
    for(h=0;h<m;h++)
       if(mtf[h]==c) break;

    if(h<m) {
      // -------- c found in mtf[h] -----------
      out[o++]= (uchar) h;   
      for(j=h;j>0;j--)            // update mtf[]
         mtf[j]=mtf[j-1];
      mtf[0]=c;
    }
    else {
      // -------- c was not in mtf[] -------
      if(m<Mtf_save) {      // initial mtf list not completely full
        mtf_start[m]=c;     // write c at position m of the initial mtf list
        out[o++]= (uchar)m; // the mtf value is m
        m++;                // the list has grown by one
      }
      else {                // initial mtf list is already full
        assert(m==Mtf_save); 
        out[o++]= (uchar) Mtf_save; // escape code
        out[o++]=c;        // remapped char
      }
      for(j=m-1;j>0;j--)    // update mtf list
        mtf[j]=mtf[j-1];
      mtf[0]=c;
    }
  }
  assert(o<=2*len);
  return o;
}


/* *********************************************************
   The current format of the prologue is the following:
         8 bits      type of compression (2=Hier, 4=Multi Table Huff)
         1 int       size of input file
         1 int       position of eof in s->bw
         1 uint16    size of a super bucket divided by 1024
         1 uchar     size of a bucket divided by 1024 (divides the previous)
	 1 uchar     size-1 mtf_list stored in each bucket  
	 1 uchar     size-1 of the compacted alphabet of the text
	 1 uchar     remapped char selected for occurrence list
	 1 int       # skipped occ of chosen_char in bwt
	 1 int       starting byte of occ-explicit list
       256 bits      boolean map of chars in the text (S = # of 1)
         S int       prefix sum of character occurrences

      for each superbucket
	 S' bytes    map of compact_alph chars occurring in THIS superbucket
	             (S' = (S+7)/8 -- it is byte aligned)  ****FLUSH****
         S int       # occ of all compact_alphabet chars in prev superbuckets

      finally:
        NB x L    starting position of each bucket in the compressed file
                  NB is the number of buckets and L is the number of bits
		  sufficient to represent that length (byte_aligned)


     -------- Body of the compressed file [byte-aligned]   -------------

     for each bucket (let R be the # of distinct chars in the superbucket)

         R 7x8val   # of occ of each char in the previous buckets of the
                    same superbucket. Each value is represented with 
                    the 7x8 encoding. This information is missing for the 
                    first bucket of each superbucket

         R  bits    map of chars appearing in this bucket. 
                    Let R' be the # of distinct chars in this bucket 
                    and L' the # of bits required to represent R'

       L' x M bits  Initial move to front list for this bucket
                    M = min(R',Mtf_save)         

	 ... bits   needed to byte-align in case ONLY of Arith-coding	    

         ??? bits   compressed bucket in mtf + rle + [Ari|Hier|Una] format 

         --- bits   ****FLUSH**** to have byte alignment

     -------- Body of the occ explicit list [byte-aligned]   -------------
     ---------------------------------------------------------------------
     --- URL: we have occ for text positions and rows ---

         ... L bits  list of positions where character ch occurs in
	             the original text. ch = character that occurs
		     close to Marked_char_freq times in the text.
  **************************************************************** */
void write_prologue(bwi_input *s)
{
  void init_bit_buffer(void);
  int int_log2(int);
  void uint_write(int);
  void bit_write(int,int);
  void bit_flush(void);
  void write7x8(int);
  bucket_lev1 sb;
  int i,len,k;

  /* ----- write file and bucket size ------ */
  init_bit_buffer();
  bit_write(8,Type_compression);
  uint_write(s->text_size);
  uint_write(s->bwt_eof_pos);

  assert(Bucket_size_lev1>>10<65536);
  assert((Bucket_size_lev1 & 0x3ff) == 0);
  bit_write(16,Bucket_size_lev1>>10);

  assert(Bucket_size_lev2>>10<256);
  assert((Bucket_size_lev2 & 0x3ff) == 0);
  bit_write(8,Bucket_size_lev2>>10);

  // ---- mtf and alphabet information
  assert(Mtf_save>0 && Mtf_save<=256);
  bit_write(8,Mtf_save-1);

  assert(s->alpha_size>0 && s->alpha_size<=256);
  bit_write(8,s->alpha_size-1);   

  // ---- write chosen_char & starting byte of occ-list
  bit_write(8,s->chosen_char);
  uint_write(s->skip);
  uint_write(0);
  
  // ---- boolean alphabet char map
  for(i=0;i<256;i++)
   if(s->bool_char_map[i]) bit_write(1,1);
   else bit_write(1,0);  

  // ---- write prefix sum of char occ
  for(i=0; i<s->alpha_size; i++)
    uint_write(s->pfx_char_occ[i]);

  // ----- process superbuckets
  for(i=0;i<Num_bucs_lev1;i++) {
    sb = s->buclist_lev1[i];

    for(k=0;k<s->alpha_size;k++)     // boolean char_map
      if(sb.bool_char_map[k]) bit_write(1,1);
      else bit_write(1,0);
    bit_flush();                    // we keep everything byte aligned

    if(i>0)                          // write prefix-occ 
      for(k=0;k<s->alpha_size;k++)
        uint_write(sb.occ[k]);
  }

  // ----- leave space for storing the start positions of buckets
  len = (int_log2(s->text_size)+7)/8;   //it's byte-aligned
  for(i=0;i<Num_bucs_lev2;i++)
    bit_write(len * 8,0);

}

/* *****************************************************************
   write the starting position (in the output file) of each one 
   of the Num_bucs_lev2 buckets. For simplicity we use 32 bits for
   each position. These values are written at the end of the prologue
   just before the beginning of the first bucket.
   It writes also the starting position of the occurrence list
   ***************************************************************** */ 
void write_susp_infos(bwi_input *s)
{
  void bit_write(int,int);
  void uint_write(int);
  void bit_flush(void);
  int int_log2(int);
  int i,offset,len;

  /* -- write starting position of occ-explicit list --*/
  // warning! the constant 19 depends on the structure of the prologue!!!
  if (fseek(Outfile,19,SEEK_SET)) {
    fprintf(stderr,"seek error on output file -write_susp_infos-\n");
    exit(1);
  }
  uint_write(Start_prologue_occ);
  bit_flush();

  // Warning: the offset heavily depends on the structure of prologue.
  //          The value of start_level2[0] has been initialized in
  //          the procedure compress_superbucket()      
  len = (int_log2(s->text_size)+7)/8;   // variable length representation
  offset = s->start_lev2[0] - Num_bucs_lev2*len;
  if (fseek(Outfile,offset,SEEK_SET)) {
    fprintf(stderr,"seek error on output file -write_susp_infos-\n");
    exit(1);
  }

  for(i=0;i<Num_bucs_lev2;i++)
    bit_write(len*8,s->start_lev2[i]);
  bit_flush();

  assert(ftell(Outfile)==(int)s->start_lev2[0]);
}
 

/* ****************************************************************
   Read input text and mark which are the characters actually used 
   initialize the variables:
      s->text[], s->text_size, s->bool_char_map[]
   **************************************************************** */
void read_text(FILE *f, bwi_input *s)
{
  void out_of_mem(char *);
  int i;

  /* --------- get file size ------- */
  if (fseek(f, 0L, SEEK_END)) {
    fprintf(stderr, "Seek error on input file -read_text-\n");
    exit(1);
  }
  s->text_size=ftell(f);
  if (s->text_size==0) {
    fprintf(stderr, "Input file file empty -read_text-\n");
    exit(1);
  }
  
  /* --- alloc memory for text (the overshoot is for suffix sorting --- */
  s->text = (char *) malloc(s->text_size+N_OVERSHOOT);
  if(s->text==NULL)
    out_of_mem("read_input");
    
  for(i=0;i<256;i++)     /* clear char_map */
    s->bool_char_map[i]=0;

#if 1
  {
    ssize_t t;
    
    // ---- read text in one sweep -----
    lseek(fileno(f),0,SEEK_SET);
    t=read(fileno(f),s->text, (size_t) s->text_size);
    if(t!=s->text_size) {
      fprintf(stderr,"Error in read() -read_input-\n");
      perror("(read_input)");
      exit(1);
    }
    // --- init boolean char map ----
    for(i=0;i<s->text_size;i++)
      s->bool_char_map[s->text[i]]=1;
  }
#else
  {
    int c;
  /* ----- read text and init s->bool_char_map --- */
  for(rewind(f), i=0;i<s->text_size;i++) {
    c = getc(f); 
    if(c==EOF) {
      fprintf(stderr,"EOF Unexpected! -read_input-\n");
      exit(1);
    }
    assert(c>=0 && c<256);
    s->text[i] = c;
    s->bool_char_map[c]=1;
  }
  c=getc(f);
  assert(c==EOF);
  }
#endif
}


/* ****************************************************************
   compute the size of the alphabet and change the text so
   that every char is in the range 0 .. size -1.
   The correspondence between old and new alphabet is 
   stored in s->char_map[]. Also initialize s->pfx_char_occ[] so that 
   in entry i contains the # of occ's of remapped chars 0 ... i-1
   **************************************************************** */
void remap_alphabet(bwi_input *s)
{
  int i,sum,temp;

  s->alpha_size=0;
  for(i=0;i<256;i++) {
    s->pfx_char_occ[i]=0;
    if(s->bool_char_map[i]) 
      s->char_map[i]=s->alpha_size++;
  }

  for(i=0;i<s->text_size;i++){
    s->text[i] = s->char_map[s->text[i]];   // Remap the text
    s->pfx_char_occ[s->text[i]]++;    //Compute character occurrences
  }

  sum=0;
  for(i=0;i<s->alpha_size;i++){
    temp = s->pfx_char_occ[i];
    s->pfx_char_occ[i] = sum;    //Compute prefix sum  of char occ
    sum += temp;
  }
}     


/* **********************************************************
   build suffix array for the text in s->text[].
   Input:
     s->text, s->text_size
   Ouput:
     s->sa
   ********************************************************** */
extern int Use_larsson_sada;    // defined in bwi.c
extern char *Safile_name;       // defined in bwi.c 
void build_sa(bwi_input *s)
{
  int scmp3(unsigned char *p, unsigned char *q, int maxl);
  void init_bit_buffer(void);
  int fbit_read(FILE *,int);
  int *larsson_sada_sufsort(uchar *, int, int);
  int *suffixsort5n(uchar *, int);
  void out_of_mem(char *s);
  int int_log2(int);  
  int i, n, pointer_size,q,r,sa_size;
  FILE *safile;
  
  /* ------------ check sa file ---------------- */
  n=0;
  safile = fopen(Safile_name,"rb");
  if(safile!=NULL) {
    fseek(safile,0L,SEEK_END);
    n=ftell(safile);
  }

  if (n==0) { 
    // ------- build sa using larsson-sada or 5n
    if(Verbose)  fprintf(stderr, " from scratch ");
    if(Use_larsson_sada) {
      if(Verbose)  fprintf(stderr, "(using ls) ... ");
      s->sa = larsson_sada_sufsort(s->text,s->text_size,s->alpha_size);
    }
    else {
      if(Verbose)  fprintf(stderr, "(using 5n) ... ");
      s->sa = suffixsort5n(s->text,s->text_size);
    }
  } 
  else {     
    // ------ read sa from file --------     
    pointer_size = int_log2(s->text_size);
    // --- compute  sa_size = s->text_size * pointer_size + 7)/8
    // --- use q and r to avoid overflow
    q = s->text_size/8; r = s->text_size % 8;
    sa_size = (q*pointer_size) + (r*pointer_size+7)/8; 
    if (n != sa_size)
      fatal_error("Invalid .sa file\n");
    if(Verbose) fprintf(stderr, " by reading it from file... ");
    // allocate space for the suffix array
    s->sa = (int *) malloc(s->text_size * sizeof(int));
    if(s->sa==NULL) out_of_mem("build_sa");
    rewind(safile);
    init_bit_buffer();
    for(i=0; i<s->text_size; i++)// read one suffix-array pointer at a time
      s->sa[i] = fbit_read(safile,pointer_size);
    fclose(safile);
  }
  // check the suffix array
#if 0
   for (i=0; i<s->text_size-1; ++i)
     if (scmp3(s->text+s->sa[i], s->text+s->sa[i+1], 
                MIN(s->text_size-s->sa[i], s->text_size-s->sa[i+1]))>=0) {
       fprintf(stderr, "Suffix array check failed at position %d\n", i);
       exit(1);
     }
#endif
}

/* ***************************************************************
   build the suffix array calling the larrson-sadakane algorithm
   ************************************************************** */
int *larsson_sada_sufsort(uchar * s, int size, int alpha_size)
{
  void out_of_mem(char *s);
  void suffixsort(int *x, int *p, int n, int k, int l);
  int *sa, *aux, i;

  sa = (int *) malloc((1+size)*sizeof(int));
  aux = (int *) malloc((1+size)*sizeof(int));
  if(!sa || !aux)
    out_of_mem("larsson_sada_sufsort");

  /* ---- copy text in auxiliary array --------- */
  for(i=0;i<size;i++)
    aux[i] = s[i];   

  /* ----- build sa ---- */
  suffixsort(aux,sa,size,alpha_size,0);
  free(aux);
  return sa+1;      /* discard first position */
}


/* *********************************************************
   compute BWT of input text.
   Input:
     s->text, s->text_size, s->sa
   Output
     s->bwt, s->bwt_eof_pos
   ********************************************************* */               
void compute_bwt(bwi_input *s)
{
  void out_of_mem(char *);
  int i,j;
  
  /* ------ alloc memory for bwt -------- */
  s->bwt = (char *) malloc(s->text_size);
  if(s->bwt==NULL)
    out_of_mem("compute_bwt");
      
  s->bwt[0] = s->text[s->text_size-1];  //L[0] = Text[n-1]
  for(i=0,j=1;i<s->text_size;i++) 
    if(s->sa[i]!=0)
      s->bwt[j++]=s->text[s->sa[i]-1];
    else
      s->bwt_eof_pos = i; // EOF is not written but its position remembered !
      
  assert(j==s->text_size);
}       


/* *****************************************************************
   compute locations of "marked" occurrences. This procedure does
   the following: 
     1) compute the desired # of marked chars (=desired_marked_chars)
     2) compute the best pair i,j such that (occ[i]/2^j) is as close
        as possible (but <= ) to desired_marked_chars
        write  i in s->chosen_char and 2^j in s->skip 
     3) scan s->bwt[] and "select" one out of s->skip occurrences of 
        s->chosen_char. For each selected occurrence write in s->loc_occ
       its position in the original text.

   I think this procedure could be improved (that is, simplified 
   and faster in doing the search with bwhuffw) using the ideas
   introduced in compute_locations_dict() and compute_locations_huffword()
   more precisely
    1) remove s->chosen_char and mark simply one row every s->skip
       (this would make the marked chars more evely distributed in the
        text).
    2) consider the row starting with a "marked char" rather than ending
       (this would simplify the code)
   **************************************************************** */
int compute_locations(bwi_input *s)
{
  int i,max,j,count,chosen_occ;
  int ch_occ[256],rescaled,skip;
  int exponent, desired_marked_chars, marked_chars;
  
  if(Marked_char_freq==0) {
    s->skip=0; s->chosen_char = 0; 
    return 0;
  }
  /* ------ compute the desired number of marked chars ------ */
  desired_marked_chars =  (int) (s->text_size * Marked_char_freq);
  if(desired_marked_chars==0)
    desired_marked_chars=1;
      
  // ---- Count occurrences for each character
  ch_occ[s->alpha_size-1]= s->text_size-s->pfx_char_occ[s->alpha_size-1];
  for(i=0;i<s->alpha_size-1;i++)
    ch_occ[i]=s->pfx_char_occ[i+1]-s->pfx_char_occ[i];
  
  // ----- select best (char,skip) pair
  for(i=0, max=-1; i<s->alpha_size; i++){
    if(i==s->bwt[0]) continue;  // Exclude bwt-first-char (see below)
    /* --- determine the number of skipped char for i */
    if (ch_occ[i] > desired_marked_chars) {
      exponent = int_log2(ch_occ[i]/desired_marked_chars);
      assert(exponent > 0);
      skip = int_pow2(exponent);
    }
    else
      skip = 1;
    /* --- check if this is the best choice seen so far --- */
    rescaled = ch_occ[i] / skip;
    if(rescaled>max && rescaled <= desired_marked_chars) {
      max = rescaled;
      s->chosen_char = i;
      s->skip = skip;
    }
  }
  assert(max > 0);
  assert(s->skip>0);
  
  if(Verbose>1) {
    for(i=0;i<256;i++)
      if(s->char_map[i]==s->chosen_char) break;
    fprintf(stderr,"Marked char is ascii %d; ", i); 
    fprintf(stderr,"one occ every %d is marked; ",s->skip);
  }

  // ------- compute number of marked chars
  chosen_occ = ch_occ[s->chosen_char];
  if(chosen_occ % s->skip)
    marked_chars = chosen_occ/s->skip + 1;
  else
    marked_chars = chosen_occ/s->skip;
  // -------- alloc s->loc_occ
  s->loc_occ = (int *) malloc(sizeof(int) * (marked_chars));

  // write the text location of the ROWS ending with ch
  for(i=1,j=0,count=0; i<s->text_size; i++) 
    {                         // bwt[0] is not the marked char (see above) 
      if (s->bwt[i] == s->chosen_char) {
	if ((count % s->skip) == 0) {
	  if (i <= s->bwt_eof_pos) { 
	    s->loc_occ[j] = (int) s->sa[i-1];
	    assert(s->text[s->loc_occ[j]-1] == s->chosen_char);
	  }
	  else {
	    s->loc_occ[j] = s->sa[i]; 
	    assert(s->text[s->loc_occ[j]-1] == s->chosen_char);
	  }
	  j++;
	} 
	count++;
      }
    }
  // j is the number of marked chars
  if(Verbose>1) 
    fprintf(stderr,"%d chars marked.\n", j); 
  assert(j == marked_chars);
  assert(count == chosen_occ);  
  return j;
}       



/* ********************************************************************
   compute locations of the word in a dictionary this is the marking 
   strategy specialized for the case in which we are compressing a dictionary. 
   The marked char is always the NULL (\0) char which
   is used in the dictionary as a word separator; therefore
   the procedure simply determines the optimal skip.
   Note that here we are interested in words rather than in chars, 
   so we write in s->loc_occ the number of words preceeding a given
   one (rather than the number of chars).
   For this reason we return in *rank_range the # of distict ranges; that is,
   the values written in s->loc_occ are in the range [0,*rank_range)

   Note that there are 2 importanto point to be noted:
   1) here the marked rows are those BEGINNING with NULL (in 
      compute_locations() we use the rows ending with s->chosen_char)
   2) for each marked row (that is row beginning with NULL and multiple
      of s->skip) we maintain the position in the text (=dictionary) of the 
      char in the FIRST column.

   Note also that the percentage Marked_char_freq here refer to the 
   total number of words, rather than to the dictionary size in chars.
   ******************************************************************** */
int compute_locations_dict(bwi_input *s, int *rank_range)
{
  int mc_pos_cmp(const void *, const void *);
  int i,j,null_occ,count;
  int desired, marked_chars;
  uchar null_remap;
  mc_pos *aux;

  // -------- easy case: no marked char
  if(Marked_char_freq==0) {
    s->skip=0; s->chosen_char = 0; *rank_range=0;
    return 0;
  }

  // -------- Make sure that NULL occurs
  if(s->bool_char_map[0]==0) 
    fatal_error("The dictionary does not contain the NULL char!\n");
  // -------- Determine the remap for NULL
  null_remap = s->char_map[0];
  s->chosen_char = null_remap;
  // -------- Determine the # of occ of null_remap
  if(null_remap==s->alpha_size-1)
    null_occ=s->text_size - s->pfx_char_occ[null_remap];
  else
    null_occ=s->pfx_char_occ[null_remap+1] - s->pfx_char_occ[null_remap];

  // -------- determine how many occ to skip
  desired = (int) (Marked_char_freq*null_occ);
  if(desired==0) desired=1;
  // compute s->skip
  if(null_occ>desired) {
    // s->skip= int_pow2(int_log2(null_occ/desired)); do we need a pow of 2?
    s->skip = 1+null_occ/(desired+1);
    assert(s->skip > 0);
  }
  else
    s->skip=1;

  // -------- compute number of marked chars
  if(null_occ % s->skip)
    marked_chars = null_occ/s->skip + 1;
  else
    marked_chars = null_occ/s->skip;

  // --- alloc s->loc_occ and auxiliary struct for text positions
  s->loc_occ = (int *) malloc(sizeof(int) * (marked_chars));
  aux = (mc_pos *) malloc(marked_chars*sizeof(mc_pos));
  if(aux==NULL || s->loc_occ==NULL) 
    out_of_mem("compute_locations_dict (aux/loc_occ)");

  // -----  determine the text position of the marked chars
  for(j=i=0; i<null_occ; i++)  {                          
    if ((i % s->skip) == 0) {
      aux[j].text_pos = s->sa[s->pfx_char_occ[null_remap] + i]; 
      assert(aux[j].text_pos>=0 && aux[j].text_pos<s->text_size);
      assert(s->text[aux[j].text_pos] == null_remap);
      aux[j].bwt_pos = j;
      j++;                 // positions written so far
    }
  }
  assert(j==marked_chars);

  // ---- sort marked chars according to their text position 
  qsort(aux, marked_chars, sizeof(mc_pos), mc_pos_cmp);

  // ---- Assign ranks to rows ----
  for(i=j=count=0;i<s->text_size;i++)
    if(s->text[i]==null_remap) {
      if(i==aux[j].text_pos)  {                // was this position marked?
	s->loc_occ[aux[j].bwt_pos] = count;    // yes! write "word number"
	if(++j==marked_chars)                  // advance to next marked pos
	  break;                               // no more marked pos   
      }
      count++;                                 // increase "word number"
    }
  assert(count<=null_occ);
  assert(j==marked_chars);
  free(aux);

  *rank_range = null_occ;   // ranks are in the range [0,null_occ)
  return marked_chars;
}


/* ********************************************************************
   PAOLO: Per il calcolo dei rank delle parole, cosi da poter
   realizzare il mapping inverso. 
   ******************************************************************** */
int compute_ranks_dict(bwi_input *s, int *rank_range)
{
  int mc_pos_cmp(const void *, const void *);
  int i,j,null_occ,count;
  int desired, marked_chars;
  uchar null_remap;
  mc_pos *aux;

  // -------- easy case: no marked char
  if(Marked_char_freq==0) {
    s->skip=0; s->chosen_char = 0; *rank_range=0;
    return 0;
  }

  // -------- Make sure that NULL occurs
  if(s->bool_char_map[0]==0) 
    fatal_error("The dictionary does not contain the NULL char!\n");
  // -------- Determine the remap for NULL
  null_remap = s->char_map[0];
  s->chosen_char = null_remap;
  // -------- Determine the # of occ of null_remap
  if(null_remap==s->alpha_size-1)
    null_occ=s->text_size - s->pfx_char_occ[null_remap];
  else
    null_occ=s->pfx_char_occ[null_remap+1] - s->pfx_char_occ[null_remap];


  // -------- determine how many occ to skip
  desired = (int) (Marked_char_freq*null_occ);
  if(desired==0) desired=1;

  // compute s->skip
  if(null_occ>desired) {
    // s->skip= int_pow2(int_log2(null_occ/desired)); do we need a pow of 2?
    s->skip = 1+null_occ/(desired+1);
    assert(s->skip > 0);
  }
  else
    s->skip=1;

  // -------- compute number of marked chars
  if(null_occ % s->skip)
    marked_chars = null_occ/s->skip + 1;
  else
    marked_chars = null_occ/s->skip;

  // --- alloc s->loc_occ and auxiliary struct for text positions
  s->loc_occ = (int *) malloc(sizeof(int) * (marked_chars));
  aux = (mc_pos *) malloc(null_occ*sizeof(mc_pos));
  if(aux==NULL || s->loc_occ==NULL) 
    fatal_error("compute_locations_dict (aux/loc_occ)");

  // -----  determine the text position of all the \0 chars
  for(i=0; i<null_occ; i++)  {                          

    // SA does not have the suffix starting with EOF, but
    // the bwt matrix has this as the first row
    aux[i].text_pos = s->sa[s->pfx_char_occ[null_remap] + i]; 
    
    assert(aux[i].text_pos>=0 && aux[i].text_pos<s->text_size);
    assert(s->text[aux[i].text_pos] == null_remap);
    
    // We deal with the fact that the first row is the one starting
    // with EOF but actually is not represented. Hence we sum 1
    aux[i].bwt_pos = EOF_shift(s->pfx_char_occ[null_remap] + i);

  }

  // ---- sort all \0 according to their text position 
  qsort(aux, null_occ, sizeof(mc_pos), mc_pos_cmp);

  // ---- Assign ranks to rows ----
  for(i=count=0,j=1;i<s->text_size;i++)
    if(s->text[i]==null_remap) {

      if (j==null_occ){   // last NULL char is surely marked
	s->loc_occ[count] = aux[j-1].bwt_pos;  
	count++;   // how many marked
      } else if(j % s->skip == 0)  {        // this position must be saved
	s->loc_occ[count] = aux[j-1].bwt_pos;    // write the "word rank"
	count++;   // how many marked
      }
	
      j++;                                  // increase the "word rank"
    }

  free(aux);

  if(Verbose>1) {
    fprintf(stderr,"Rank: Marked %d NULLs out of %d ", count,null_occ); 
    fprintf(stderr,"(one every %d is marked)\n", s->skip);
  }

  if(count != marked_chars)
    fatal_error("No correct marking in compute_ranks_dict!\n");

  *rank_range = null_occ;   // ranks are in the range [0,null_occ)
  return count;
}



/* *****************************************************************
   This procedure works as compute_locations_() but is specialized for
   the case in which we are processing the output of the 7-bit huffword 
   compressor (recall that the 8th bit is set to 1 for the byte at 
   the beginning of a word). 
   This procedure has been derived from compute_locations_dict()
   so look there for comments. The only difference is that instead of
   consider for marking the rows starting with NULL here we consider
   the rows starting with a char >=128 (>=bit8_set after remapping).
   We mark a fraction (Marked_char_freq) of these rows. Note that there
   is no longer the concept of chosen char (s->chosen_char). We think this
   is more robust and ensure a better "equipartion" of the marked char
   in the text (recall the problem with jdk13c in the FM-index).

   For each marked char we store the number of words preceding it.
   For this reason we return in *word_range the # of distinct words; that is,
   the values written in s->loc_occ are in the range [0,*word_range)

   Note that the percentage Marked_char_freq here refer to the 
   number of words, rather than to the total size of the body in chars.
   **************************************************************** */
int compute_locations_huffword(bwi_input *s, int *word_range)
{
  int mc_pos_cmp(const void *, const void *);
  int i,j,count,words_occ, bit8_set, marked_chars, desired;
  mc_pos *aux;

  // ---- easy case: no marked chars
  if(Marked_char_freq==0) {
    s->skip=0; s->chosen_char = 0; *word_range=0;
    return 0;
  }

  // --- no single char is marked here! 
  // --- we "mark" the chars >= 128 (the beginning of a word)
  s->chosen_char = 255;

  /* ---- compute codeword of the first char>=128 ------- */
  for(i=0,bit8_set=0;i<128;i++) 
    if(s->bool_char_map[i]) bit8_set++;
  if(bit8_set==s->alpha_size) 
    fatal_error("Invalid body, no chars with the 8th bit set!\n");
  else if(bit8_set>s->alpha_size)
    fatal_error("Invalid s->bool_char_map[]!\n");

  // ---- get total number of words ----------
  words_occ = s->text_size - s->pfx_char_occ[bit8_set];

  // -------- determine how many occ to skip
  desired = (int) (Marked_char_freq*words_occ);
  if(desired==0) desired=1;
  // --------  compute s->skip
  if(words_occ>desired) {
    s->skip = 1+words_occ/(desired+1);
    assert(s->skip > 0);
  }
  else
    s->skip=1;

  // -------- compute number of marked chars
  if(words_occ % s->skip)
    marked_chars = words_occ/s->skip + 1;
  else
    marked_chars = words_occ/s->skip;

  // --- alloc s->loc_occ and auxiliary struct for text positions
  s->loc_occ = (int *) malloc(sizeof(int) * (marked_chars));
  aux = (mc_pos *) malloc(marked_chars*sizeof(mc_pos));
  if(aux==NULL || s->loc_occ==NULL) 
    out_of_mem("compute_locations_huffword (aux/loc_occ)");

  // -----  determine the text position of the marked chars
  for(j=i=0; i<words_occ; i++)  {                          
    if ((i % s->skip) == 0) {
      aux[j].text_pos = s->sa[s->pfx_char_occ[bit8_set] + i]; 
      assert(aux[j].text_pos>=0 && aux[j].text_pos<s->text_size);
      assert(s->text[aux[j].text_pos] >= bit8_set);
      aux[j].bwt_pos = j;
      j++;                 // positions written so far
    }
  }
  assert(j==marked_chars);

  // ---- sort marked chars according to their text position 
  qsort(aux, marked_chars, sizeof(mc_pos), mc_pos_cmp);
  for(i=j=count=0;i<s->text_size;i++)
    if(s->text[i]>=bit8_set) {
      if(i==aux[j].text_pos)  {                // was this position marked?
	s->loc_occ[aux[j].bwt_pos] = count;    // yes! write "word number"
	if(++j==marked_chars)                  // advance to next marked pos
	  break;                               // no more marked pos   
      }
      count++;                                 // increase "word number"
    }
  assert(count<=words_occ);
  assert(j==marked_chars);
  free(aux);

  if(Verbose>1) {
    fprintf(stderr,"Marked %d BegOfWord out of %d ", marked_chars,words_occ); 
    fprintf(stderr,"(one every %d is marked)\n", s->skip);
  }
  *word_range = words_occ;   // words are in the range [0,words_occ)
  return marked_chars;
}


/* ---- comparison function to sort occ in order of apperance in T --- */
int mc_pos_cmp(const void *a, const void *b)
{
  mc_pos *aa,*bb;
  
  aa = (mc_pos *) a;
  bb = (mc_pos *) b;
  return ( aa->text_pos -  bb->text_pos); 
}


/* ********************************************************************
   rle+compression of a string: Fenwick;s three-level model 
   input
     int len         size of mtf sequence
     uchar *in         input mtf sequence
     int alpha_size   size of the alphabet
   output
     the compressed string is written in the output file 
   note 
     the mtf rank is coded as follows: if <=2 then the first level
     codes the rank using 2 bits. In case (>2 and <=9) then we 
     write 11 as escape code and go to the second level where the rank
     is coded using 3 bits as (rank - 3). If rank > 9, then again an 
     escape code is output as 111, and (rank - 10) is represented
     using enough bits. One important point is that, when the rank
     is equal to Mtf_save then the next value denotes a character and
     thus it is coded in binary using a proper number of bits.
   ******************************************************************* */
void rle_hierarchical(uchar *in, int len, int alpha_size)
{
  int int_log2(int);
  void bit_write(int,int);
  int bits_x_char,i,z;
  uchar c;
  int mtf_size;

  mtf_size = MIN(Mtf_save,alpha_size);
  bits_x_char = int_log2(alpha_size);
  z=-1;                       // # of pending zeroes (-1) 
  for(i=0; i<len;i++) {
    assert(in[i]<alpha_size);

    if(in[i]==0) {
	 
      if(++z==255) {
	 bit_write(2,0);     // unary code for 0 
         bit_write(8,255);   // write 255 using 8 bits
	 z=-1;
      }
    }
    else {
      /* ----- check if there are pending zeores ---- */
      if(z>=0) {
	bit_write(2,0);     // unary code for 0 
	bit_write(8,z);     // write z using 8 bits
        z=-1;
      }
      // ---- write a nonzero mtf rank ----- 
      if(in[i]<=2) {
	bit_write(2,in[i]);   // binary coding
      }
      if((in[i]>2) && (in[i] <= 9)){
	bit_write(2,3);   // escape level 1
	bit_write(3,in[i]-3);   // binary coding of second level
      }
      if((in[i] > 9) && (in[i] <= mtf_size)){
	bit_write(5,31);           // escape code level 2
	bit_write(int_log2(mtf_size-10),in[i]-10);
      }
      if(in[i] == mtf_size){
	c=in[++i];                      // get actual char
        assert(c<alpha_size); 
	bit_write(bits_x_char,c);      // write remapped char
      }
    }
  }
  // ---- there could be some pending zeroes
  if(z>=0) {
     bit_write(2,0);     // unary code for 0 
     bit_write(8,z);     // 255 using 8 bits
  }
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









