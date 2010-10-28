/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   bwt-based compression and indexing
   Routines for run-length encoding and compression:
   compr_routines.c  -----  P. Ferragina & G. Manzini, 10 June 2000

      rle_hierarchical();
   
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "common.h"



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












































