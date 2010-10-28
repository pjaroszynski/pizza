
#include "common.h"


extern double Cache_percentage;  // size of cache (% wrt uncompressed size)
extern int Report_occ;       // if !=0 report the postions of the occurrences   
extern int Report_position;       // if !=0 report the postions of the occurrences   
extern int Locate_occ;       // if !=0 compute the postions of the occurrences   
extern int Display_occ;      // if !=0 display the text sourronding each occurrence
extern int Oneline_report;   // report time occ and startrow in a single line


int fm_search(char *pattern, int pat_len, char *infile_name)
{
  void fatal_error(char *);
  void read_prologue(bwi_out *);

  bwi_out s_main;
  bwi_out *s = &s_main; 
  int num_occ,sp,ep;

  Report_occ = 0;
  Locate_occ = 0;
  Display_occ = 0;
  Oneline_report = 0;
  Type_mem_ops = EXT_MEM;
  Cache_percentage = 0;
 
  if (check_bwi_suffix(infile_name) == 0)
    fatal_error("The file name must end with .bwi -main-\n");  

  my_open_file(infile_name);

  init_bwi_cache();

  read_prologue(s); 

  num_occ = bwsearch(s, pattern, pat_len, &sp, &ep);

  my_fclose(Infile);

  return num_occ;
}



int fm_getid(char *url,char *infile_name)
{
  int url_id;
  void get_id_byurl(char *url, int *Url_id);
 

  Type_mem_ops = EXT_MEM;
  Cache_percentage = 0;
  Report_position=0;
  
  if (check_bwi_suffix(infile_name) == 0)
    fatal_error("The file name must end with .bwi -main-\n");  
  
  my_open_file(infile_name);
  
  init_bwi_cache();
  
  get_id_byurl(url, &url_id);
  
  my_fclose(Infile);
  
  return(url_id);
}


char *fm_geturl(int url_id,char *infile_name)
{
  char *url_txt;
  void get_url_byid(int Url_id, char **url);
 

  Type_mem_ops = EXT_MEM;
  Cache_percentage = 0;
  Report_position=0;
  
  if (check_bwi_suffix(infile_name) == 0)
    fatal_error("The file name must end with .bwi -main-\n");  
  
  my_open_file(infile_name);
  
  init_bwi_cache();
  
  get_url_byid(url_id, &url_txt);
  
  my_fclose(Infile);
  
  return(url_txt);
}
