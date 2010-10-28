
#include "common.h"

int main(int argc, char *argv[])
{

  char *fm_geturl(int url_id,char *infile_name);
  int fm_search(char *pattern, int pat_len, char *infile_name);
  int fm_getid(char *url,char *infile_name);

  char *pattern = strdup(argv[1]);
  char *file =strdup(argv[2]);
  int urlid1 = atoi(argv[3]);
  char *url;
  int n,urlid2;

  if(argc<4)
    fatal_error("pochi parametri\n");

  printf("pattern= %s\n",pattern);

  n = fm_search(pattern,strlen(pattern),file);
  printf("occorrenze = %d\n",n);

  url = fm_geturl(urlid1,file);
  printf("url = %s for id=%d\n",url,urlid1);
  
  urlid2 = fm_getid(url,file);
  printf("url_id= %d for url= %s\n",urlid2,url);


  return 0;
}

