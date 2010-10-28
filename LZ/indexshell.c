
#include "lzindex.h"


static void txtload(char *fname, byte **text, ulong *n)
 {
    int f;
    struct stat sdata;
 
    /* read the file into memory */
    if (stat(fname,&sdata) != 0) {
       fprintf (stderr,"Cannot stat file %s\n",fname);
          exit(1);
    }
 
    *n = sdata.st_size;
    *text = malloc(*n+1);
 
    f = open (fname,O_RDONLY);
    if (f == -1) {
       fprintf (stderr,"Cannot open file %s\n",fname);
       exit(1);
    }
 
    if (read(f,*text,*n) != *n) {
       fprintf (stderr,"Cannot read file %s\n",fname);
       exit(1);                                                              
    }
 
    /* text terminator */
    (*text)[*n] = 0;
     
    close(f);
 }                                                                            

main (int argc, char **argv) 
 { 
    byte *text;
    lzindex *I = malloc(sizeof(lzindex));
    uint i,level;
    ulong size,n;
     

    if (argc != 2) {
       fprintf (stderr,"Usage: %s <file> indexes <file> and produces files "
                       "<file>.*,\n       you can then delete <file>\n",argv[0]);
       exit (1);
    }
    // load text
    printf("Loading %s...\n",argv[1]);
    txtload(argv[1],&text,&n);

    // build index
    printf("Building index...\n");
    build_index(text, n, NULL, (void **)&I);
    // free text
    printf("Freeing text...\n");
    free(text);
    printf("Size of the Text: %d bytes\n", n);
    index_size(I,&size);
    printf("Size of LZ-index: %u bytes\n", size);
    printf("The size of LZ-index is %f times the text size\n", (float)size/n);
    // save index into disk
    printf("Saving index...\n");
    save_index(I, argv[1]);
    // free structures
    printf("Freeing index\n");
    free_index(I);

    return 0;
 }
