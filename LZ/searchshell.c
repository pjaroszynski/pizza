

#include "lzindex.h"

main (int argc, char **argv)
 { 
    byte patt[1024], *snippet_text;
    lzindex *I;
    uint i,level;
    ulong numocc, *occ, *snippet_lengths, pos, j,
          length, numc;

    if (argc != 2) { 
       fprintf (stderr,"Usage: %s <file> searches index <file>.*\n",
                          argv[0]);
       exit (1);
    }

    // load text
    printf ("Loading index %s...\n",argv[1]);
    load_index(argv[1],(void **) &I);
    // answer queries
    printf ("Ready to answer queries\n");
    while (true) {
       printf ("pattern = ");
       for (i=0;((patt[i]=getchar())!='\n')&& (i < 1023);i++);
       if (i==0) continue;
       patt[i] = 0;
       if (!strcmp(patt,"END")) break;
       if ((patt[0]=='c') && (patt[1]==':')) {
          strcpy(patt,patt+2); 
          length = strlen(patt); 
          if (length==0) continue;
          count(I, patt, length, &numocc);
       }
       else if ((patt[0]=='p') && (patt[1]==':')) {
               strcpy(patt,patt+2); 
               length = strlen(patt);
               if (length==0) continue;
               locate(I, patt, length, &occ, &numocc);
               for (i = 0; i < numocc; i++)
                  printf("%i\n", occ[i]);
               free(occ); 
            }
            else {
               if ((strlen(patt) >= 2) && (patt[0]=='s') && (patt[1]==':')) strcpy(patt,patt+2);                      
               length = strlen(patt);
               if (length==0) continue;
               numc = 20;
               display(I, patt, length, numc, &numocc, 
                       &snippet_text, &snippet_lengths);
               for (i = 0, pos = 0; i < numocc; i++, pos = i*(length+2*numc)) {
                  printf("(***Ocurrence %i***)\n", i+1);
                  for (j = 0; j < snippet_lengths[i]; j++, pos++) 
                     printf("%c",snippet_text[pos]);
                  printf("\n");
               }
               free(snippet_lengths);
               free(snippet_text);
            }
       printf ("%s: %i occurrences\n",patt,numocc);
    }
    // free structures
    printf ("Freeing index\n");
    free_index(I);
    return 0;
 }
