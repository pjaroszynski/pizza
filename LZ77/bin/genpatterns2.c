
// Extracts random patterns from a file

#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
static int Seed;
#define ACMa 16807
#define ACMm 2147483647
#define ACMq 127773
#define ACMr 2836
#define hi (Seed / ACMq)
#define lo (Seed % ACMq)

static int fst = 1;

	/*
	 * returns a random integer in 0..top-1 
	 */

int
aleat (unsigned int top)
{
	long test;
	struct timeval t;
	if (fst)
	{
		gettimeofday (&t, NULL);
		Seed = t.tv_sec * t.tv_usec;
		fst = 0;
	}
	{
		Seed = ((test =
			 ACMa * lo - ACMr * hi) > 0) ? test : test + ACMm;
		return ((double) Seed) * top / ACMm;
	}
}

unsigned char* readFile(char const* filename, unsigned int* length){
    unsigned char* str = NULL;    
    FILE* file = NULL;
    int read;
    int len;
    *length = 0; 
    file = fopen(filename,"r");
    /*Check for validity of the file.*/
    if(file == NULL){
        fprintf(stderr,"ERROR: Can't open file: %s\n",filename);
        return NULL;
    }
    /*Calculate the file length in bytes. This will be the length of the source string.*/
    fseek(file, 0, SEEK_END);
    len = (unsigned int)ftell(file);
    fseek(file, 0, SEEK_SET);
    //str = (unsigned char*)malloc((len+1)*sizeof(char));
    str = (unsigned char*)malloc(len+1);
    if(str == NULL){
        fprintf(stderr,"ERROR: readFile: Out of memory when trying to alloc %d bytes\n",len+1);
        return NULL;
    }
    read = fread(str, sizeof(char), len, file);
    fclose(file);
    if(read<len){
        fprintf(stderr,"ERROR: readFile: Read less chars than size of file\n");
        return NULL;
    }
    *length = len;
    return str;
}

unsigned char getChar(unsigned char* sigma_chars, unsigned int sigma){
    return sigma_chars[aleat(sigma-1)];
}
unsigned char* getSigma(char* filename, unsigned int* sigma){
    unsigned int count[256];
    for(unsigned int i=0;i<256;i++)count[i]=0;
    unsigned int len;
    unsigned char* text = readFile(filename, &len);
    for(unsigned int i=0;i<len;i++){
        count[text[i]]+=1;
    }
    *sigma=0;
    for(unsigned int i=0;i<256;i++)if(count[i]!=0)*sigma=*sigma+1;
    unsigned char* sigma_chars = (unsigned char*)malloc(*sigma*sizeof(unsigned char));
    *sigma=0;
    for(unsigned int i=1;i<256;i++)if(count[i]!=0){
        sigma_chars[*sigma] = i;
        *sigma=*sigma+1;
    }
    return sigma_chars;
}

main (int argc, char **argv)
{
	int n, m, J, t;
	struct stat sdata;
	FILE *ifile, *ofile;
	unsigned char *buff;
	unsigned char *forbid, *forbide = NULL;

	if (argc < 5)
	{
		fprintf (stderr,
			 "Usage: genpatterns <file> <length> <number> <patterns file> <forbidden>\n"
			 "  randomly extracts <number> substrings of length <length> from <file>,\n"
			 "  avoiding substrings containing characters in <forbidden>.\n"
			 "  The output file, <patterns file> has a first line of the form:\n"
			 "    # number=<number> length=<length> file=<file> forbidden=<forbidden>\n"
			 "  and then the <number> patterns come successively without any separator.\n"
			 "  <forbidden> uses \\n, \\t, etc. for nonprintable chracters or \\cC\n"
			 "  where C is the ASCII code of the character written using 3 digits.\n\n");
		exit (1);
	}

	if (stat (argv[1], &sdata) != 0)
	{
		fprintf (stderr, "Error: cannot stat file %s\n", argv[1]);
		fprintf (stderr, " errno = %i\n", errno);
		exit (1);
	}
	n = sdata.st_size;

    unsigned int sigma;
    unsigned char* sigma_chars = getSigma(argv[1],&sigma);
	m = atoi (argv[2]);
	if ((m <= 0) || (m > n))
	{
		fprintf (stderr,
			 "Error: length must be >= 1 and <= file length"
			 " (%i)\n", n);
		exit (1);
	}

	J = atoi (argv[3]);
	if (J < 1)
	{
		fprintf (stderr, "Error: number of patterns must be >= 1\n");
		exit (1);
	}

	if (argc > 5) {
	} else
		forbid = NULL;

	ifile = fopen (argv[1], "r");
	if (ifile == NULL)
	{
		fprintf (stderr, "Error: cannot open file %s for reading\n", argv[1]);
		fprintf (stderr, " errno = %i\n", errno);
		exit (1);
	}

	buff = (unsigned char *) malloc (n);
	if (buff == NULL)
	{
		fprintf (stderr, "Error: cannot allocate %i bytes\n", n);
		fprintf (stderr, " errno = %i\n", errno);
		exit (1);
	}

	if (fread (buff, n, 1, ifile) != 1)
	{
		fprintf (stderr, "Error: cannot read file %s\n", argv[1]);
		fprintf (stderr, " errno = %i\n", errno);
		exit (1);
	}
	fclose (ifile);

	ofile = fopen (argv[4], "w");
	if (ofile == NULL)
	{
		fprintf (stderr, "Error: cannot open file %s for writing\n",
			 argv[4]);
		fprintf (stderr, " errno = %i\n", errno);
		exit (1);
	}

	if (fprintf (ofile, "# number=%i length=%i file=%s forbidden=%s\n",
		     J, m, argv[1],
		     forbid == NULL ? "" : (char *) forbid) <= 0)
	{
		fprintf (stderr, "Error: cannot write file %s\n", argv[4]);
		fprintf (stderr, " errno = %i\n", errno);
		exit (1);
	}

	for (t = 0; t < J; t++)
	{
		int j, l;
        for (l = 0; l < m; l++){
            buff[l] = getChar(sigma_chars,sigma);
        }
		for (l = 0; l < m; l++)
			if (putc (buff[j + l], ofile) != buff[j + l])
			{
				fprintf (stderr,
					 "Error: cannot write file %s\n",
					 argv[4]);
				fprintf (stderr, " errno = %i\n", errno);
				exit (1);
			}
	}

	if (fclose (ofile) != 0)
	{
		fprintf (stderr, "Error: cannot write file %s\n", argv[4]);
		fprintf (stderr, " errno = %i\n", errno);
		exit (1);
	}

	fprintf (stderr, "File %s successfully generated\n", argv[4]);
	free(forbide);
	exit (1);
}
