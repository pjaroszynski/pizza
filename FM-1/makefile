CC = gcc

#these are for testing
#CFLAGS = -g -W -Wall -O2

#these are for maximum speed
CFLAGS = -g -O3 -fomit-frame-pointer -W -Wall -DDEBUG=0

# these are for profiling 
#CFLAGS = -g -pg -a -O2 -W -Wall -Winline -DDEBUG=0


.PHONY : all 
all : bwi bwsearch sa bwrand bwxml bwhuffw prova-lib


bwi: bwi.c fmindex.a
	$(CC) $(CFLAGS) -o bwi bwi.c fmindex.a -lm

bwsearch: bwsearch.c fmindex.a
	$(CC) $(CFLAGS) -o bwsearch bwsearch.c fmindex.a -lm

sa: common.h sa.o basic_routines.o qsufsort.o copy0xh2.o sacopy.h
	$(CC) $(CFLAGS) -o sa sa.o basic_routines.o qsufsort.o copy0xh2.o -lm


bwrand: bwrand.c fmindex.a
	$(CC) $(CFLAGS) -o bwrand bwrand.c fmindex.a -lm

bwxml:  bwxml.c fmindex.a
	$(CC) $(CFLAGS) -o bwxml bwxml.c fmindex.a -lm

bwhuffw:  bwhuffw.c fmindex.a
	$(CC) $(CFLAGS) -o bwhuffw bwhuffw.c fmindex.a -lm


prova-lib: prova-lib.c fmindex.a
	$(CC) $(CFLAGS) -o prova-lib prova-lib.c fmindex.a -lm


# archive containing the fm-index software
fmindex.a: basic_routines.o bwrand_lib.o compr_main.o decompr_routines.o \
	multihuf.o bwsearch_lib.o copy0xh2.o decompr_main.o \
	search_main.o common.h bzlib_private.h bwhuffw_lib.o huffman.o \
	qsufsort.o fmindex_api.o
	ar rc fmindex.a basic_routines.o bwrand_lib.o compr_main.o \
		decompr_routines.o multihuf.o fmindex_api.o \
		bwsearch_lib.o copy0xh2.o decompr_main.o search_main.o \
		bzlib_private.h bwhuffw_lib.o huffman.o qsufsort.o



# pattern rule for all objects files
%.o: %.c *.h
	 $(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f *.o *.a

tarfile:
	tar zcvf bwi.tgz makefile *.c *.h *template \
                         bwh bwhc bwhd bwhc1 scripts/*corpus scripts/*dict \
                         scripts/*.awk tools/*.c tools/dic2txt






