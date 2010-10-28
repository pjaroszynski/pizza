CC=g++
CFLAGS=-W -Wall -pedantic -pedantic-errors -O9 -DNDEBUG

all: index

# pattern rule for all objects files
%.o: %.c *.h
	$(CC) -c $(CFLAGS) $< -o $@

%.o: %.cpp *.h
	$(CC) -c $(CFLAGS) $< -o $@

index: SSA.o bitarray.o HuffAlphabetRank.o Huffman_Codes.o bitrankw32int.o ds/ds.o ds/globals.o ds/helped.o ds/shallow.o ds/deep2.o ds/blind2.o 
	ar rc SSA.a $^
clean:
	rm -f *~ *.o ; cd ds ; make -f Makefile clean ; cd ..
cleanall:
	rm -f *~ *.o *.a; cd ds ; make -f Makefile cleanall ; cd ..
tarfile:
	tar zcvf SSA2_v1.1.tgz ?akefile *.c *.h *.cpp ds README COPYRIGHT
