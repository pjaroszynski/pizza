CXX = g++
CXXFLAGS = -O2

OBJS_random_dna = random_dna.o

all: random_dna

%.o : %.cc *.hh
	$(CXX) $(CXXFLAGS) -c $< -o $@

random_dna: $(OBJS_random_dna)
	$(CXX) $(CXXFLAGS) $^ -o $@

data: random_dna
	./random_dna 50 > data/dna.50M
	./random_dna 50 > data/dna.100M
	./random_dna 50 > data/dna.200M

.PHONY: data
