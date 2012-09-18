COPT =  -g -Wall -Wno-sign-compare -O  -fPIC -fexceptions #-DNDEBUG
CC    = g++


# Default target
all: maps2bin make_simple_coverage make_gcbins doc_walker make_exp gem2masks make_masked_windows split_gem make_fasta_random lower_resolution gcContent
	@echo "compilation done"


BGZF.o : bamtools/BGZF.cpp
	$(CC) -c $(COPT) bamtools/BGZF.cpp -o BGZF.o
BamReader.o : bamtools/BGZF.o bamtools/BamReader.cpp
	$(CC) -c $(COPT) bamtools/BamReader.cpp -o BamReader.o
maps2bin: BGZF.o BamReader.o maps2bin.o
	$(CC) maps2bin.o BGZF.o BamReader.o -lz -o maps2bin
maps2bin.o : maps2bin.cpp
	$(CC) -c $(COPT) maps2bin.cpp -o maps2bin.o
make_simple_coverage: make_simple_coverage.o
	$(CC) make_simple_coverage.o -o make_simple_coverage
make_simple_coverage.o : make_simple_coverage.cpp
	$(CC) -c $(COPT) make_simple_coverage.cpp -o make_simple_coverage.o
make_gcbins: make_gcbins.o
	$(CC) make_gcbins.o -o make_gcbins
make_gcbins.o : make_gcbins.cpp
	$(CC) -c $(COPT) make_gcbins.cpp -o make_gcbins.o
doc_walker: doc_walker.o
	$(CC) doc_walker.o -o doc_walker
doc_walker.o : doc_walker.cpp
	$(CC) -c $(COPT) doc_walker.cpp -o doc_walker.o
make_exp: make_exp.o
	$(CC) make_exp.o -o make_exp
make_exp.o : make_exp.cpp
	$(CC) -c $(COPT) make_exp.cpp -o make_exp.o
gem2masks: gem2masks.o
	$(CC) gem2masks.o -o gem2masks
gem2masks.o : gem2masks.cpp
	$(CC) -c $(COPT) gem2masks.cpp -o gem2masks.o
make_masked_windows: make_masked_windows.o
	$(CC) make_masked_windows.o -o make_masked_windows
make_masked_windows.o : make_masked_windows.cpp
	$(CC) -c $(COPT) make_masked_windows.cpp -o make_masked_windows.o
split_gem: BGZF.o BamReader.o split_gem.o
	$(CC) split_gem.o BGZF.o BamReader.o -lz -o split_gem
split_gem.o : split_gem.cpp
	$(CC) -c $(COPT) split_gem.cpp -o split_gem.o
make_fasta_random: make_fasta_random.o
	$(CC) make_fasta_random.o -o make_fasta_random
make_fasta_random.o : make_fasta_random.cpp
	$(CC) -c $(COPT) make_fasta_random.cpp -o make_fasta_random.o
lower_resolution: lower_resolution.o
	$(CC) lower_resolution.o -o lower_resolution
lower_resolution.o : lower_resolution.cpp
	$(CC) -c $(COPT) lower_resolution.cpp -o lower_resolution.o
gcContent: gcContent.o
	$(CC) gcContent.o -o gcContent
gcContent.o : gcContent.cpp
	$(CC) -c $(COPT) gcContent.cpp -o gcContent.o

# Target deleting unwanted files
clean:
	rm -f *.o *~ 
