CXX=g++
AR=ar cr
CPPFLAGS=-std=c++11 -Ofast
LDFLAGS=-L lib/ -lgzstream -lpthread -lz
INC=-I include/ 

GZPATH=$(CURDIR)/src/gzstream/
LDIR=$(CURDIR)/lib/
BDIR=$(CURDIR)/bin/
SRC=$(CURDIR)/src/


all: gzstream
	$(CXX) -o $(BDIR)/add_dosage $(SRC)/add_dosage.cpp $(CPPFLAGS) $(INC) $(LDFLAGS)

gzstream:
	cd $(GZPATH) && make && cp libgzstream.a $(LDIR)/ 

clean:
	rm $(BDIR)/* && rm $(LDIR)/*
	cd $(GZPATH) && make clean
