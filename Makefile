#
# Makefile for non-Microsoft compilers
#

## Linux  (uncomment the 2 lines below for compilation on Linux)
CXXFLAGS += -std=c++98 -Wall -ggdb
LDFLAGS += -L/sw/lib/gcc4.4/lib/x86_64/

OBJLIBS = blas.o lapack.o

## CygWin (uncomment the 2 lines below for compilation on CygWin)
#CXXFLAGS += -Wall
#LDFLAGS += 

MAIN = influmax

all: $(MAIN) generate_nets dominator_trees find_sets

opt: CXXFLAGS += -O4
opt: LDFLAGS += -O4 -L/sw/lib/gcc4.4/lib/x86_64/
opt: $(MAIN) generate_nets dominator_trees

linux: LDFLAGS = -lrt -L/usr/lib/gcc/x86_64-linux-gnu/4.4/
linux: $(MAIN) generate_nets dominator_trees

opt_linux: CXXFLAGS += -O4
opt_linux: LDFLAGS = -lrt -O4 /usr/lib/gcc/x86_64-linux-gnu/4.4/
opt_linux: $(MAIN) generate_nets dominator_trees

# COMPILE
$(MAIN): $(MAIN).cpp cascinf.cpp Snap.o lapack.o blas.o dnchbv.o dgexpv.o dgmatv.o dgpadm.o
	g++ $(LDFLAGS) -o $(MAIN) $(MAIN).cpp cascinf.cpp dgraph.cpp dgraph_snca.cpp dgraph_slt.cpp dgraph_sdom.cpp dgraph_lt.cpp dgraph_iter.cpp Snap.o blas.o lapack.o dnchbv.o dgexpv.o dgmatv.o dgpadm.o -lgfortran  -I./glib -I./snap

dominator_trees: dominator_trees.cpp cascinf.cpp Snap.o lapack.o blas.o dnchbv.o dgexpv.o dgmatv.o dgpadm.o
	g++ $(LDFLAGS) -o dominator_trees dominator_trees.cpp cascinf.cpp dgraph.cpp dgraph_snca.cpp dgraph_slt.cpp dgraph_sdom.cpp dgraph_lt.cpp dgraph_iter.cpp Snap.o blas.o lapack.o dnchbv.o dgexpv.o dgmatv.o dgpadm.o -lgfortran  -I./glib -I./snap

find_sets: find_sets.cpp cascinf.cpp Snap.o lapack.o blas.o dnchbv.o dgexpv.o dgmatv.o dgpadm.o
	g++ $(LDFLAGS) -o find_sets find_sets.cpp cascinf.cpp dgraph.cpp dgraph_snca.cpp dgraph_slt.cpp dgraph_sdom.cpp dgraph_lt.cpp dgraph_iter.cpp Snap.o blas.o lapack.o dnchbv.o dgexpv.o dgmatv.o dgpadm.o -lgfortran  -I./glib -I./snap

generate_nets: generate_nets.cpp cascinf.cpp Snap.o lapack.o blas.o dnchbv.o dgexpv.o dgmatv.o dgpadm.o
	g++ $(LDFLAGS) -o generate_nets generate_nets.cpp cascinf.cpp dgraph.cpp dgraph_snca.cpp dgraph_slt.cpp dgraph_sdom.cpp dgraph_lt.cpp dgraph_iter.cpp Snap.o blas.o lapack.o dnchbv.o dgexpv.o dgmatv.o dgpadm.o -lgfortran  -I./glib -I./snap

lapack.o:
	gfortran -static -c -m64 lapack.f
	
blas.o:
	gfortran -static -c -m64 blas.f
	
dnchbv.o:
	gfortran -static -c -m64 dnchbv.f

dgexpv.o:
	gfortran -static -c -m64 dgexpv.f

dgmatv.o:
	gfortran -static -c -m64 dgmatv.f

dgpadm.o:
	gfortran -static -c -m64 dgpadm.f
	
Snap.o: 
	g++ -c $(CXXFLAGS) ./snap/Snap.cpp -I./glib -I./snap
		
clean:
	rm -f   $(MAIN) generate_nets find_sets dominator_trees Snap.o  $(MAIN).exe generate_nets.exe find_sets.exe dominator_trees.exe
