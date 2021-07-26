###########################################################################################
##
## target platform mahuika assumes:
## ml GCC
## ml Eigen
## ml SUNDIALS
## ml h5pp
##
## notes 
##    source code formating: clang-format -i ./src/*
##
###########################################################################################
## C++
CPP_DEFS = -std=c++17 -fopenmp -DNDEBUG
CPP = g++
CC = gcc

###########################################################################################
## NOTE: uncomment block for ONE of debug, profile or optimised

## debug
#CPP_FLAGS =  $(CPP_DEFS) -m64 -Og -Wall -Wextra

## profile
#CPP_FLAGS =  $(CPP_DEFS) -m64 -O3 -pg 

## optimised
CPP_FLAGS =  $(CPP_DEFS) -m64 -O3 -march=native -mtune=native -Wall -Wpedantic -g
#CPP_FLAGS =  $(CPP_DEFS) -m64 -O0 -Wall -Wpedantic -g
#CPP_FLAGS =  $(CPP_DEFS) -m64 -Og -Wextra -march=native -mtune=native -Wall -Wpedantic -g

## C
CFLAGS = -O3 -g

###########################################################################################
## NOTE: uncomment lines for ONE platform
CPP_INCLUDES = -I$(SRC)/inih/cpp
#CPP_INCLUDES = -I/usr/include  # mahuika 
#CPP_INCLUDES = -I/usr/local/include/eigen3 -I/usr/local/Cellar/boost/1.75.0_2/include -I/usr/local/Cellar/boost/1.75.0_2/lib # iMac

#CPP_LIBS = -lsundials_cvodes
#CPP_LIBS = -l/usr/local/Cellar/boost/1.75.0_2/lib  # iMac
#CPP_LIBS = -lboost_system
CPP_LIBS = -lstdc++ -lstdc++fs -lm -lsundials_cvode -lspdlog -lfmt -lhdf5_hl -lhdf5

###########################################################################################

EXE = ../run/mini-gland-sim
#OBJS = CVode.o
OBJS = main.o utils.o cMiniGland.o cDuct.o cLTree.o cSCell.o cCMesh.o ini.o INIReader.o
SRC = ./src

all: $(EXE)

ini.o: $(SRC)/inih/ini.c $(SRC)/inih/ini.h
	$(CC) $(CFLAGS) $(CPP_INCLUDES) -o $@ -c $<

INIReader.o: $(SRC)/inih/cpp/INIReader.cpp $(SRC)/inih/cpp/INIReader.h $(SRC)/inih/ini.h
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

main.o: $(SRC)/main.cpp $(SRC)/global_defs.hpp $(SRC)/cMiniGland.hpp $(SRC)/githash.h
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

utils.o: $(SRC)/utils.cpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cMiniGland.o: $(SRC)/cMiniGland.cpp $(SRC)/cMiniGland.hpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cDuct.o: $(SRC)/cDuct.cpp $(SRC)/cDuct.hpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cLTree.o: $(SRC)/cLTree.cpp $(SRC)/cLTree.hpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cSCell.o: $(SRC)/cSCell.cpp $(SRC)/cSCell.hpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

cCMesh.o: $(SRC)/cCMesh.cpp $(SRC)/cCMesh.hpp $(SRC)/utils.hpp $(SRC)/global_defs.hpp
	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

#cCVode.o: $(SRC)/cCVode.cpp $(SRC)/cCVode.hpp
#	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

#LSODA.o: $(SRC)/libsoda/LSODA.cpp $(SRC)/libsoda/LSODA.h $(SRC)/libsoda/helper.h
#	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

#cLSODA.o: $(SRC)/cLSODA.cpp $(SRC)/cLSODA.hpp $(SRC)/libsoda/LSODA.h
#	$(CPP) $(CPP_FLAGS) $(CPP_INCLUDES) -o $@ -c $<

.PHONY: $(SRC)/githash.h
$(SRC)/githash.h: ../.git/HEAD ../.git/index
	echo "#ifndef GIT_HASH" > $@ && \
	echo "#define GIT_HASH \"$(shell git describe --tags --always --dirty)\"" >> $@ && \
	echo "#endif" >> $@

$(EXE): $(OBJS)
	$(CPP) $(CPP_FLAGS) $(CPP_LIBS) -o $@ $+

clean:
	rm -f *.o githash.h
