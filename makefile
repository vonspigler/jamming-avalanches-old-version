SRC = SRC/
INC = SRC/INC/
CMP = COMP/
OPT_BEG = -std=c++11 -Wall -Wextra
OPT_END = -O3 #-lgsl -lgslcblas -lm
CXX = mpic++

all: MAIN.bin

MAIN.bin: COMP/other.o COMP/vectors.o COMP/cells.o COMP/inc.o $(INC)other.h $(INC)vectors.h $(INC)cells.h $(INC)inc.h COMP/MAIN.o;\
  $(CXX) $(OPT_BEG) COMP/other.o COMP/vectors.o COMP/cells.o COMP/inc.o COMP/MAIN.o -o $@ $(OPT_END)

COMP/MAIN.o: $(SRC)MAIN.cpp $(INC)other.h $(INC)vectors.h $(INC)cells.h $(INC)inc.h;\
  $(CXX) $(OPT_BEG) $(SRC)MAIN.cpp -c -o $@ $(OPT_END)

COMP/other.o: $(INC)other.cpp $(INC)other.h $(INC)inc.h;\
  $(CXX) $(OPT_BEG) $(INC)other.cpp -c -o $@ $(OPT_END)

COMP/cells.o: $(INC)cells.cpp $(INC)cells.h $(INC)vectors.h $(INC)inc.h;\
  $(CXX) $(OPT_BEG) $(INC)cells.cpp -c -o $@ $(OPT_END)

COMP/vectors.o: $(INC)vectors.cpp $(INC)vectors.h;\
  $(CXX) $(OPT_BEG) $(INC)vectors.cpp -c -o $@ $(OPT_END)

COMP/inc.o: $(INC)inc.cpp $(INC)inc.h;\
  $(CXX) $(OPT_BEG) $(INC)inc.cpp -c -o $@ $(OPT_END)

.PHONY: clean
clean:;\
  rm -f *~ MAIN.bin *.o fit.log;\
  rm -f COMP/*
