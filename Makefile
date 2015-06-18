ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lMinuit
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = gcc -fPIC
LD            = gcc
LDFLAGS       = -g
SOFLAGS       = -shared -O3

CXXFLAGS       = $(ROOTCFLAGS) -fPIC
INCLUDE_FLAGS  = 
LDLIBS         = $(ROOTLIBS)
GLIBS          = $(ROOTGLIBS)

EXE           = cfit

INC 	      = cfit.h

SRC	      = cfit.cxx

OBJS          = cfit.o

LIB           = libCFIT.so

all: 	      $(LIB)

$(LIB):	      $(INC) $(SRC)
	      @echo "####### Generating dictionary"
	      @rootcint -f cfitDict.cxx -c -p $(CXXFLAGS) \
	      $(INCLUDE_FLAGS) -I. $(INC) LinkDef.h

	      @echo "####### Building library $(LIB)"
	      @$(CXX) $(SOFLAGS) $(CXXFLAGS) $(ROOTLIBS) $(INCLUDE_FLAGS) -I. $(SRC) \
	      cfitDict.cxx -o $(LIB) $(ROOTLIBS)
	      
	      @echo  "####### Removing generated dictionary"
	      @rm -f cfitDict.cxx cfitDict.h
	      @rm -f *.o

clean:
	      @rm -f $(OBJS) $(EXE) cfitDict.cxx cfitDict.h $(LIB)
