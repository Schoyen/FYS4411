# Define the force removal command. This lets the clean-recipe remove the
# build-directory with all its subfolders
RM=rm -r
MK=mkdir -p

# Define the target directories
SDIR=src
IDIR=include
BDIR=build
ODIR=$(BDIR)/obj

# Define the compiler to use
CXX=g++

# Define compiler flags
#
#  -I$(IDIR) => Look in include-directory for header files
#  -g => Debug flag, this lets us debug with Valgrind and gdb
#  -Wall => Include all compiler warnings
#  -std=c++11 => Compiler as C++11
#  -O3 => Optimize as much as possible
CXXFLAGS=-I$(IDIR) -g -Wall -std=c++11 -O3

# Define libraries to use
LIBS=-larmadillo

# Add headers as seen in $(IDIR)
_HEADERS=simple_vmc.h
# Tell make to look in $(IDIR) for the header-files defined in $(_HEADERS)
HEADERS=$(patsubst %, $(IDIR)/%, $(_HEADERS))

# Add object files to be generated
_OBJECTS=simple_vmc.o
# Tell make to look in $(ODIR) for the object files defined $(_OBJECTS)
OBJECTS=$(patsubst %, $(ODIR)/%, $(_OBJECTS))

$(ODIR)/%.o: $(SDIR)/%.cpp $(HEADERS) | build
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(LIBS)

simple_vmc.x: $(OBJECTS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

.PHONY: clean build

build:
	$(MK) $(BDIR) $(ODIR)

clean:
	$(RM) $(BDIR) simple_vmc.x
