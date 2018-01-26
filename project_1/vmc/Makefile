# Define the force removal command. This lets the clean-recipe remove the
# build-directory with all its subfolders
RM=rf -f

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

