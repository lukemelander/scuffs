# Build SCUFFS (Structured Compressible Upwind Flow Finite-volume Solver)
#   make            optimised build -> bin/scuffs
#   make debug      bounds checking + backtraces
#   make clean      remove build/ and bin/
# Override the compiler or flags on the command line, e.g. make FC=ifx FFLAGS=-O3

FC         = gfortran
FFLAGS     = -O2
DEBUGFLAGS = -O0 -g -fbacktrace -fcheck=all -Wall

SRCDIR   = src
BUILDDIR = build
BINDIR   = bin
TARGET   = $(BINDIR)/scuffs

SRCS = $(SRCDIR)/util.f90 $(SRCDIR)/solver.f90 $(SRCDIR)/main.f90
OBJS = $(patsubst $(SRCDIR)/%.f90,$(BUILDDIR)/%.o,$(SRCS))

all: $(TARGET)

$(TARGET): $(OBJS) | $(BINDIR)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

# Module files (.mod) are written to and read from build/
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -J$(BUILDDIR) -c $< -o $@

# Module dependencies: util -> solver -> main
$(BUILDDIR)/solver.o: $(BUILDDIR)/util.o
$(BUILDDIR)/main.o: $(BUILDDIR)/util.o $(BUILDDIR)/solver.o

$(BUILDDIR) $(BINDIR):
	mkdir -p $@

debug: FFLAGS = $(DEBUGFLAGS)
debug: clean all

clean:
	rm -rf $(BUILDDIR) $(BINDIR)

.PHONY: all debug clean
