CC=g++
CFLAGS=`root-config --cflags --libs` -Wall
INC=-ISrc/include

#Debug with OPT=-g
#Reduce compilation time and be able to debug with OPT=-O0
ifdef OPT
OPTM=$(OPT)
else
OPTM=-O3
endif

SRCDIR := Src
SRC_BUILDDIR := $(SRCDIR)/build
PHYSDIR := Physics
PHYS_BUILDDIR := $(PHYSDIR)/build

HEADER := $(shell find $(SRCDIR)/include -type f -name '*.h')
SOURCES := $(shell find $(SRCDIR) -type f -name '*.cpp')
PHYSICS := $(shell find $(PHYSDIR) -type f -maxdepth 1 -name '*.cpp')

OBJS := $(patsubst $(SRCDIR)/%,$(SRC_BUILDDIR)/%,$(SOURCES:.cpp=.o))
PHYS := $(patsubst $(PHYSDIR)/%,$(PHYS_BUILDDIR)/%,$(PHYSICS:.cpp=.o)) 
MAIN := $(shell find Physics -type f -maxdepth 1 -name '*.cpp')
EXE := $(patsubst Physics/%,%,$(MAIN:.cpp=))

all:$(EXE)
.PHONY : all

$(EXE): $(PHYS) $(OBJS) 
	$(CC) -o $@ Physics/$@.cpp $(OBJS) $(OPTM) $(INC) $(CFLAGS)


$(PHYS_BUILDDIR)/%.o: $(PHYSDIR)/%.cpp
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(SRC_BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<


clean:
	$(RM) -r $(EXE) $(OBJS) $(PHYS) *.dSYM
