BUILD = ../build
VPATH = .:$(BUILD)

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h



LD = gfortran
CF90 = gfortran
CF77 = gfortran

FFLAGS = -c -g -fdefault-real-8 -C -u -Xlist -O2 -std=gnu -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -Wunused -Werror -O3

F77FLG =
F90FLG =
LINKFLAGS =
CPPARCH =
CPPMODEL  =

CPPFLAGS = -traditional -P $(CPPARCH) $(CPPMODEL)

LIBS =    ./libenkfanalysis.a ./libsampling.a -lblas -lfftw3 -llapack

BINDIR = $(HOME)/bin

CPP = /usr/bin/cpp

# Rules for running cpp and updating files in $(BUILD) directory
.H.h:
	@mkdir -p $(BUILD)
	@rm -f ./$(BUILD)/$*.h
	cat MODEL.CPP $*.H | $(CPP) $(CPPFLAGS) > ./$(BUILD)/$*.h


.F90.o:
	@mkdir -p $(BUILD)
	@rm -f ./$(BUILD)/$*.f90
	@cat MODEL.CPP $*.F90 | $(CPP) $(CPPFLAGS) > ./$(BUILD)/$*.f90
	cd $(BUILD) ; $(CF90) $(FFLAGS) $(F90FLG) -o $*.o $*.f90

.F.o:
	@mkdir -p $(BUILD)
	@rm -f ./$(BUILD)/$*.f
	@cat MODEL.CPP $*.F | $(CPP) $(CPPFLAGS) > ./$(BUILD)/$*.f
	cd $(BUILD) ; $(CF77) $(FFLAGS) $(F77FLG) -o $*.o $*.f



include source.files
include target.mk

INCFILES =$(INCSOURCE:.H=.h)
F90FILES =$(F90SOURCE:.F90=.f90)
F77FILES =$(F77SOURCE:.F=.f)
OBJECTS = $(F90SOURCE:.F90=.o) $(F77SOURCE:.F=.o)


all: new $(TARGET) install

libs: libsampling.a libenkfanalysis.a

libsampling.a:
	cd ../../EnKF_sampling/lib; make -f makefile.ge BUILD=../../EnKF_analysis/build

libenkfanalysis.a:
	cd ../../EnKF_analysis/lib; make -f makefile.ge BUILD=../../EnKF_analysis/build

$(TARGET): $(INCFILES) $(OBJECTS)
	@mkdir -p $(BUILD)
	cd ./$(BUILD) ; $(LD) $(LINKFLAGS) -o $(TARGET) $(OBJECTS) $(LIBDIR) $(LIBS)

install:
	cp $(BUILD)/$(TARGET) $(BINDIR)

clean:
	cd $(BUILD) ; rm -rf *.o *.mod *.f90 rii_files lib*.a

new: source depend

source:
	../bin/mksource.sh > source.files

depend:
	../bin/mkdepend_linux.pl | sort -u > depends.file

tags: $(F90FILES) $(F77FILES) $(INCFILES)
	ctags *.F90
	#f90tags.sh

include depends.file
