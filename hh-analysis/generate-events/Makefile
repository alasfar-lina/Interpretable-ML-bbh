FC = gfortran
FFLAGS  = -ffixed-line-length-none -fno-align-commons  -ffast-math -fno-automatic
ROOT = /users/pep/alasfarl/HEP_tools/RootTuple-1.0.0/lib

LIBROOT= $(ROOT)/libRootTuple.so
#COLLIER
COLPATH= /users/pep/alasfarl/HEP_tools/COLLIER-1.2.5
export LD_RUN_PATH:=$(COLPATH)
export LD_RUN_PATH:=$(ROOT)
INCC = $(COLPATH)/modules
LIBCOLL = $(COLPATH)/libcollier.so
AR = ar rcs

LHAPDF_CONFIG= lhapdf-config
CXXFLAGS+= $(shell $(LHAPDF_CONFIG) --cxxflags)
LIBSLHAPDF= -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF -lstdc++

LIBS=$(LIBCOLL)
LIBS+=$(LIBSLHAPDF)
LIBS+=$(LIBROOT)
OBJSXSGG = \
	pythia-6.4.28.o hh.o  alphas.o  xsgg.o prog.o


OBJSXSQQ = \
	pythia-6.4.28.o hh.o  alphas.o prog.o xsqq.o

.f.o:
	$(FC) -c $(FFLAGS) -I$(INCC)  $*.f


hhgg: $(OBJSXSGG)
		$(FC) -L$(ROOT) $(FFLAGS) -I$(INCC) $(LIBS) $(OBJSXSGG)  -o run $(CXXFLAGS)

hhqq: $(OBJSXSQQ)
		$(FC) -L$(ROOT) $(FFLAGS) -I$(INCC) $(LIBS) $(OBJSXSQQ)  -o run $(CXXFLAGS)		
clean:
	rm -f *.o
