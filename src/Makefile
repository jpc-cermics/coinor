SHELL = /bin/sh

SCIDIR=../../..
SCIDIR1=..\..\..

LIBRARY=libcoinor.a

OBJS= nspcoinmp.o nspclp.o nspcoinor-IN.o nspcoinread.o

include $(SCIDIR)/Makefile.incl

COIN_FLAGS=`/usr/bin/pkg-config gtk+-2.0 coinmp --cflags`
COIN_LIBS=`/usr/bin/pkg-config gtk+-2.0 coinmp --libs`

CFLAGS = $(CC_OPTIONS) $(COIN_FLAGS)
CXXFLAGS = $(CC_OPTIONS) $(COIN_FLAGS)

# extra libraries needed for linking 
# it is mandatory on win32 to give this extra argument.
OTHERLIBS=$(COIN_LIBS)

include $(SCIDIR)/config/Makeso.incl

all:: Makefile.mak

Makefile.mak	: Makefile
	$(SCIDIR)/scripts/Mak2VCMak Makefile

Makefile.amk	: Makefile
	$(SCIDIR)/scripts/Mak2ABSMak Makefile

distclean:: clean 

clean:: 
	@$(RM) *.o *.lo 
	@$(RM) -r */.libs 

