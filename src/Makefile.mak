#-----------------------------
# generated from Makefile: DO NOT EDIT
# -----------------------------
SHELL = /bin/sh

SCIDIR=../../..
SCIDIR1=..\..\..

LIBRARY=libcoinor.lib

OBJS= nspcoinmp.obj nspclp.obj nspcoinor-IN.obj nspcoinread.obj

include $(SCIDIR)/Makefile.incl.mak

COIN_FLAGS=`/usr/bin/pkg-config gtk+-2.0 coinmp --cflags`
COIN_LIBS=`/usr/bin/pkg-config gtk+-2.0 coinmp --libs`

CFLAGS = $(CC_OPTIONS) $(COIN_FLAGS)
CXXFLAGS = $(CC_OPTIONS) $(COIN_FLAGS)

# extra libraries needed for linking 
# it is mandatory on win32 to give this extra argument.
OTHERLIBS=$(COIN_LIBS)

include $(SCIDIR)/config/Makeso.incl



Makefile.mak	: Makefile
	$(SCIDIR)/scripts/Mak2VCMak Makefile

Makefile.libmk	: Makefile
	$(SCIDIR)/scripts/Mak2ABSMak Makefile

distclean:: clean 

clean:: 
	@$(RM) *.obj *.lo 
	@$(RM) -r */.libs 

