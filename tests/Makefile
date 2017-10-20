SHELL = /bin/sh

.SUFFIXES: .dia .tst .vl

TOP=..

include $(TOP)/Path.incl 
#include $(TOP)/$(SCIDIR)/Makefile.incl

NSP=$(TOP)/$(SCIDIR)/bin/nsp

TESTS= linprog.dia clpsparse.dia examples.dia alfonsi.dia

VL=$(TESTS:.dia=.vl) 

tests: $(TESTS) 

vl  : $(VL)

.tst.dia :
	@echo -n "testing $*: "
	@$(NSP) -ns -nw -e "exec ../src/loader.sce;nsptest('$*.tst');quit" 

#export G_DEBUG=gc-friendly;export G_SLICE=always-malloc;valgrind --suppressions=tests/valgrind.supp --leak-check=full  bin/scilex -nw 

.tst.vl :
	@echo -n "valgrind $*: "
	@setenv G_DEBUG gc-friendly; setenv G_SLICE always-malloc;valgrind --suppressions=valgrind.supp --leak-check=full $(TOP)/$(SCIDIR)/bin/scilex -ns -nw -e "exec ../src/loader.sce;exec('$*.tst');quit" >& $*.vl
	@echo  " $*.vl"

distclean:
	@echo cleaning tests
	@$(RM) $(TESTS) 

clean:	distclean

