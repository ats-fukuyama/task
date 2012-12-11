#GNU Makefile to build glf
#
# 03/15/10 C. Ludescher-Furth
-include ../share/Make.local
include ../share/Make.flags


libs: 
	@echo make libglf.a
	(cd src; ${MAKE} libs )

all: libs exec

exec:
	@echo make glftest
	(cd test; ${MAKE} exec )

clean:
	(cd src; $(MAKE) clean)
	(cd test; $(MAKE) clean)

realclean: clean
	(cd src; $(MAKE) realclean)
	(cd test; $(MAKE) realclean)


install:
	(cd src; $(MAKE) install)




