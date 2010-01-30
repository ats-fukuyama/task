# makefile to build NCLASS_PT

SHELL     = /bin/sh
NAME      = NCLASS_PT
INCLUDES  =
SOURCES   = nclass_pt_dr.f \
            nclass_mod.f \
            rarray_copy.f \
            rarray_sum.f \
            rarray_zero.f \
            u_erf.f \
            u_lu_backsub.f \
            u_lu_decomp.f \
            write_mod.f
OBJECTS   = $(SOURCES:.f=.o)
LIBRARIES = 
DEBUG     = -g -C
OPTIMIZE  = 
LISTING   = 
SPECIAL   = -qfixed
#Compile line for HYDRA at DIII-D
COMPILE   = f90 -c $(OPTIMIZE) $(LISTING) $(SPECIAL)
#Load line for HYDRA at DIII-D
LOAD      = f90 -o $(NAME)                $(SPECIAL)
ARCHIVE   = ar rv
CATALOG   = ar tv
.SUFFIXES : .for
.for.f:
	 mv $(*F).for $(*F).f
.f.o:
	$(COMPILE) $< -o $@
################################################################################
$(NAME):    $(OBJECTS) $(LIBRARIES)
	$(LOAD)    $(OBJECTS) $(LIBRARIES)
################################################################################
nclass_pt_dr.o:    pamx_mi.inc \
                   pamx_ms.inc \
                   pamx_mz.inc \
                   nclass_mod.f \
                   rarray_copy.f \
                   rarray_sum.f \
                   rarray_zero.f \
                   write_mod.f
################################################################################
nclass_mod.o:      pamx_mi.inc \
                   pamx_ms.inc \
                   pamx_mz.inc \
                   rarray_zero.f \
                   u_erf.f \
                   u_lu_backsub.f \
                   u_lu_decomp.f
