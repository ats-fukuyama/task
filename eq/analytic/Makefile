LIBS=-L /usr/local/lib -lgsp -lgdp -L /usr/X11R6/lib -lX11

fc=pgf77

SRC=analytic.f

OBJ=$(SRC:.f=.o)

FLAGS=-O2 -fast

analytic : $(OBJ)
	$(fc) -o $@ $(OBJ) $(FLAGS) $(LIBS) 

.f.o :
	$(fc) -c $(FLAGS) $< -o $@

check : 
	ftnchek -array=0 -usage=no-com-block-unused $(SRC)| less

clean :
	-rm *~ *.o
