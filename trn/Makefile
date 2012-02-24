all:
	(cd ../lib; make)
	(cd ../plx; make)
	(cd ../pl; make)
	(cd ../eq; make)
	(cd trcore; make trcomm.o trlib.o)
	(cd trmodels; make)
	(cd trgraphics; make)
	(cd trcore; make)

clean:
	(cd trcore; make clean)
	(cd trmodels; make clean)
	(cd trgraphics; make clean)

veryclean:
	(cd trcore; make veryclean)
	(cd trmodels; make veryclean)
	(cd trgraphics; make veryclean)
