#
# If pkg-config isn't installed on your system, comment the following lines and
# set the fasterac flags as indicated with your own paths:
#
#FASTERAC_CFLAGS = -I/../Downloads/fasterac-2.17.1/fasterac-2.17/include
#FASTERAC_LIBS   = -L/../Downloads/fasterac-2.17.1/fasterac-2.17/libs
#
FASTERAC_CFLAGS = $(shell pkg-config  --cflags libfasterac)
FASTERAC_LIBS   = $(shell pkg-config  --libs   libfasterac)

ROOT_CFLAGS     = $(shell root-config --cflags)
ROOT_LIBS       = $(shell root-config --libs)

GSL_CFLAGS      = $(shell gsl-config --cflags)
GSL_LIBS        = $(shell gsl-config --libs)

CC        = g++
CFLAGS    = ${FASTERAC_CFLAGS} ${ROOT_CFLAGS} ${GSL_CFLAGS}
LIBS      = ${FASTERAC_LIBS}   ${ROOT_LIBS} ${GSL_LIBS}

SRCEXE    = $(shell ls *.C)
EXE       = $(SRCEXE:.C=)

all : $(EXE)

$(EXE): $(SRCEXE)
	${CC} $@.C -o $@ ${CFLAGS} ${LIBS}

clean :
	rm -f *.o
	rm -f $(EXE)


