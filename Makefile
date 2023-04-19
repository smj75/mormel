F77=f77

CFLAGS=-O -C -Wall -Wsurprising
LFLAGS=

# EXECUTABLES

all : mormel

mormel : mormel.o lib_veltem.o lib_mel.o lib_k03.o lib_mb88.o lib_res.o lib_comp.o lib_td.o Makefile
	${F77} ${LFLAGS} -o mormel mormel.o lib_veltem.o lib_mel.o lib_k03.o lib_mb88.o lib_res.o lib_comp.o lib_td.o

clean : Makefile
	/bin/rm -f core a.out mormel mormel.o \
          lib_veltem lib_veltem.o lib_mel lib_mel.o lib_td.o \
          lib_k03.o lib_mb88.o lib_res.o katz katz.o blobx blobx.o \
          lib_comp.o *.O
         


# OBJECTS

mormel.o : mormel.f
	${F77} ${CFLAGS} -c wj.f
lib_veltem.o : lib_veltem.f
	${F77} ${CFLAGS} -c lib_veltem.f
lib_mel.o : lib_mel.f
	${F77} ${CFLAGS} -c lib_mel.f
lib_k03.o : lib_k03.f
	${F77} ${CFLAGS} -c lib_k03.f
lib_mb88.o : lib_mb88.f
	${F77} ${CFLAGS} -c lib_mb88.f
lib_res.o : lib_res.f
	${F77} ${CFLAGS} -c lib_res.f
lib_comp.o : lib_comp.f
	${F77} ${CFLAGS} -c lib_comp.f
lib_td.o : lib_td.f
	${F77} ${CFLAGS} -c lib_td.f
#OBJS=${MAIN:%.f=%.o} ${LIBS:%.f=%.o} 
 
