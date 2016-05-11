CC          = gcc
CFLAGS      = 
LIB         = -lhdf5_hl -lhdf5 -lz -lm 

INCHDF5   = -I/usr/local/include/ -I/home/bustamsn/Software/hdf5-1.8.10/include
LIBHDF5   = -L/usr/local/lib/ 

all: NLLSF
     
NLLSF: NLLSF.c
	$(CC) $(CFLAGS) -o $@ NLLSF.c $(INCHDF5) $(LIBHDF5) $(LIB)

clean: 
	rm -f *.h5 *.o \
        h5_rdwt 