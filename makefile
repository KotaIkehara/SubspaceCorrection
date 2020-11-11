#--------------------------------------------
#Hokudai Subsystem A
#--------------------------------------------
FC		= ifort
CC		= icc
CCFLAG	= -O3 -qopenmp -ipo -xCORE-AVX512
BLAS	= -mkl=parallel

LD		= $(CC)
LDFLAG  = $(CCFLAG) $(BLAS) -lm

MAIN	= main
LOAD	= a.out

OBJS	= main.o

.PHONY : clean
.PRECIOUS : $(LOAD)
.SUFFIXES : .o .c .f .f90


all : $(LOAD)

clean :
	rm -f *.o *.out *.mod

$(LOAD) : $(OBJS)
	rm -f $(LOAD)
	$(LD) -o $(LOAD) $(OBJS) $(LDFLAG) $(DEF)
	@echo "$(LOAD) is now up-to-date"


.c.o : 
	$(CC) $(CCFLAG) -c $< $(DEF)

.f.o : 
	$(FC) $(CCFLAG) -c $< $(DEF)

.f90.o : 
	$(FC) $(CCFLAG) -c $< $(DEF)

