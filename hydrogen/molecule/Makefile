FC = gfortran
FFLAGS = -Wall -Wextra -march=native -O3 -ffast-math
LDFLAGS =
LIBS = 

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

OBJS =
OBJS += h2.o
OBJS += qmc.o

all: qmc

qmc:	$(OBJS)
		$(LINK) -o $@ $^ $(LIBS)

%.o: %.f90
		$(COMPILE) -o $@ -c $<

clean:
		$(RM) sums $(OBJS) *.mod
