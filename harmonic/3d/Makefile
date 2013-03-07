FC = gfortran
FFLAGS = -Wall -Wextra -march=native -O3
LDFLAGS =
LIBS = -llapack
COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)
OBJS =
OBJS += qmc.o
all: qmc
qmc:	$(OBJS)
		$(LINK) -o $@ $^ $(LIBS)
%.o: %.f90
		$(COMPILE) -o $@ -c $<
clean:
		$(RM) sums $(OBJS) *.mod
