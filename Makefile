# ----- PROGRAM -----

PROG = combine.exe


# ----- COMPILER -----

# Linux
FC = gfortran
#FC = /usr/bin/f77

# Cygwin
#FC = /usr/bin/f77


# ----- FLAGS -----

# Linux
FFLAGS = -O2 

# Debug Linux
#FFLAGS = -O0 -march=pentium4 -fbounds-check -g -Wall

# Cygwin
#FFLAGS = -O2 -march=nocona -mfpmath=sse,387 -static

# Debug Cygwin
#FFLAGS = -O0 -march=nocona -fbounds-check -g -Wall


# ----- ACTIONS -----

RM = /bin/rm -f


# ----- OBJECTS -----

OBJ = combine.o


# ----- RULES -----

all:	$(PROG)

$(PROG):	$(OBJ)
		$(FC) $(FFLAGS) $(OBJ) -o $(PROG)
		
clean:
	$(RM) *.o $(PROG)
