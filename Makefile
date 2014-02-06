##
## Makefile for ProtoCellSelect
##  
## Made by Bart Markvoort
##

##############################
# Complete this to make it ! #
##############################
NAME 	= protocell
SRC	= main.cpp simbox.cpp exception.cpp fileio.cpp fileioutils.cpp random.cpp protocell.cpp 
INCL  	= simbox.h exception.h fileio.h fileioutils.h random.h protocell.h

################
# Optional add #
################
INCPATH   = #-I.           # path of include file
OBJOPT  = -O3 -Wall	 # option for obj
FSTOPT  = -ffast-math -ffinite-math-only
EXEOPT  = -O3 -Wall        # option for exe (-lefence ...)
LIBPATH   = #-L.           # path for librairies ... 
FLAGS	= -DDEBUG
LIBFLAGS = -shared

#####################
# Macro Definitions #
#####################
CC 	= g++
MAKE 	= make
SHELL	= /bin/sh
OBJS 	= $(SRC:.cpp=.o)
RM 	= /bin/rm -f 	
COMP	= gzip -9v
UNCOMP	= gzip -df
STRIP	= strip

CFLAGS  = $(OBJOPT) $(FSTOPT) $(FLAGS) $(IPATH)
LDFLAGS = $(EXEOPT) $(FSTOPT) $(FLAGS) $(LPATH)

.SUFFIXES: .h.Z .c.Z .h.gz .c.gz .c.z .h.z 

########################
# Compile Instructions #
########################

all:	$(NAME)

$(NAME): $(OBJS) $(SRC) $(MAINFILE) $(INCL)  
	$(CC) $(OBJS) $(LDFLAGS) -o $(NAME) 
#	$(STRIP) ./$(NAME) # if you debug ,don't strip ...

depend:
	gcc $(IPATH) -MM $(SRC) 

clean:
	-$(RM) $(OBJS) *~

fclean:
	-$(RM) $(NAME)

comp: clean
	$(COMP) $(INCL) $(SRC)

new:
	make clean
	make

ucomp: 
	$(UNCOMP) $(SRC) $(INCL)

.c.Z.c .h.Z.h .c.gz.c .h.gz.h .c.z.c .h.z.h :
	 -$(UNCOMP) $<

.cpp.o:
	$(CC) $(CFLAGS) -c $< 

