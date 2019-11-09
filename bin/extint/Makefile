CC       = cc
CFLAGS   = -ansi -O6 -fomit-frame-pointer -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic -fno-schedule-insns -fschedule-insns2 -fstrict-aliasing -m64
CLINKERS = -lm -static

#############
PROGRAMS = extint
EXT_OBJS = extint.o structures.o
 
#############
all:	$(PROGRAMS)

extint:	$(EXT_OBJS) structures.h
	$(CC) $(CFLAGS) -o $@ $(EXT_OBJS) $(CLINKERS)

#############
clean-obj:
	rm -f *.o

clean: clean-obj
	rm -f $(PROGRAMS)
