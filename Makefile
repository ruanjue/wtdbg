#VERSION=1.0.0
#MINOR_VER=20160621

CC=gcc

ifeq (0, ${MAKELEVEL})
TIMESTAMP=$(shell date +"%s")
endif

ifeq (1, ${DEBUG})
CFLAGS=-g3 -W -Wall -O0 -DTIMESTAMP=$(TIMESTAMP) -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -mpopcnt -mssse3
else
CFLAGS=-g3 -W -Wall -O4 -DTIMESTAMP=$(TIMESTAMP) -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -mpopcnt -mssse3
endif

INSTALLDIR=/usr/local/bin
GLIBS=-lm -lpthread
GENERIC_SRC=mem_share.h string.h file_reader.h file_reader.c bitvec.h hashset.h sort.h list.h dna.h thread.h

PROGS=wtdbg wtdbg-cns

all: $(PROGS)

wtdbg: $(GENERIC_SRC) wtdbg.c dmo.h
	$(CC) $(CFLAGS) -o wtdbg file_reader.c wtdbg.c $(GLIBS)

wtdbg-cns: $(GENERIC_SRC) wtdbg-cns.c kswx.h ksw.h ksw.c dbgcns.h dagcns.h queue.h general_graph.h
	$(CC) $(CFLAGS) -o wtdbg-cns file_reader.c wtdbg-cns.c ksw.c $(GLIBS)

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out $(PROGS)

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out

install:
	cp -fvu $(PROGS) $(INSTALLDIR)
