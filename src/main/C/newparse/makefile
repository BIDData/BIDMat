
RELOPTS=-O2 -DNDEBUG -Wno-deprecated
DBGOPTS=-O2 -DDEBUG -g
RELLINK=
DBGLINK=
CC_OPTS=$(DBGOPTS)
LINK_OPTS=$(DBGLINK)
CC=gcc
CPP=g++
LD=g++
include gmake.options
.SUFFIXES: 
.SUFFIXES: .c .cpp .o .exe .lxc .flex

EXES=xmltweet.exe xmlwiki.exe
OBJS=newparse.o gzstream.o

.SECONDARY: xmltweet.lxc

all: $(EXES)

.flex.lxc: 
	flex -o $@ $< 

.o.exe: 
	$(LD) -o $@ $(OBJS) $< $(LINK_OPTS) -lz

.cpp.o: 
	$(CPP) $(CC_OPTS) -o $@ -c $<

.lxc.o: 
	$(CPP) $(CC_OPTS) -DYY_NO_UNISTD_H -o $@ -c -x c++ $<

$(EXES): $(OBJS)

newparse.o: newparse.h gzstream.h

gzstream.o: gzstream.h

debug: clean
	echo -e CC_OPTS=$(DBGOPTS)\\nLINK_OPTS=$(DBGLINK) > gmake.options

release: clean
	echo -e CC_OPTS=$(RELOPTS)\\nLINK_OPTS=$(RELLINK) > gmake.options

gmake.options:
	touch gmake.options

clean:
	rm -f $(EXES) *.o *.lxc
