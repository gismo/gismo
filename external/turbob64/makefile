# powturbo  (c) Copyright 2016-2019
# Linux: "export CC=clang" "export CXX=clang". windows mingw: "set CC=gcc" "set CXX=g++" or uncomment the CC,CXX lines
CC ?= gcc
CXX ?= g++

#CC=powerpc64le-linux-gnu-gcc
#DEBUG=-DDEBUG -g
#DEFS=$(DEBUG)
#CFLAGS=$(DEBUG)
#------- OS/ARCH -------------------
ifneq (,$(filter Windows%,$(OS)))
  OS := Windows
  CC=gcc
  CXX=g++
  ARCH=x86_64
else
  OS := $(shell uname -s)
  ARCH := $(shell uname -m)

ifneq (,$(findstring aarch64,$(CC)))
  ARCH = aarch64
else ifneq (,$(findstring powerpc64le,$(CC)))
  ARCH = ppc64le
endif
endif

ifeq ($(ARCH),ppc64le)
  CFLAGS=-mcpu=power9 -mtune=power9
  MSSE=-D__SSSE3__
else ifeq ($(ARCH),aarch64)
  CFLAGS+=-march=armv8-a 
ifneq (,$(findstring clang, $(CC)))
  CFLAGS+=-march=armv8-a -falign-loops -fomit-frame-pointer
else
  CFLAGS+=-march=armv8-a 
endif
  MSSE=-march=armv8-a
else ifeq ($(ARCH),$(filter $(ARCH),x86_64 ppc64le))
  CFLAGS=-march=native
  MSSE=-mssse3
endif

ifeq (,$(findstring clang, $(CC)))
DEFS+=-falign-loops
endif
#$(info ARCH="$(ARCH)")

ifeq ($(OS),$(filter $(OS),Linux GNU/kFreeBSD GNU OpenBSD FreeBSD DragonFly NetBSD MSYS_NT Haiku))
LDFLAGS+=-lrt
endif
ifeq ($(STATIC),1)
LDFLAGS+=-static
endif

#FPIC=-fPIC

all: tb64app 
#libtb64.so

ifeq ($(NCHECK),1)
DEFS+=-DNB64CHECK
else
ifeq ($(FULLCHECK),1)
DEFS+=-DB64CHECK
endif
endif

turbob64c.o: turbob64c.c
	$(CC) -O3 $(MARCH) $(DEFS) $(FPIC) -fstrict-aliasing  $< -c -o $@ 

tb64app.o: tb64app.c
	$(CC) -O3 $(DEFS) $< -c -o $@ 

turbob64d.o: turbob64d.c
	$(CC) -O3 $(MARCH) $(DEFS) $(FPIC) -fstrict-aliasing $< -c -o $@ 

turbob64sse.o: turbob64sse.c
	$(CC) -O3 $(MSSE) $(DEFS) $(FPIC) -fstrict-aliasing $< -c -o $@ 

turbob64avx.o: turbob64sse.c
	$(CC) -O3 $(DEFS) $(FPIC) -march=corei7-avx -mtune=corei7-avx -mno-aes -fstrict-aliasing $< -c -o turbob64avx.o 

turbob64avx2.o: turbob64avx2.c
	$(CC) -O3 $(FPIC) $(DEFS) -march=haswell -fstrict-aliasing -falign-loops $< -c -o $@ 

turbob64avx512.o: turbob64avx512.c
	$(CC) -O3 $(FPIC) -march=skylake-avx512 -mavx512vl -fstrict-aliasing -falign-loops $< -c -o $@ 

_tb64.o: _tb64.c
	$(CC) -O3 $(FPIC) -I/usr/include/python2.7 $< -c -o $@ 

LIB=turbob64c.o turbob64d.o turbob64sse.o
ifeq ($(ARCH),x86_64)
LIB+=turbob64avx.o turbob64avx2.o
endif

ifeq ($(BASE64),1)
include xtb64make
endif

ifeq ($(AVX512),1)
DEFS+=-DUSE_AVX512
LIB+=turbob64avx512.o
endif

#_tb64.so: _tb64.o
#	gcc -shared $^ -o $@

libtb64.so: $(LIB)
	gcc -shared $^ -o $@
	cp libtb64.so ~/.local/lib/
	./python/tb64/build.py
	cp _tb64.so ~/.local/lib/

tb64app: $(LIB) tb64app.o 
	$(CC) -O3 $(LIB) tb64app.o $(LDFLAGS) -o tb64app

tb64bench: $(LIB) tb64bench.o 
	$(CC) -O3 $(LIB) tb64bench.o $(LDFLAGS) -o tb64bench

tb64test: $(LIB) tb64test.o 
	$(CC) -O3 $(LIB) tb64test.o $(LDFLAGS) -o tb64test
	
	
.c.o:
	$(CC) -O3 $(CFLAGS) $(MARCH) $< -c -o $@

clean:
	@find . -type f -name "*\.o" -delete -or -name "*\~" -delete -or -name "core" -delete -or -name "tb64app" -delete -or -name "_tb64.so" -delete -or -name "_tb64.c" -delete -or -name "xlibtb64.so"  -delete -or -name "libtb64.a"

cleanw:
	del /S *.o

