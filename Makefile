COMPILER    = gnu

program = PI

source = \
main.c \
matrix.c \
init.c

obj = $(source:.c=.o)

# Standard Flags
CFLAGS := -std=gnu99 -Wall

# Linker Flags
LDFLAGS = -lm

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
  LDFLAGS += 
  CFLAGS += -Ofast -ffast-math
endif

# intel Compiler
ifeq ($(COMPILER),intel)
  CC = icc
  LDFLAGS += 
  CFLAGS += -O3 -xhost -ansi-alias -no-prec-div -DINTEL -vec-report6
endif

$(program): $(obj) PI_header.h
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj) data.dat
run:
	./PI
graph:
	gnuplot graph.gp
edit:
	vim -p $(source) PI_header.h
