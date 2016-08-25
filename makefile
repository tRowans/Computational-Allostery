CC = gcc
# Compiler flags: all warnings + debugger meta-data
CFLAGS = -g -std=c99 -I/usr/include

# External libraries: only math in this example
LIBS = -L/usr/local/lib -L/usr/lib -lm 

# Pre-defined macros for conditional compilation
DEFS = -DDEBUG_FLAG -DEXPERIMENTAL=0

# The final executable program file, i.e. name of our program
BIN = main

# Object files from which $BIN depends
OBJS = lib.o sysutils.o rwutils.o sest.o mcutils.o avo.o drutils.o

# This default rule compiles the executable program
$(BIN): $(OBJS) $(BIN).c
	$(CC) $(CFLAGS) $(DEFS) $(OBJS) $(BIN).c -o avo $(LIBS)

# This rule compiles each module into its object file
%.o: %.c lib.h
	$(CC) -c $(CFLAGS) $(DEFS) $< -o $@ $(LIBS)

clean:
	rm -f *~ *.o $(BIN)



