CC = clang
CFLAGS = -Wall -std=c99 -pedantic

all: A3

clean:
	rm -f *.o *.so A3

_molecule.so: molecule_wrap.o
	$(CC) molecule_wrap.o -shared -lmol -lpython3.7m -dynamiclib -L/home/undergrad/1/irfan/Cis2750/A2 -L/usr/lib/python3.7/config-3.7m-x86_64-linux-gnu -o _molecule.so

molecule_wrap.o: molecule_wrap.c
	$(CC) $(CFLAGS) -c molecule_wrap.c -fPIC -I/usr/include/python3.7m -o molecule_wrap.o

libmol.so: mol.o
	$(CC) mol.o -shared -o libmol.so

mol.o: mol.c mol.h
	$(CC) $(CFLAGS) -c mol.c -fPIC -o mol.o

molecule.py:
	swig -python molecule.i

A3: libmol.so _molecule.so molecule.py
	$(CC) $(CFLAGS) -L. -lmol -lm -o A3