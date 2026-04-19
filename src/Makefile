#
CC=gcc
COPTS= -g -O0 -Wall

hsh_table:	mesh.o hsh_table.o eigen.o libmesh6.o
	$(CC) $(COPTS) -o hsh_table.exe mesh.o  hsh_table.o eigen.o libmesh6.o -lm

insertion: mesh.o insertion.o eigen.o libmesh6.o
	$(CC) $(COPTS) -o insertion.exe mesh.o insertion.o eigen.o libmesh6.o -lm

compression : mesh.o compression.o eigen.o libmesh6.o
	$(CC) $(COPTS) -o compression.exe mesh.o compression.o eigen.o libmesh6.o -lm

libmesh6.o :	libmesh6.c libmesh6.h
	$(CC) -c $(COPTS)  -I. libmesh6.c

mesh.o :	mesh.c mesh.h libmesh6.h
	$(CC) -c $(COPTS)  -I. mesh.c

eigen.o :	eigen.c mesh.h
	$(CC) -c $(COPTS)  -I. eigen.c

hsh_table.o : hsh_table.c mesh.h
	$(CC) -c $(COPTS)  -I. hsh_table.c

insertion.o : insertion.c mesh.h
	$(CC) -c $(COPTS)  -I. insertion.c

compression.o : compression.c mesh.h
	$(CC) -c $(COPTS)  -I. compression.c
clean :
	-rm -f *.exe *.o *.mesh *.solb *.sol
