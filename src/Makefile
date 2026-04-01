#
CC=gcc
COPTS= -g -O0 -Wall

mesh:	mesh.o main_mesh.o eigen.o libmesh6.o
	$(CC) $(COPTS) -o mesh.exe mesh.o  main_mesh.o eigen.o libmesh6.o -lm

comp_test : mesh.o comp_test.o eigen.o libmesh6.o
	$(CC) $(COPTS) -o comp_test.exe mesh.o comp_test.o eigen.o libmesh6.o -lm

libmesh6.o :	libmesh6.c libmesh6.h
	$(CC) -c $(COPTS)  -I. libmesh6.c

mesh.o :	mesh.c mesh.h libmesh6.h
	$(CC) -c $(COPTS)  -I. mesh.c

eigen.o :	eigen.c mesh.h
	$(CC) -c $(COPTS)  -I. eigen.c

main_mesh.o :	main_mesh.c mesh.h
	$(CC) -c $(COPTS)  -I. main_mesh.c

comp_test.o : comp_test.c mesh.h
	$(CC) -c $(COPTS)  -I. comp_test.c
clean :
	-rm mesh.exe comp_test.exe comp_test.o mesh.o main_mesh.o eigen.o libmesh6.o 
