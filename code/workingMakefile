CC = g++

main: main.o utility.o
	$(CC) main.o utility.o -o main -g
main.o: main.cpp utility.h
	$(CC) main.cpp -c -g
utility.o: utility.cpp utility.h
	$(CC) utility.cpp -c -g
