CC = g++

main: main.o utility.o ga.o
	$(CC) main.o utility.o ga.o -o main -g
main.o: main.cpp utility.h ga.h
	$(CC) main.cpp -c -g
utility.o: utility.cpp utility.h
	$(CC) utility.cpp -c -g
ga.o: ga.cpp ga.h
	$(CC) ga.cpp -c -g
