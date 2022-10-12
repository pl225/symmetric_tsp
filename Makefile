FLAG = -O3 -m64 -Wall

CHRISTOFIDES=christofides/

OBJ = $(CHRISTOFIDES)Graph.o $(CHRISTOFIDES)Matching.o $(CHRISTOFIDES)BinaryHeap.o $(CHRISTOFIDES)TSPLIB_parser.o

main.cpp:
	g++ $(FLAG) -c main.cpp -o main.o

build: main.o
	g++ $(FLAG) main.o $(OBJ) -o exe