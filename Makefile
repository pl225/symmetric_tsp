FLAG = -O3 -m64 -Wall

CHRISTOFIDES=christofides/

OBJ = $(CHRISTOFIDES)Graph.o $(CHRISTOFIDES)Matching.o $(CHRISTOFIDES)BinaryHeap.o $(CHRISTOFIDES)TSPLIB_parser.o $(CHRISTOFIDES)Christofides.o $(CHRISTOFIDES)MST.o

main.cpp:
	g++ $(FLAG) -c main.cpp -o main.o

lagrangeano_principal.cpp:
	g++ $(FLAG) -c lagrangeano_principal.cpp -o lagrangeano_principal.o

build: main.o lagrangeano_principal.o
	g++ $(FLAG) main.o lagrangeano_principal.o $(OBJ) -o exe