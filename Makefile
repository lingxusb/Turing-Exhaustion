###########################################
#Makefile for Turing-topo programs
###########################################

INC = -I "/usr/include" -I "/usr/local/include"
LIB = -lgsl -lgslcblas -lm

CC = gcc

PRG=turing
OBJ=main.o dynfuc.o find_robust_topology.o generate_net.o get_sample.o make_text.o solveroot.o tools.o turing_criteria.o
SOR=dynfuc.c find_robust_topology.c generate_net.c get_sample.c main.c make_text.c solveroot.c tools.c turing_criteria.c

$(PRG):$(OBJ)
		$(CC) -c $(SOR) $(LIB)
		$(CC) -o $(PRG) $(OBJ) $(LIB) 
			
%.o: %.c
		$(CC) $(INC) -c $<


clean:
	rm -f $(OBJ) $(PRG)
