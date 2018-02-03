###########################################
#Makefile for Turing-topo programs
###########################################

INC = -I "/usr/include" -I "/usr/local/include"
LIB = -L /usr/local/lib -lgsl -lgslcblas -lm

CC = gcc

PRG=turing
OBJ=main.o dynfuc.o find_robust_topology.o generate_net.o get_sample.o make_text.o solveroot.o tools.o turing_criteria.o

$(PRG):$(OBJ)
		$(CC) $(LIB) -o $(PRG) $(OBJ)
			
%.o: %.c
		$(CC) $(INC) -c $<


clean:
	rm -f $(OBJ) $(PRG)
