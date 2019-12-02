cplex = /opt/ibm/ILOG/CPLEX_Studio128/cplex
concert = /opt/ibm/ILOG/CPLEX_Studio128/concert
compile = g++ -DIL_STD
cflags  = -c -O3
xcflags  = $(cflags) -isystem$(cplex)/include -isystem$(concert)/include
lflags  = -L$(cplex)/lib/x86-64_linux/static_pic -L$(concert)/lib/x86-64_linux/static_pic -lconcert -lilocplex -lcplex -lpthread -ldl

mbvt : main.o
	$(compile) -o mbvt main.o $(lflags)

main.o : main.cpp globais.hpp
	$(compile) main.cpp $(xcflags)

clean:
	rm -f mbvt *.o
