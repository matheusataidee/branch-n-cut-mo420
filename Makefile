cplex = /opt/ibm/ILOG/CPLEX_Studio129/cplex
concert = /opt/ibm/ILOG/CPLEX_Studio129/concert
compile = g++ -DIL_STD
cflags  = -c -O3
xcflags  = $(cflags) -isystem$(cplex)/include -isystem$(concert)/include -Iinc/
lflags  = -L$(cplex)/lib/x86-64_linux/static_pic -L$(concert)/lib/x86-64_linux/static_pic -lconcert -lilocplex -lcplex -lpthread -ldl

mbvt : main.o utils.o
	$(compile) -o bnc main.o utils.o $(lflags)

main.o : src/main.cpp inc/globais.hpp
	$(compile) src/main.cpp $(xcflags)

utils.o : src/utils.cpp inc/utils.hpp
	$(compile) src/utils.cpp $(xcflags)

clean:
	rm -f bnc *.o
