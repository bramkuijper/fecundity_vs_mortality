xnumsolve : generate_cpp.py numsolve.cpp auxiliary.h
	./generate_cpp.py > numsolve2.cpp
	g++ -Wall -O3 -ggdb -o xnumsolve numsolve2.cpp -lm -lrt -lgsl -lgslcblas


clean :
	rm -rf numsolve
