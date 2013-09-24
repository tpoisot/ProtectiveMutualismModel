prg = model
src = src.cpp

$(prg): $(src)
	g++ -Wall $(src) -o $(prg) -lgsl -lgslcblas -lm -O3 -DHAVE_INLINE -std=c++11

run: $(prg)
	./$(prg)
