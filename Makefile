CXX = g++
LFLAGS = -static -lm -lboost_program_options


.PHONY: all

all: 
	$(CXX) src/*.cpp -o JoRGS $(LFLAGS)

.PHONY: clean

clean:
	rm -f JoRGS