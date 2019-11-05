INC=-I./Tools/c++
CXXFLAGS=-std=c++11 -O3 -g
LIBFLAGS=-pthread

all: clean  fastXML_predict


fastXML_predict:
	g++ -o fastXML_predict $(CXXFLAGS) $(INC) fastXML_predict.cpp fastXML.cpp ./Tools/c++/mat.cpp $(LIBFLAGS)

clean:
	rm -f  fastXML_predict

