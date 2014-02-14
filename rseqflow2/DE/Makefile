all: DE_Image

DE_Image: DE_Image.o 
	g++ -o DE_Image DE_Image.o -Wall -Wno-deprecated -g -lgd

DE_Image.o: DE_Image.cpp
	g++ -c DE_Image.cpp -Wall -Wno-deprecated -g -lgd

clean: 
	rm -f DE_Image

