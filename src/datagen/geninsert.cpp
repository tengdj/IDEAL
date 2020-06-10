#include "../geometry/MyPolygon.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <vector>
#include <string.h>
#include <fstream>
#include <iostream>
using namespace std;


int main(int argc, char** argv) {
	ifstream infile;
	infile.open(argv[1], ios::in | ios::binary);
	int index = 0;
	while(!infile.eof()){
		MyPolygon * poly = MyPolygon::read_polygon_binary_file(infile);
		if(!poly){
			continue;
		}
		poly->setid(++index);
		printf("%d|\"",poly->getid());
		poly->print_without_return(false);
		printf("\"\n");
		delete poly;
	}
	infile.close();
}
