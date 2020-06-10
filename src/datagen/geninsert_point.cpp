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
	int index = 0;
	string input_line;
	while(getline(cin,input_line)){
		Point *p = Point::read_one_point(input_line);
		if(!p){
			continue;
		}
		printf("%d|\"",++index);
		p->print_without_return();
		printf("\"\n");
		delete p;
	}
}
