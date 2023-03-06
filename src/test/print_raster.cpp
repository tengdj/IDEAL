#include "MyPolygon.h"
#include "query_context.h"


int main(){
	query_context ctx;
	vector<MyPolygon *> source = load_binary_file("/gisdata/ideal/idl/has_child.idl",ctx);

	for(MyPolygon *p:source){
		if(p->boundary->num_vertices>10000){
			p->rasterization(100);
			p->print();
			p->get_rastor()->print();
			return 0;
		}
	}

	return 0;
}
