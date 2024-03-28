#include "Ideal.h"
#include "query_context.h"

int main() {
    query_context ctx;
    vector<Ideal *> source =
        load_binary_file("/home/qmh/data/has_child.idl", ctx);

    for (Ideal* p : source) {
        p->rasterization(100);
        p->MyPolygon::print();
        p->MyRaster::print();
    }

    return 0;
}
