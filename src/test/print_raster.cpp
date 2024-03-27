#include "Ideal.h"
#include "query_context.h"

int main() {
    query_context ctx;
    vector<Ideal *> source =
        load_binary_file("/home/qmh/data/has_child.idl", ctx);

    for (Ideal* p : source) {
        printf("BOX: lowx=%lf lowy=%lf, highx=%lf highy=%lf\n", p->mbr->low[0], p->mbr->low[1], p->mbr->high[0], p->mbr->high[1]);
        p->rasterization(100);
        p->MyPolygon::print();
        p->MyRaster::print();
    }

    return 0;
}
