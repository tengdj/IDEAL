#include "Ideal.h"
#include "query_context.h"

int main(int argc, char** argv) {
    query_context global_ctx;
	global_ctx = get_parameters(argc, argv);
	global_ctx.query_type = QueryType::contain;

	global_ctx.source_polygons = load_polygons_from_path(global_ctx.source_path.c_str(),global_ctx);
		
	preprocess(&global_ctx);


    return 0;
}
