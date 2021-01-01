/*
 * Poly2Tri Copyright (c) 2009-2010, Poly2Tri Contributors
 * http://code.google.com/p/poly2tri/
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * * Neither the name of Poly2Tri nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include <cstdlib>
#include <time.h>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <signal.h>

using namespace std;

#include "../triangulate/poly2tri.h"
#include "../geometry/MyPolygon.h"
#include <map>
using namespace p2t;

void segfault_sigaction(int signal, siginfo_t *si, void *arg)
{
    exit(0);
}

int main(int argc, char* argv[]){

	struct sigaction sa;

	memset(&sa, 0, sizeof(struct sigaction));
	sigemptyset(&sa.sa_mask);
	sa.sa_sigaction = segfault_sigaction;
	sa.sa_flags   = SA_SIGINFO;

	sigaction(SIGSEGV, &sa, NULL);

	query_context global_ctx;
	global_ctx = get_parameters(argc, argv);
	global_ctx.query_type = QueryType::contain;

//	global_ctx.sort_polygons = true;
//	global_ctx.source_polygons = MyPolygon::load_binary_file(global_ctx.source_path.c_str(),global_ctx);
//
//	ofstream of;
//	of.open(global_ctx.target_path, ios::out | ios::binary);
//	size_t numpolygons = global_ctx.source_polygons.size();
//	of.write((char *)&numpolygons, sizeof(numpolygons));
//	for(int i=0;i<numpolygons;i++){
//		unsigned int si = global_ctx.source_polygons[i]->offset+sizeof(numpolygons)+sizeof(unsigned int)*numpolygons;
//		of.write((char *)&si, sizeof(unsigned int));
//	}
//	of.close();

//	return 0;

	MyPolygon *polygon = MyPolygon::load_binary_file_single(global_ctx.source_path.c_str(),global_ctx, global_ctx.vpr);
	vector<TrPoint*> polyline;
	for(int i=0;i<polygon->boundary->num_vertices-1;i++){
		polyline.push_back(new TrPoint(polygon->boundary->x[i],polygon->boundary->y[i]));
	}
	struct timeval start = get_cur_time();
	CDT* cdt = new CDT(polyline);
	cdt->Triangulate();
	vector<Triangle*> triangles = cdt->GetTriangles();
	logt("triangulation",start);
	cout<<"processed successful,"<<global_ctx.vpr<<endl;



	//  cout<<"MULTIPOLYGON(";
	//  for(int i=0;i<triangles.size();i++){
	//	  if(i>0){
	//		  printf(",");
	//	  }
	//	  Triangle *t = triangles[i];
	//	  printf("((%f %f,%f %f,%f %f,%f %f))",
	//			  t->point(0)->x,t->point(0)->y,
	//			  t->point(1)->x,t->point(1)->y,
	//			  t->point(2)->x,t->point(2)->y,
	//			  t->point(0)->x,t->point(0)->y
	//	  );
	//  }
	//  cout<<")"<<endl;





  // Cleanup
  
  
  return 0;
}
