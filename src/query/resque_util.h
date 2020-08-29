/*
 * resque_util.h
 *
 *  Created on: Aug 28, 2020
 *      Author: teng
 */

#ifndef SRC_QUERY_RESQUE_UTIL_H_
#define SRC_QUERY_RESQUE_UTIL_H_

#include <geos/geom/PrecisionModel.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/operation/distance/DistanceOp.h>
#include <geos/geom/Point.h>
#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>
#include <geos/opBuffer.h>


using namespace geos;
using namespace geos::io;
using namespace geos::geom;
using namespace geos::operation::buffer;
using namespace geos::operation::distance;


using namespace std;

class resque_queue_element{
public:
	int id = 0;
	vector<Geometry *> geoms;
	resque_queue_element(int i){
		id = i;
	}
	~resque_queue_element(){
		for(Geometry *p:geoms){
			delete p;
		}
		geoms.clear();
	}
};



#endif /* SRC_QUERY_RESQUE_UTIL_H_ */
