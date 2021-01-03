/*
 * GEOSTool.h
 *
 *  Created on: Sep 4, 2020
 *      Author: teng
 */

#ifndef SRC_GEOMETRY_GEOSTOOL_H_
#define SRC_GEOMETRY_GEOSTOOL_H_

#include <geos/geom/PrecisionModel.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/operation/distance/DistanceOp.h>
#include <geos/geom/Point.h>
#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>
#include <geos/opBuffer.h>
#include "../geometry/MyPolygon.h"


using namespace geos;
using namespace geos::io;
using namespace geos::geom;
using namespace geos::operation::buffer;
using namespace geos::operation::distance;


using namespace std;
void process_geometries(query_context *ctx, vector<unique_ptr<Geometry>> &dest, bool procee_target = false);
void process_points(query_context *ctx, vector<unique_ptr<Geometry>> &dest);


#endif /* SRC_GEOMETRY_GEOSTOOL_H_ */
