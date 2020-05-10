/* This contains headers (dependencies) of RESQUE spatial processing engine */
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>
#include <cstdlib>
#include <getopt.h>
#include <time.h>
#include <vector>
#include <string.h>

#include <geos/geom/PrecisionModel.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/operation/distance/DistanceOp.h>
#include <geos/geom/Point.h>
#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>
#include <geos/opBuffer.h>
#include <spatialindex/SpatialIndex.h>

#include "util.h"

// Constants used for building the R-tree
#define FillFactor 0.9
#define IndexCapacity 10 
#define LeafCapacity 50
#define COMPRESS true

const int OSM_SRID = 4326;

