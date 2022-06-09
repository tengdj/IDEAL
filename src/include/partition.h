/*
 * partition.hpp
 *
 *  Created on: Jun 8, 2022
 *      Author: teng
 */

#ifndef SRC_PARTITION_PARTITION_HPP_
#define SRC_PARTITION_PARTITION_HPP_

#include "MyPolygon.h"

typedef enum {
	STR = 0,
	SLC,
	BOS,
	HC,
	FG,
	QT,
	BSP
}PARTITION_TYPE;

vector<box *> genschema_str(vector<box *> &geometries, size_t part_num);
vector<box *> genschema_slc(vector<box *> &geometries, size_t part_num);
vector<box *> genschema_bos(vector<box *> &geometries, size_t part_num);
vector<box *> genschema_hc(vector<box *> &geometries, size_t part_num);

vector<box *> genschema_fg(vector<box *> &geometries, size_t part_num);
vector<box *> genschema_qt(vector<box *> &geometries, size_t part_num);
vector<box *> genschema_bsp(vector<box *> &geometries, size_t part_num);

vector<box *> genschema(vector<box *> &geometries, size_t part_num, PARTITION_TYPE type);

PARTITION_TYPE parse_partition_type(const char *type);
#endif /* SRC_PARTITION_PARTITION_HPP_ */
