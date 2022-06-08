/*
 * partition.hpp
 *
 *  Created on: Jun 8, 2022
 *      Author: teng
 */

#ifndef SRC_PARTITION_PARTITION_HPP_
#define SRC_PARTITION_PARTITION_HPP_

#include "MyPolygon.h"

vector<Pixel *> genschema_str(vector<Pixel *> &geometries, size_t part_num);
vector<Pixel *> genschema_slc(vector<Pixel *> &geometries, size_t part_num);
vector<Pixel *> genschema_bos(vector<Pixel *> &geometries, size_t part_num);
vector<Pixel *> genschema_hc(vector<Pixel *> &geometries, size_t part_num);

vector<Pixel *> genschema_fg(vector<Pixel *> &geometries, size_t part_num);
vector<Pixel *> genschema_qt(vector<Pixel *> &geometries, size_t part_num);
vector<Pixel *> genschema_bsp(vector<Pixel *> &geometries, size_t part_num);



#endif /* SRC_PARTITION_PARTITION_HPP_ */
