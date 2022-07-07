/*
 * Tile.cpp
 *
 *  Created on: Jul 4, 2022
 *      Author: teng
 */

#include "partition.h"



Tile::Tile(){
	pthread_mutex_init(&lk, NULL);
	id = 0;
}

Tile::Tile(box b){
	pthread_mutex_init(&lk, NULL);
	id = 0;
	low[0] = b.low[0];
	low[1] = b.low[1];
	high[0] = b.high[0];
	high[1] = b.high[1];

}

Tile::~Tile(){
	targets.clear();
	objects.clear();
}

void Tile::lock(){
	pthread_mutex_lock(&lk);
}

void Tile::unlock(){
	pthread_mutex_unlock(&lk);
}

void Tile::merge(Tile *n){
	if(n->objects.size()==0&&n->targets.size()==0){
		return;
	}
	lock();
	update(*n);
	objects.insert(objects.end(),n->objects.begin(), n->objects.end());
	targets.insert(targets.end(),n->targets.begin(), n->targets.end());
	unlock();
}

bool Tile::insert(MyPolygon *obj, bool update_box){
	lock();
	if(update_box){
		update(*obj->getMBB());
	}
	objects.push_back(obj);
	unlock();
	return true;
}

bool Tile::insert_target(Point *tgt){
	lock();
	targets.push_back(tgt);
	unlock();
	return true;
}

void Tile::build_index(){
	lock();
	struct timeval start = get_cur_time();
	for(MyPolygon *p:objects){
		tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
	}
	indexing_latency = get_time_elapsed(start);
	unlock();
}

size_t Tile::conduct_query(){
	size_t found = 0;
	lock();
	struct timeval start = get_cur_time();
	for(Point *p:targets){
		vector<MyPolygon *> result = lookup(p);
		for(MyPolygon *poly:result){
			found += poly->contain(*p);
		}
		result.clear();
	}
	querying_latency = get_time_elapsed(start);
	unlock();
	return found;
}

bool Tile::lookup_tree(MyPolygon *obj, void *arg){
	vector<MyPolygon *> *results = (vector<MyPolygon *> *)arg;
	results->push_back(obj);
	return true;
}

vector<MyPolygon *> Tile::lookup(box *b){
	vector<MyPolygon *> results;
	tree.Search(b->low, b->high, lookup_tree, (void *)&results);
	return results;
}

vector<MyPolygon *> Tile::lookup(Point *p){
	vector<MyPolygon *> results;
	tree.Search((double *)p, (double *)p, lookup_tree, (void *)&results);
	return results;
}

void print_tiles(vector<Tile *> &boxes){
	MyMultiPolygon *cboxes = new MyMultiPolygon();
	for(Tile *p:boxes){
		MyPolygon *m = MyPolygon::gen_box(*p);
		cboxes->insert_polygon(m);
	}
	cboxes->print();
	delete cboxes;
}

double skewstdevratio(vector<Tile *> &tiles, int tag){
	if(tiles.size()==0){
		return 0;
	}
	size_t total = 0;
	for(Tile *t:tiles){
		size_t obj_num = 0;
		switch(tag){
		case 0:
			obj_num += t->objects.size();
			break;
		case 1:
			obj_num += t->targets.size();
			break;
		case 2:
			obj_num += t->objects.size();
			obj_num += t->targets.size();
			break;
		default:
			assert(false&&"can only be 0 1 2");
		}
		total += obj_num;
	}
	double avg = 1.0*total/tiles.size();
	double st = 0.0;
	for(Tile *t:tiles){
		size_t obj_num = 0;
		switch(tag){
		case 0:
			obj_num += t->objects.size();
			break;
		case 1:
			obj_num += t->targets.size();
			break;
		case 2:
			obj_num += t->objects.size();
			obj_num += t->targets.size();
			break;
		default:
			assert(false&&"can only be 0 1 2");
		}
		st += (obj_num-avg)*(obj_num-avg)/tiles.size();
	}
	return sqrt(st)/avg;
}
