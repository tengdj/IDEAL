/*
 * rasterization.cpp
 *
 *  Created on: Jun 2, 2020
 *      Author: teng
 */

#include "MyPolygon.h"
#include <math.h>
#include <stack>

void MyPolygon::init_partition(const int dimx, const int dimy){
	assert(mbb);
	const double start_x = mbb->low[0];
	const double start_y = mbb->low[1];
	for(double i=0;i<=dimx;i++){
		vector<Pixel> v;
		for(double j=0;j<=dimy;j++){
			Pixel m;
			m.id[0] = i;
			m.id[1] = j;
			m.low[0] = i*step_x+start_x;
			m.high[0] = (i+1.0)*step_x+start_x;
			m.low[1] = j*step_y+start_y;
			m.high[1] = (j+1.0)*step_y+start_y;
			v.push_back(m);
		}
		partitions.push_back(v);
	};
}

void MyPolygon::evaluate_edges(const int dimx, const int dimy){
	// normalize
	assert(mbb);
	const double start_x = mbb->low[0];
	const double start_y = mbb->low[1];

	for(int i=0;i<get_num_vertices()-1;i++){
		double x1 = getx(i);
		double y1 = gety(i);
		double x2 = getx(i+1);
		double y2 = gety(i+1);

		int cur_startx = (x1-start_x)/step_x;
		int cur_endx = (x2-start_x)/step_x;
		int cur_starty = (y1-start_y)/step_y;
		int cur_endy = (y2-start_y)/step_y;

		if(cur_startx==dimx){
			cur_startx--;
		}
		if(cur_endx==dimx){
			cur_endx--;
		}

		int minx = min(cur_startx,cur_endx);
		int maxx = max(cur_startx,cur_endx);


		if(cur_starty==dimy){
			cur_starty--;
		}
		if(cur_endy==dimy){
			cur_endy--;
		}
		// should not happen for normal cases
		if(cur_startx>=dimx||cur_endx>=dimx||cur_starty>=dimy||cur_endy>=dimy){
			cout<<"xrange\t"<<cur_startx<<" "<<cur_endx<<endl;
			cout<<"yrange\t"<<cur_starty<<" "<<cur_endy<<endl;
			printf("xrange_val\t%f %f\n",(x1-start_x)/step_x, (x2-start_x)/step_x);
			printf("yrange_val\t%f %f\n",(y1-start_y)/step_y, (y2-start_y)/step_y);
		}
		assert(cur_startx<dimx);
		assert(cur_endx<dimx);
		assert(cur_starty<dimy);
		assert(cur_endy<dimy);


		partitions[cur_startx][cur_starty].status = BORDER;
		partitions[cur_endx][cur_endy].status = BORDER;

		//in the same pixel
		if(cur_startx==cur_endx&&cur_starty==cur_endy){
			continue;
		}

		if(y1==y2){
			//left to right
			if(cur_startx<cur_endx){
				for(int x=cur_startx;x<cur_endx;x++){
					partitions[x][cur_starty].leave(y1,RIGHT,i);
					partitions[x+1][cur_starty].enter(y1,LEFT,i);
				}
			}else { // right to left
				for(int x=cur_startx;x>cur_endx;x--){
					partitions[x][cur_starty].leave(y1, LEFT,i);
					partitions[x-1][cur_starty].enter(y1, RIGHT,i);
				}
			}
		}else if(x1==x2){
			//bottom up
			if(cur_starty<cur_endy){
				for(int y=cur_starty;y<cur_endy;y++){
					partitions[cur_startx][y].leave(x1, TOP,i);
					partitions[cur_startx][y+1].enter(x1, BOTTOM,i);
				}
			}else { //border[bottom] down
				for(int y=cur_starty;y>cur_endy;y--){
					partitions[cur_startx][y].leave(x1, BOTTOM,i);
					partitions[cur_startx][y-1].enter(x1, TOP,i);
				}
			}
		}else{
			// solve the line function
			double a = (y1-y2)/(x1-x2);
			double b = (x1*y2-x2*y1)/(x1-x2);

			int x = cur_startx;
			int y = cur_starty;
			while(x!=cur_endx||y!=cur_endy){
				bool passed = false;
				double yval = 0;
				double xval = 0;
				int cur_x = 0;
				int cur_y = 0;
				//check horizontally
				if(x!=cur_endx){
					if(cur_startx<cur_endx){
						xval = ((double)x+1)*step_x+start_x;
					}else{
						xval = (double)x*step_x+start_x;
					}
					yval = xval*a+b;
					cur_y = (yval-start_y)/step_y;
					//printf("y %f %d\n",(yval-start_y)/step_y,cur_y);
					if(cur_y>max(cur_endy, cur_starty)){
						cur_y=max(cur_endy, cur_starty);
					}
					if(cur_y<min(cur_endy, cur_starty)){
						cur_y=min(cur_endy, cur_starty);
					}
					if(cur_y==y){
						passed = true;
						// left to right
						if(cur_startx<cur_endx){
							partitions[x++][y].leave(yval,RIGHT,i);
							partitions[x][y].enter(yval,LEFT,i);
						}else{//right to left
							partitions[x--][y].leave(yval,LEFT,i);
							partitions[x][y].enter(yval,RIGHT,i);
						}
					}
				}
				//check vertically
				if(y!=cur_endy){
					if(cur_starty<cur_endy){
						yval = (y+1)*step_y+start_y;
					}else{
						yval = y*step_y+start_y;
					}
					xval = (yval-b)/a;
					int cur_x = (xval-start_x)/step_x;
					//printf("x %f %d\n",(xval-start_x)/step_x,cur_x);
					if(cur_x>max(cur_endx, cur_startx)){
						cur_x=max(cur_endx, cur_startx);
					}
					if(cur_x<min(cur_endx, cur_startx)){
						cur_x=min(cur_endx, cur_startx);
					}
					if(cur_x==x){
						passed = true;
						if(cur_starty<cur_endy){// bottom up
							partitions[x][y++].leave(xval, TOP,i);
							partitions[x][y].enter(xval, BOTTOM,i);
						}else{// top down
							partitions[x][y--].leave(xval, BOTTOM,i);
							partitions[x][y].enter(xval, TOP,i);
						}
					}
				}
				// for debugging, should never happen
				if(!passed){
					this->print(true);
					cout<<"dim\t"<<dimx<<" "<<dimy<<endl;
					printf("val\t%f %f\n",(xval-start_x)/step_x, (yval-start_y)/step_y);
					cout<<"curxy\t"<<x<<" "<<y<<endl;
					cout<<"calxy\t"<<cur_x<<" "<<cur_y<<endl;
					cout<<"xrange\t"<<cur_startx<<" "<<cur_endx<<endl;
					cout<<"yrange\t"<<cur_starty<<" "<<cur_endy<<endl;
					printf("xrange_val\t%f %f\n",(x1-start_x)/step_x, (x2-start_x)/step_x);
					printf("yrange_val\t%f %f\n",(y1-start_y)/step_y, (y2-start_y)/step_y);
				}
				assert(passed);
			}
		}
	}

	for(vector<Pixel> &rows:partitions){
		for(Pixel &p:rows){
			if(p.crosses.size()>0){
				p.status = BORDER;
			}
		}
	}
}

Pixel *MyPolygon::get_closest_pixel(Point p){
	assert(this->is_grid_partitioned());
	int pixx = this->get_pixel_x(p.x);
	int pixy = this->get_pixel_y(p.y);
	if(pixx<0){
		pixx = 0;
	}
	if(pixx>=partitions.size()){
		pixx = partitions.size()-1;
	}
	if(pixy<0){
		pixy = 0;
	}
	if(pixy>=partitions[0].size()){
		pixy = partitions[0].size()-1;
	}
	return &partitions[pixx][pixy];

}


vector<vector<Pixel>> MyPolygon::partition(int vpr){
	assert(vpr>0);

	getMBB();
	const double start_x = mbb->low[0];
	const double start_y = mbb->low[1];
	const double end_x = mbb->high[0];
	const double end_y = mbb->high[1];

	double multi = abs((end_y-start_y)/(end_x-start_x));
	int dimx = std::pow((get_num_vertices()/vpr)/multi,0.5);
	int dimy = dimx*multi;

	if(dimx==0){
		dimx = 1;
	}
	if(dimy==0){
		dimy = 1;
	}

	return partition(dimx, dimy);;
}

vector<vector<Pixel>> MyPolygon::partition(int dimx, int dimy){
	pthread_mutex_lock(&partition_lock);
	if(is_grid_partitioned()){
		pthread_mutex_unlock(&partition_lock);
		return partitions;
	}
	getMBB();
	const double start_x = mbb->low[0];
	const double start_y = mbb->low[1];
	const double end_x = mbb->high[0];
	const double end_y = mbb->high[1];

	step_x = (end_x-start_x)/dimx;
	step_y = (end_y-start_y)/dimy;

	if(step_x<0.00001){
		step_x = 0.00001;
		dimx = (end_x-start_x)/step_x+1;
	}
	if(step_y<0.00001){
		step_y = 0.00001;
		dimy = (end_y-start_y)/step_y+1;
	}

	init_partition(dimx, dimy);

	// edge crossing
	evaluate_edges(dimx, dimy);


	//scanline rendering
	for(int y=1;y<dimy;y++){
		bool isin = false;
		for(int x=0;x<dimx;x++){
			if(partitions[x][y].status!=BORDER){
				if(isin){
					partitions[x][y].status = IN;
				}
				continue;
			}
			for(cross_info &c:partitions[x][y].crosses){
				if(c.direction==BOTTOM){
					isin = !isin;
				}
			}
		}
	}

	pthread_mutex_unlock(&partition_lock);
	return partitions;
}


vector<vector<Pixel>> MyPolygon::partition_with_query(int vpr){
	assert(vpr>0);
	pthread_mutex_lock(&partition_lock);
	if(is_grid_partitioned()){
		pthread_mutex_unlock(&partition_lock);
		return partitions;
	}
	// normalize
	getMBB();
	const double start_x = mbb->low[0];
	const double start_y = mbb->low[1];
	const double end_x = mbb->high[0];
	const double end_y = mbb->high[1];

	double multi = abs((end_y-start_y)/(end_x-start_x));
	int dimx = std::pow((get_num_vertices()/vpr)/multi,0.5);
	int dimy = multi*dimx;


	step_x = (end_x-start_x)/dimx;
	step_y = (end_y-start_y)/dimy;

	if(step_x<0.00001){
		step_x = 0.00001;
		dimx = (end_x-start_x)/step_x+1;
	}
	if(step_y<0.00001){
		step_y = 0.00001;
		dimy = (end_y-start_y)/step_y+1;
	}

	for(double i=0;i<=dimx;i++){
		vector<Pixel> v;
		for(double j=0;j<=dimy;j++){
			Pixel m;
			m.id[0] = i;
			m.id[1] = j;
			m.low[0] = i*step_x+start_x;
			m.high[0] = (i+1.0)*step_x+start_x;
			m.low[1] = j*step_y+start_y;
			m.high[1] = (j+1.0)*step_y+start_y;
			v.push_back(m);
			Point p1(m.low[0],m.low[1]);
			Point p2(m.high[0],m.high[1]);
			bool in1 = contain(p1,NULL);
			bool in2 = contain(p2,NULL);;
			if(in1&&in2){
				m.status = IN;
			}else if(!in1&&!in2){
				m.status = OUT;
			}else{
				m.status = BORDER;
			}
		}
		partitions.push_back(v);
	}

	pthread_mutex_unlock(&partition_lock);
	return partitions;
}


QTNode *MyPolygon::partition_qtree(const int vpr){

	pthread_mutex_lock(&partition_lock);
	if(this->is_qtree_partitioned()){
		pthread_mutex_unlock(&partition_lock);
		return qtree;
	}
	int num_boxes = get_num_vertices()/vpr;
	int level = 1;
	while(pow(4,level)<num_boxes){
		level++;
	}
	int dimx = pow(2,level);
	int dimy = dimx;

	this->partition(dimx, dimy);
	int box_count = 4;
	int cur_level = 1;

	qtree = new QTNode(*(this->getMBB()));
	std::stack<QTNode *> ws;
	qtree->split();
	qtree->push(ws);

	vector<QTNode *> level_nodes;

	//breadth first traverse
	query_context qc;
	while(box_count<num_boxes||!ws.empty()){

		if(cur_level>level){
			dimx *= 2;
			dimy *= 2;
			level = cur_level;
			this->reset_grid_partition();
			this->partition(dimx, dimy);
		}

		while(!ws.empty()){
			QTNode *cur = ws.top();
			ws.pop();
			bool inside = this->contain_try_partition(&cur->mbb, &qc);

			if(!qc.raster_checked_only){
				level_nodes.push_back(cur);
			}else if(inside){
				cur->interior = true;
			}else{
				cur->exterior = true;
			}
		}


		for(QTNode *n:level_nodes){
			if(box_count<num_boxes){
				n->split();
				n->push(ws);
				box_count += 3;
			}else{
				break;
			}
		}
		cur_level++;
		level_nodes.clear();
	}
	//assert(box_count>=num_boxes);

	pthread_mutex_unlock(&partition_lock);

	return qtree;
}



vector<vector<Pixel>> MyPolygon::decode_partition(char *data){
	assert(data);
	vector<vector<Pixel>> partitions;

	int dimx = data[0];
	int dimy = data[1];
	if(dimx==0||dimy==0||dimx>255||dimy>255){
		return partitions;
	}
	double *meta = (double *)(data+2);
	double start_x = meta[0], start_y = meta[1], step_x = meta[2], step_y = meta[3];

	int index = 2+4*8;
	for(int i=0;i<dimx;i++){
		char cur = data[index++];
		int shift = 0;
		vector<Pixel> row;
		for(int j=0;j<dimy;j++){
			Pixel p;

			p.status = (PartitionStatus)((cur>>shift)&3);
			p.id[0] = i;
			p.id[1] = j;
			p.low[0] = i*step_x+start_x;
			p.low[1] = j*step_y+start_y;
			p.high[0] = (i+1)*step_x+start_x;
			p.high[1] = (j+1)*step_y+start_y;
			row.push_back(p);
			if(shift==6){
				cur = data[index++];
				shift=0;
			}else{
				shift+=2;
			}
		}
		partitions.push_back(row);
	}

	return partitions;
}

char *MyPolygon::encode_partition(vector<vector<Pixel>> partitions){
	int dimx = partitions.size();
	if(dimx==0){
		return NULL;
	}
	int dimy = partitions[0].size();
	if(dimy==0){
		return NULL;
	}
	for(vector<Pixel> &rows:partitions){
		 if(rows.size()!=dimy){
			 return NULL;
		 }
	}
	assert(dimx<=255&&dimy<=255);

	char *data = new char[((dimy+3)/4*dimx)+2+4*8];
	data[0] = (char)dimx;
	data[1] = (char)dimy;
	double *meta = (double *)(data+2);
	meta[0] = partitions[0][0].low[0];
	meta[1] = partitions[0][0].low[1];
	meta[2] = partitions[0][0].high[0]-partitions[0][0].low[0];
	meta[3] = partitions[0][0].high[1]-partitions[0][0].low[1];
	int index = 2+4*8;
	for(int i=0;i<dimx;i++){
		char cur = 0;
		int shift = 0;
		for(int j=0;j<dimy;j++){
			cur |= ((uint8_t)partitions[i][j].status)<<shift;
			if(shift==6){
				data[index++]=cur;
				cur = 0;
				shift=0;
			}else{
				shift+=2;
			}
		}
		if(shift>0){
			data[index++]=cur;
		}
	}
	return data;
}

void *partition_unit(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	//log("thread %d is started",ctx->thread_id);
	int local_count = 0;
	while(ctx->next_batch(1)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			struct timeval start = get_cur_time();

			if(ctx->use_grid){
				gctx->source_polygons[i]->partition(ctx->vpr);
			}
			if(ctx->use_qtree){
				gctx->source_polygons[i]->partition_qtree(ctx->vpr);
			}
			double latency = get_time_elapsed(start);
			int num_vertices = gctx->source_polygons[i]->get_num_vertices();
			//ctx->report_latency(num_vertices, latency);
			if(latency>10000||num_vertices>200000){
				logt("partition %d vertices",start,num_vertices);
			}
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}

void process_partition(query_context *gctx){

	gctx->index = 0;
	size_t former = gctx->target_num;
	gctx->target_num = gctx->source_polygons.size();
	struct timeval start = get_cur_time();
	pthread_t threads[gctx->num_threads];
	query_context ctx[gctx->num_threads];
	for(int i=0;i<gctx->num_threads;i++){
		ctx[i] = *gctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = gctx;
	}

	for(int i=0;i<gctx->num_threads;i++){
		pthread_create(&threads[i], NULL, partition_unit, (void *)&ctx[i]);
	}

	for(int i = 0; i < gctx->num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	//collect partitioning status
	size_t num_partitions = 0;
	for(MyPolygon *poly:gctx->source_polygons){
		num_partitions += poly->get_num_partitions();
	}

	logt("partitioned %d polygons with %ld average partitions", start,
			gctx->source_polygons.size(),
			num_partitions/gctx->source_polygons.size());
	gctx->index = 0;
	gctx->query_count = 0;
	gctx->target_num = former;
}
