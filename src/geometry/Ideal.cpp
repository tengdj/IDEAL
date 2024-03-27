#include "../include/Ideal.h"

void Ideal::add_edge(int idx, int start, int end){
	edge_sequences[idx] = make_pair(start, end - start  + 1);
}

void Ideal::init_edge_sequences(int num_edge_seqs){
	len_edge_sequences = num_edge_seqs;
	edge_sequences = new pair<uint16_t, uint16_t>[num_edge_seqs];
}

void Ideal::process_pixels_null(int x, int y){
	offset[(x+1)*(y+1)] = len_edge_sequences;
	for(int i = (x+1)*(y+1)-1; i >= 0; i --){
		if(show_status(i) != BORDER){
			offset[i] = offset[i + 1]; 
		}
	}
}

void Ideal::process_crosses(map<int, vector<cross_info>> edges_info){
	int num_edge_seqs = 0;
	for(auto ei : edges_info){
		num_edge_seqs += ei.second.size();
	}
	init_edge_sequences(num_edge_seqs);

	int idx = 0;
	int edge_count = 0;
	for(auto info : edges_info){
		auto pix = info.first;
		auto crosses = info.second; 
		if(crosses.size() == 0) return;
		
		if(crosses.size() % 2 == 1){
			crosses.push_back(cross_info((cross_type)!crosses[crosses.size()-1].type, crosses[crosses.size()-1].edge_id));
		}
		
		assert(crosses.size()%2==0);

		// Initialize based on crosses.size().
		int start = 0;
		int end = crosses.size() - 1;
		set_offset(pix, idx);

		if(crosses[0].type == LEAVE){
			assert(crosses[end].type == ENTER);
			add_edge(idx ++, 0, crosses[0].edge_id);
			add_edge(idx ++, crosses[end].edge_id, boundary->num_vertices - 2);
			start ++;
			end --;
		}

		for(int i = start; i <= end; i++){
			assert(crosses[i].type == ENTER);
			//special case, an ENTER has no pair LEAVE,
			//happens when one edge crosses the pair
			if(i == end || crosses[i + 1].type == ENTER){
				add_edge(idx ++, crosses[i].edge_id, crosses[i].edge_id);
			}else{
				add_edge(idx ++, crosses[i].edge_id, crosses[i+1].edge_id);
				i++;
			}
		}
	}
}

void Ideal::process_intersection(map<int, vector<double>> intersection_info, Direction direction){
	int num_nodes = 0;
	for(auto i : intersection_info){
		num_nodes += i.second.size();
	}
	if(direction == HORIZONTAL){
		horizontal->init_intersection_node(num_nodes);
		horizontal->set_num_crosses(num_nodes);
		int idx = 0;
		for(auto info : intersection_info){
			auto h = info.first;
			auto nodes = info.second;
			
			sort(nodes.begin(), nodes.end());

			horizontal->set_offset(h, idx);

			for(auto node : nodes){
				horizontal->add_node(idx, node);
				idx ++;
			}
		}
		horizontal->set_offset(dimy, idx);
	}else{
		vertical->init_intersection_node(num_nodes);
		vertical->set_num_crosses(num_nodes);

		int idx = 0;
		for(auto info : intersection_info){
			auto h = info.first;
			auto nodes = info.second;
			
			sort(nodes.begin(), nodes.end());

			vertical->set_offset(h, idx);

			for(auto node : nodes){
				vertical->add_node(idx, node);
				idx ++;
			}
		}
		vertical->set_offset(dimx, idx);		
	}

}

void Ideal::init_pixels(){
	assert(mbr);
    offset = new uint16_t[(dimx+1)*(dimy+1) + 1];    // +1 here is to ensure that pointer[num_pixels] equals len_edge_sequences, so we don't need to make a special case for the last pointer.
	horizontal = new Grid_line(dimy);
	vertical = new Grid_line(dimx);
}

void Ideal::evaluate_edges(){
	map<int, vector<double>> horizontal_intersect_info;
	map<int, vector<double>> vertical_intersect_info;
	map<int, vector<cross_info>> edges_info;
	
	// normalize
	assert(mbr);
	const double start_x = mbr->low[0];
	const double start_y = mbr->low[1];

	for(int i=0;i<boundary->num_vertices-1;i++){
		double x1 = boundary->p[i].x;
		double y1 = boundary->p[i].y;
		double x2 = boundary->p[i+1].x;
		double y2 = boundary->p[i+1].y;

		int cur_startx = (x1-start_x)/step_x;
		int cur_endx = (x2-start_x)/step_x;
		int cur_starty = (y1-start_y)/step_y;
		int cur_endy = (y2-start_y)/step_y;

		if(cur_startx==dimx+1){
			cur_startx--;
		}
		if(cur_endx==dimx+1){
			cur_endx--;
		}

		int minx = min(cur_startx,cur_endx);
		int maxx = max(cur_startx,cur_endx);

		if(cur_starty==dimy+1){
			cur_starty--;
		}
		if(cur_endy==dimy+1){
			cur_endy--;
		}
		// todo should not happen for normal cases
		if(cur_startx>dimx||cur_endx>dimx||cur_starty>dimy||cur_endy>dimy){
			cout<<"xrange\t"<<cur_startx<<" "<<cur_endx<<endl;
			cout<<"yrange\t"<<cur_starty<<" "<<cur_endy<<endl;
			printf("xrange_val\t%f %f\n",(x1-start_x)/step_x, (x2-start_x)/step_x);
			printf("yrange_val\t%f %f\n",(y1-start_y)/step_y, (y2-start_y)/step_y);
			assert(false);
		}
		assert(cur_startx<=dimx);
		assert(cur_endx<=dimx);
		assert(cur_starty<=dimy);
		assert(cur_endy<=dimy);

		//in the same pixel
		if(cur_startx==cur_endx&&cur_starty==cur_endy){
			continue;
		}

		if(y1==y2){
			//left to right
			if(cur_startx<cur_endx){
				for(int x=cur_startx;x<cur_endx;x++){
					vertical_intersect_info[x + 1].push_back(y1);
					edges_info[get_id(x, cur_starty)].push_back(cross_info(LEAVE, i));
					edges_info[get_id(x + 1, cur_starty)].push_back(cross_info(ENTER, i));
					set_status(get_id(x, cur_starty), BORDER);
					set_status(get_id(x+1, cur_starty), BORDER);
				}
			}else { // right to left
				for(int x=cur_startx;x>cur_endx;x--){
					vertical_intersect_info[x].push_back(y1);
					edges_info[get_id(x, cur_starty)].push_back(cross_info(LEAVE, i));
					edges_info[get_id(x - 1, cur_starty)].push_back(cross_info(ENTER, i));
					set_status(get_id(x, cur_starty), BORDER);
					set_status(get_id(x-1, cur_starty), BORDER);
				}
			}
		}else if(x1==x2){
			//bottom up
			if(cur_starty<cur_endy){
				for(int y=cur_starty;y<cur_endy;y++){
					horizontal_intersect_info[y + 1].push_back(x1);
					edges_info[get_id(cur_startx, y)].push_back(cross_info(LEAVE, i));
					edges_info[get_id(cur_startx, y + 1)].push_back(cross_info(ENTER, i));
					set_status(get_id(cur_startx, y), BORDER);
					set_status(get_id(cur_startx, y+1), BORDER);
				}
			}else { //border[bottom] down
				for(int y=cur_starty;y>cur_endy;y--){
					horizontal_intersect_info[y].push_back(x1);
					edges_info[get_id(cur_startx, y)].push_back(cross_info(LEAVE, i));
					edges_info[get_id(cur_startx, y - 1)].push_back(cross_info(ENTER, i));
					set_status(get_id(cur_startx, y), BORDER);
					set_status(get_id(cur_startx, y-1), BORDER);
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
							vertical_intersect_info[x + 1].push_back(yval);
							set_status(get_id(x, y), BORDER);
							edges_info[get_id(x ++, y)].push_back(cross_info(LEAVE, i));
							edges_info[get_id(x, y)].push_back(cross_info(ENTER, i));
							set_status(get_id(x, y), BORDER);
						}else{//right to left
							vertical_intersect_info[x].push_back(yval);
							set_status(get_id(x, y), BORDER);
							edges_info[get_id(x --, y)].push_back(cross_info(LEAVE, i));
							edges_info[get_id(x, y)].push_back(cross_info(ENTER, i));
							set_status(get_id(x, y), BORDER);
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
							horizontal_intersect_info[y + 1].push_back(xval);
							set_status(get_id(x, y), BORDER);
							edges_info[get_id(x, y ++)].push_back(cross_info(LEAVE, i));
							edges_info[get_id(x, y)].push_back(cross_info(ENTER, i));
							set_status(get_id(x, y), BORDER);
						}else{// top down
							horizontal_intersect_info[y].push_back(xval);
							set_status(get_id(x, y), BORDER);
							edges_info[get_id(x, y --)].push_back(cross_info(LEAVE, i));
							edges_info[get_id(x, y)].push_back(cross_info(ENTER, i));
							set_status(get_id(x, y), BORDER);
						}
					}
				}
				// for debugging, should never happen
				if(!passed){
					boundary->print();
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


	process_crosses(edges_info);
	process_intersection(horizontal_intersect_info, HORIZONTAL);
	process_intersection(vertical_intersect_info, VERTICAL);
	process_pixels_null(dimx, dimy);
}

void Ideal::scanline_reandering(){
	const double start_x = mbr->low[0];
	const double start_y = mbr->low[1];

	for(int y = 1; y < dimy; y ++){
		bool isin = false;
		uint16_t i = horizontal->get_offset(y), j = horizontal->get_offset(y + 1);
		for(int x = 0; x < dimx; x ++){
			if(show_status(get_id(x, y)) != BORDER){
				if(isin){
					set_status(get_id(x, y), IN);
				}else{
					set_status(get_id(x, y), OUT);
				}
				continue;
			}
			int pass = 0;
			while(i < j && horizontal->get_intersection_nodes(i) <= start_x + step_x * (x + 1)){
				pass ++;
				i ++;
			}
			if(pass % 2 == 1) isin = !isin;

		}
	}
}

void Ideal::rasterization(){

	//1. create space for the pixels
	init_pixels();

	//2. edge crossing to identify BORDER pixels
	evaluate_edges();

	//3. determine the status of rest pixels with scanline rendering
	scanline_reandering();
}

void Ideal::rasterization(int vpr){
	assert(vpr > 0);
	pthread_mutex_lock(&ideal_partition_lock);
    init_raster(boundary->num_vertices / vpr);
    rasterization();
	pthread_mutex_unlock(&ideal_partition_lock);
}