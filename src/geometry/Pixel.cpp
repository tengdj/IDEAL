/*
 * Pixel.cpp
 *
 *  Created on: May 26, 2020
 *      Author: teng
 */


#include "MyPolygon.h"

bool print_debug = false;


bool Pixel::overwrites(cross_info &enter1, cross_info &leave1, cross_info &enter2, cross_info &leave2){

	if(enter2.direction==LEFT&&leave2.direction==LEFT){
		if(enter2.vertex<leave2.vertex){
			if(enter1.direction==LEFT&&enter1.vertex>leave2.vertex){
				if(leave1.direction==LEFT){
					if(leave1.vertex<enter2.vertex){
						return true;
					}
				}else{
					return true;
				}
			}
			if(enter1.direction==RIGHT&&leave1.direction==LEFT&&leave1.vertex<enter2.vertex){
				return true;
			}
			if(enter1.direction==TOP&&leave1.direction==BOTTOM){
				return true;
			}
			if(enter1.direction==TOP&&leave1.direction==LEFT&&leave1.vertex<enter2.vertex){
				return true;
			}
		}
		return false;
	}else if(enter2.direction==LEFT&&leave2.direction==TOP){
		// right bottom
		if(enter1.direction==RIGHT&&leave1.direction==LEFT&&
				leave1.vertex<enter2.vertex){
			return true;
		}
		if(enter1.direction==TOP&&leave1.direction==LEFT&&
				leave1.vertex<enter2.vertex&&enter1.vertex>leave2.vertex){
			return true;
		}
		if(enter1.direction==TOP&&leave1.direction==BOTTOM&&
				enter1.vertex>leave2.vertex){
			return true;
		}
		return false;
	}else if(enter2.direction==LEFT&&leave2.direction==RIGHT){
		// bottom
		if(enter1.direction==RIGHT&&leave1.direction==LEFT&&
				enter1.vertex<leave2.vertex&&leave1.vertex<enter2.vertex){
			return true;
		}
		return false;

	}else if(enter2.direction==LEFT&&leave2.direction==BOTTOM){
		return false;
	}else if(enter2.direction==TOP&&leave2.direction==LEFT){
		return false;
	}else if(enter2.direction==TOP&&leave2.direction==TOP){
		if(enter2.vertex<leave2.vertex){
			if(enter1.direction==TOP&&enter1.vertex>leave2.vertex){
				if(leave1.direction==TOP){
					if(leave1.vertex<enter2.vertex){
						return true;
					}
				}else{
					return true;
				}
			}
			if(enter1.direction==BOTTOM&&leave1.direction==TOP&&leave1.vertex<enter2.vertex){
				return true;
			}
			if(enter1.direction==RIGHT&&leave1.direction==LEFT){
				return true;
			}
			if(enter1.direction==RIGHT&&leave1.direction==TOP&&leave1.vertex<enter2.vertex){
				return true;
			}
		}
		return false;
	}else if(enter2.direction==TOP&&leave2.direction==RIGHT){
		//left bottom
		if(enter1.direction==RIGHT&&leave1.direction==TOP&&
				leave1.vertex<enter2.vertex&&enter1.vertex<leave2.vertex){
			return true;
		}
		if(enter1.direction==RIGHT&&leave1.direction==LEFT&&
				enter1.vertex<leave2.vertex){
			return true;
		}
		if(enter1.direction==BOTTOM&&leave1.direction==TOP&&
				leave1.vertex<enter2.vertex){
			return true;
		}
		return false;
	}else if(enter2.direction==TOP&&leave2.direction==BOTTOM){
		// left
		if(enter1.direction==BOTTOM&&leave1.direction==TOP&&
				enter1.vertex<leave2.vertex&&leave1.vertex<enter2.vertex){
			return true;
		}
		return false;
	}else if(enter2.direction==RIGHT&&leave2.direction==LEFT){
		// top
		if(enter1.direction==LEFT&&leave1.direction==RIGHT&&
				enter1.vertex>leave2.vertex&&leave1.vertex>enter2.vertex){
			return true;
		}
		return false;
	}else if(enter2.direction==RIGHT&&leave2.direction==TOP){
		return false;
	}else if(enter2.direction==RIGHT&&leave2.direction==RIGHT){
		if(enter2.vertex>leave2.vertex){
			if(enter1.direction==RIGHT&&enter1.vertex<leave2.vertex){
				if(leave1.direction==RIGHT){
					if(leave1.vertex>enter2.vertex){
						return true;
					}
				}else{
					return true;
				}
			}
			if(enter1.direction==LEFT&&leave1.direction==RIGHT&&leave1.vertex>enter2.vertex){
				return true;
			}
			if(enter1.direction==BOTTOM&&leave1.direction==TOP){
				return true;
			}
			if(enter1.direction==BOTTOM&&leave1.direction==RIGHT&&leave1.vertex>enter2.vertex){
				return true;
			}
		}
		return false;
	}else if(enter2.direction==RIGHT&&leave2.direction==BOTTOM){
		// left top
		if(enter1.direction==BOTTOM&&leave1.direction==RIGHT&&
				leave1.vertex>enter2.vertex&&enter1.vertex<leave2.vertex){
			return true;
		}
		if(enter1.direction==BOTTOM&&leave1.direction==TOP&&
				enter1.vertex<leave2.vertex){
			return true;
		}
		if(enter1.direction==LEFT&&leave1.direction==RIGHT&&
				leave1.vertex>enter2.vertex){
			return true;
		}
		return false;
	}else if(enter2.direction==BOTTOM&&leave2.direction==LEFT){
		// right top
		if(enter1.direction==LEFT&&leave1.direction==BOTTOM&&
				leave1.vertex>enter2.vertex&&enter1.vertex>leave2.vertex){
			return true;
		}
		if(enter1.direction==LEFT&&leave1.direction==RIGHT&&
				enter1.vertex>leave2.vertex){
			return true;
		}
		if(enter1.direction==TOP&&leave1.direction==BOTTOM&&
				leave1.vertex>enter2.vertex){
			return true;
		}
		return false;
	}else if(enter2.direction==BOTTOM&&leave2.direction==TOP){
		// right
		if(enter1.direction==TOP&&leave1.direction==BOTTOM&&
				enter1.vertex>leave2.vertex&&leave1.vertex>enter2.vertex){
			return true;
		}
		return false;
	}else if(enter2.direction==BOTTOM&&leave2.direction==RIGHT){
		return false;
	}else if(enter2.direction==BOTTOM&&leave2.direction==BOTTOM){
		if(enter2.vertex>leave2.vertex){
			if(enter1.direction==BOTTOM&&enter1.vertex<leave2.vertex){
				if(leave1.direction==BOTTOM){
					if(leave1.vertex>enter2.vertex){
						return true;
					}
				}else{
					return true;
				}
			}
			if(enter1.direction==TOP&&leave1.direction==BOTTOM&&leave1.vertex>enter2.vertex){
				return true;
			}
			if(enter1.direction==LEFT&&leave1.direction==RIGHT){
				return true;
			}
			if(enter1.direction==LEFT&&leave1.direction==BOTTOM&&leave1.vertex>enter2.vertex){
				return true;
			}
		}
		return false;
	}else{
		assert(false);
		return false;
	}

}

void Pixel::process_enter_leave(){
	//if(id[0]==16&&id[1]==18)
	if(crosses.size()>0&&print_debug)
	{
		cout<<id[0]<<" "<<id[1]<<" "<<crosses.size()<<endl;
		for(cross_info &ci:crosses){
			cout<<" "<<direction_str[ci.direction];
		}
		cout<<endl;
		for(cross_info &ci:crosses){
			printf(" %f",ci.vertex);
		}
		cout<<endl;

		if(status==BORDER){
			for(int k=0;k<4;k++){
				if(border[k]==IN){
					cout<<" "<<direction_str[k];
				}
			}
			cout<<endl;
		}
	}
	if(crosses.size()==0){
		return;
	}
	assert(crosses.size()%2==0);
	// the first pixel is left first and then entered last.
	if(crosses[0].type==LEAVE){
		cross_info ci = crosses[crosses.size()-1];
		assert(ci.type==ENTER);
		crosses.insert(crosses.begin(), ci);
		crosses.pop_back();
	}
	// in the first round, mark all the crossed boundary as
	// BORDER edge
	for(int i=0;i<crosses.size()-1;i+=2){
		Direction enter_d= crosses[i].direction;
		Direction leave_d = crosses[i+1].direction;
		double enter_val = crosses[i].vertex;
		double leave_val = crosses[i+1].vertex;
		border[enter_d] = BORDER;
		border[leave_d] = BORDER;
	}
	// determining if any boundary is inside the polygon
	for(int i=0;i<crosses.size()-1;i+=2){
		Direction enter_d= crosses[i].direction;
		Direction leave_d = crosses[i+1].direction;
		double enter_val = crosses[i].vertex;
		double leave_val = crosses[i+1].vertex;

		// check
		bool overwritten = false;
		for(int j=0;j<crosses.size()-1;j+=2){
			if(j==i){
				continue;
			}
			if(overwrites(crosses[j],crosses[j+1],crosses[i],crosses[i+1])){
				overwritten = true;
				break;
			}
		}
		if(overwritten){
			continue;
		}

		// totally 16 possible combinations
		if(enter_d==LEFT&&leave_d==LEFT){
			if(leave_val>enter_val){
				if(border[RIGHT]==OUT){
					border[RIGHT]=IN;
				}
				if(border[TOP]==OUT){
					border[TOP]=IN;
				}
				if(border[BOTTOM]==OUT){
					border[BOTTOM]=IN;
				}
			}
		}else if(enter_d==LEFT&&leave_d==TOP){
			if(border[RIGHT]==OUT){
				border[RIGHT]=IN;
			}
			if(border[BOTTOM]==OUT){
				border[BOTTOM]=IN;
			}
		}else if(enter_d==LEFT&&leave_d==RIGHT){
			if(border[BOTTOM]==OUT){
				border[BOTTOM]=IN;
			}
		}else if(enter_d==LEFT&&leave_d==BOTTOM){

		}else if(enter_d==TOP&&leave_d==LEFT){

		}else if(enter_d==TOP&&leave_d==TOP){
			if(leave_val>enter_val){
				if(border[RIGHT]==OUT){
					border[RIGHT]=IN;
				}
				if(border[LEFT]==OUT){
					border[LEFT]=IN;
				}
				if(border[BOTTOM]==OUT){
					border[BOTTOM]=IN;
				}
			}
		}else if(enter_d==TOP&&leave_d==RIGHT){
			if(border[LEFT]==OUT){
				border[LEFT]=IN;
			}
			if(border[BOTTOM]==OUT){
				border[BOTTOM]=IN;
			}
		}else if(enter_d==TOP&&leave_d==BOTTOM){
			if(border[LEFT]==OUT){
				border[LEFT]=IN;
			}
		}else if(enter_d==RIGHT&&leave_d==LEFT){
			if(border[TOP]==OUT){
				border[TOP]=IN;
			}
		}else if(enter_d==RIGHT&&leave_d==TOP){

		}else if(enter_d==RIGHT&&leave_d==RIGHT){
			if(leave_val<enter_val){
				if(border[TOP]==OUT){
					border[TOP]=IN;
				}
				if(border[LEFT]==OUT){
					border[LEFT]=IN;
				}
				if(border[BOTTOM]==OUT){
					border[BOTTOM]=IN;
				}
			}
		}else if(enter_d==RIGHT&&leave_d==BOTTOM){
			if(border[LEFT]==OUT){
				border[LEFT]=IN;
			}
			if(border[TOP]==OUT){
				border[TOP]=IN;
			}
		}else if(enter_d==BOTTOM&&leave_d==LEFT){
			if(border[RIGHT]==OUT){
				border[RIGHT]=IN;
			}
			if(border[TOP]==OUT){
				border[TOP]=IN;
			}
		}else if(enter_d==BOTTOM&&leave_d==TOP){
			if(border[RIGHT]==OUT){
				border[RIGHT]=IN;
			}
		}else if(enter_d==BOTTOM&&leave_d==RIGHT){

		}else if(enter_d==BOTTOM&&leave_d==BOTTOM){
			if(leave_val<enter_val){
				if(border[TOP]==OUT){
					border[TOP]=IN;
				}
				if(border[LEFT]==OUT){
					border[LEFT]=IN;
				}
				if(border[RIGHT]==OUT){
					border[RIGHT]=IN;
				}
			}
		}else{
			assert(false);
		}
	}
	this->status = BORDER;
}


void Pixel::enter(double val, Direction d){
	if(print_debug){
		cout<<direction_str[d];
		cout<<" enter "<<id[0]<<" "<<id[1]<<endl;
	}
	cross_info ci;
	ci.type = ENTER;
	ci.vertex = val;
	ci.direction = d;
	crosses.push_back(ci);
}



void Pixel::leave(double val, Direction d){
	if(print_debug){
		cout<<direction_str[d];
		cout<<" leave "<<id[0]<<" "<<id[1]<<endl;
	}
	cross_info ci;
	ci.type = LEAVE;
	ci.vertex = val;
	ci.direction = d;
	crosses.push_back(ci);
}


MyPolygon *Pixel::to_polygon(){
	return 	MyPolygon::gen_box(low[0],low[1],high[0],high[1]);
};

