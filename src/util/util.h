#ifndef TENG_UTIL_H_
#define TENG_UTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <vector>
#include <thread>
#include <iostream>
#include <stdarg.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/syscall.h>
#include <time.h>
#include <sstream>
#include <math.h>

using namespace std;

namespace{

#define TENG_RANDOM_NUMBER 0315
#define OSM_SRID 4326

// some utility function

const double degree_per_kilometer_latitude = 360.0/40076.0;

inline double degree_per_kilometer_longitude(double latitude){
	double absla = abs(latitude);
	if(absla==90){
		absla = 89;
	}
	assert(absla<90);
	return 360.0/(sin((90-absla)/90)*40076);
}

inline int double_to_int(double val){
	int vi = (int)val;
	if(abs(1.0*(vi+1)-val)<0.00000001){
		vi++;
	}
	return vi;
}

inline bool is_number(char ch){
	return ch=='-'||ch=='.'||(ch<='9'&&ch>='0')||ch=='e';
}

inline double read_double(const char *input, size_t &offset){
	char tmp[100];
	while(!is_number(input[offset])){
		offset++;
	}
	int index = 0;
	while(is_number(input[offset])){
		tmp[index++] = input[offset++];
	}
	tmp[index] = '\0';
	return atof(tmp);
}
inline void skip_space(const char *input, size_t &offset){
	while(input[offset]==' '||input[offset]=='\t'||input[offset]=='\n'){
		offset++;
	}
}

inline struct timeval get_cur_time(){
	struct timeval t1;
	gettimeofday(&t1, NULL);
	return t1;
}
inline double get_time_elapsed(struct timeval &t1, bool update_start = false){
	struct timeval t2;
    double elapsedTime;
	gettimeofday(&t2, NULL);
	// compute and print the elapsed time in millisec
	elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
	elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
	if(update_start){
		t1 = get_cur_time();
	}
	return elapsedTime;
}

inline string time_string(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	struct tm *nowtm;
	char tmbuf[64];
	char buf[64];
	nowtm = localtime(&tv.tv_sec);
	strftime(tmbuf, sizeof tmbuf, "%H:%M:%S", nowtm);
	sprintf(buf,"%s.%04ld", tmbuf, tv.tv_usec/1000);
	return string(buf);
}

static pthread_mutex_t print_lock;
inline void logt(const char *format, struct timeval &start, ...){
	pthread_mutex_lock(&print_lock);
	va_list args;
	va_start(args, start);
	char sprint_buf[200];
	int n = vsprintf(sprint_buf, format, args);
	va_end(args);
	fprintf(stderr,"%s thread %ld:\t%s", time_string().c_str(), syscall(__NR_gettid),sprint_buf);

	double mstime = get_time_elapsed(start, true);
	if(mstime>1000){
		fprintf(stderr," takes %f s\n", mstime/1000);
	}else{
		fprintf(stderr," takes %f ms\n", mstime);
	}
	fflush(stderr);

	pthread_mutex_unlock(&print_lock);
}

inline void log(const char *format, ...){
	pthread_mutex_lock(&print_lock);
	va_list args;
	va_start(args, format);
	char sprint_buf[200];
	int n = vsprintf(sprint_buf, format, args);
	va_end(args);
	fprintf(stderr,"%s thread %ld:\t%s\n", time_string().c_str(), syscall(__NR_gettid),sprint_buf);
	fflush(stderr);
	pthread_mutex_unlock(&print_lock);
}

inline void log_stdout(const char *format, ...){
	pthread_mutex_lock(&print_lock);
	va_list args;
	va_start(args, format);
	char sprint_buf[200];
	int n = vsprintf(sprint_buf, format, args);
	va_end(args);
	fprintf(stdout,"%s thread %ld:\t%s\n", time_string().c_str(), syscall(__NR_gettid),sprint_buf);
	fflush(stdout);
	pthread_mutex_unlock(&print_lock);
}

inline void log(){
	pthread_mutex_lock(&print_lock);
	fprintf(stdout,"%s thread %ld:\tterry is good\n", time_string().c_str(),syscall(__NR_gettid));
	fflush(stdout);
	pthread_mutex_unlock(&print_lock);
}

inline int get_rand_number(int max_value){
	return rand()%max_value+1;
}

inline double get_rand_double(){
	return rand()/(double)RAND_MAX;
}

inline bool get_rand_sample(int rate){
	return rand()%100<rate;
}

inline bool tryluck(float possibility){
	assert(possibility<=1&&possibility>=0);
	return possibility>=1.0||(rand()*1.0)/RAND_MAX<possibility;
}

inline bool is_dir(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISDIR(buf.st_mode);
}

inline bool is_file(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISREG(buf.st_mode);
}

inline void list_files(const char *path, std::vector<string> &f_list){
	if(is_file(path)){
		f_list.push_back(std::string(path));
		return;
	}
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir (path)) != NULL) {
		/* print all the files and directories within directory */
		while ((ent = readdir (dir)) != NULL) {
			if(strcmp(ent->d_name,"..")==0||
					strcmp(ent->d_name,".")==0){
				continue;
			}
			std::stringstream ss;
			ss<<path<<"/"<<ent->d_name;
			string spath;
			ss>>spath;
			list_files(spath.c_str(), f_list);
		}
		closedir (dir);
	}
}

inline long file_size(const char *file){
	struct stat stat_buf;
	int rc = stat(file, &stat_buf);
	return rc == 0 ? stat_buf.st_size : -1;
}

inline long file_size(std::vector<string> &f_list){
	long size = 0;
	for(string s:f_list){
		long ls = file_size(s.c_str());
		if(ls>0){
			size += ls;
		}
	}
	return size;
}

inline bool file_exist(const char *path) {
  struct stat buffer;
  return (stat(path, &buffer) == 0);
}

inline int get_num_threads(){
	return std::thread::hardware_concurrency();
}

inline string read_line(){
	string input_line;
	getline(std::cin, input_line);
	return input_line;
}


inline void tokenize( const std::string& str, std::vector<std::string>& result,
	const std::string& delimiters = " ,;:\t",
	const bool keepBlankFields=false,
	const std::string& quote="\"\'"
	){
    // clear the vector
    if (!result.empty()){
    	result.clear();
    }

    // you must be kidding
    if (delimiters.empty())
	return ;

    std::string::size_type pos = 0; // the current position (char) in the string
    char ch = 0; // buffer for the current character

    char current_quote = 0; // the char of the current open quote
    bool quoted = false; // indicator if there is an open quote
    std::string token;  // string buffer for the token
    bool token_complete = false; // indicates if the current token is
    // read to be added to the result vector
    std::string::size_type len = str.length();  // length of the input-string

    // for every char in the input-string
	while(len > pos){
		// get the character of the string and reset the delimiter buffer
		ch = str.at(pos);

		bool add_char = true;
		if ( false == quote.empty()){
			// if quote chars are provided and the char isn't protected
			if (std::string::npos != quote.find_first_of(ch)){
				if (!quoted){
					quoted = true;
					current_quote = ch;
					add_char = false;
				} else {
					if (current_quote == ch){
						quoted = false;
						current_quote = 0;
						add_char = false;
					}
				}
			}
		}

		if (!delimiters.empty()&&!quoted){
			// if ch is delemiter
			if (std::string::npos != delimiters.find_first_of(ch)){
				token_complete = true;
				// don't add the delimiter to the token
				add_char = false;
			}
		}

		// add the character to the token
		if (add_char){
			token.push_back(ch);
		}

		// add the token if it is complete
		// if ( true == token_complete && false == token.empty() )
		if (token_complete){
			if (token.empty())
			{
			if (keepBlankFields)
				result.push_back("");
			}
			else
			result.push_back( token );
			token.clear();
			token_complete = false;
		}
		++pos;
    } // while
    // add the final token
    if ( false == token.empty() ) {
    	result.push_back( token );
    } else if(keepBlankFields && std::string::npos != delimiters.find_first_of(ch) ){
    	result.push_back("");
    }
}


inline void remove_slash(string &str){
	if(str.at(str.size() - 1) == '/'){
		str = str.substr(0, str.size() - 1);
	}
}

}
#endif
