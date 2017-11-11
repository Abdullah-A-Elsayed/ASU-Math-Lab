#include <iostream>
#include <string>
#include <fstream>
#include "matrix.h"
using namespace std;
int main(int argc, char** argv)
{
	if(argc < 2){
		//matrix::run("example.m");
		//to be changed ...
		cout<<"please enter your commands line by line"<<endl;
		cout<<"enter 'run' to see results."<<endl;
		string command;
		getline(cin,command);
		ofstream file ("temp.m");
		while(command!="run"){
			file<<command<<endl;
			getline(cin,command);
		}
		cout<<endl;
		matrix::run("temp.m");
		file.close();
		remove("temp.m");
	}
	else{ // passing file as a parameter
		string fpath = argv[1];
		matrix::run(fpath);
	}
	return 0;
}
