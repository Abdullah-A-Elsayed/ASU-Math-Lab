#include <iostream>
#include <string>
#include <fstream>
#include "matrix.h"
#include <map>
#include <algorithm>
#include <ctype.h>
using namespace std;
int main(int argc, char** argv)
{ 	
if(argc < 2){//no file passed to the program
		/***********************************taking commands***************************************/
		//cout<<"Enter your commands ..."<<endl; 
		//string command;
		//map<const string, matrix> matrices;
		//while(1){
		//	getline(cin,command);
		//	try{
		//		matrix::run_old_command(command,matrices);
		//	}
		//	catch(string e){ cout<<e<<endl;}
		//}
		/*************************************testing area***************************************/
	
}
else{ // passing file as a command argument
	string fpath = argv[1]; //file name
	//matrix::run(fpath); //phase1 function
	matrix::run_adv(fpath);
}
//	system ("pause");
	return 0;
}
