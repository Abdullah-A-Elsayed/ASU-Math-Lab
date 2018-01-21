#include <iostream>
#include <string>
#include <fstream>
#include "matrix.h"
#include <map>
#include <algorithm>
#include <ctype.h>
#include <math.h>
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
	map<const string, matrix> matrices;
	matrix a = matrix("1 5 6; 8 9 6; 2 5 1");
	matrix b = matrix("1 5 6; 8 9 6; 2 5 2");
	matrix c = matrix("1 5 6; 8 9 6; 2 5 3");
	matrices["abdullah"]=a;
	matrices["b"]=b;
	matrices["c"]=c;
	matrix whole;
	whole.fill_matrix_adv("abdullah;b",matrices);
	matrix powered = matrix::Rand(5,5);
	//whole.print_matrix();
	powered.print_matrix();
	/*--------------------------------------end of testing area------------------------------*/
}
else{ // passing file as a command argument
	string fpath = argv[1]; //file name
	//matrix::run(fpath); //phase1 function
	matrix::run_adv(fpath);
}
//	system ("pause");
	return 0;
}
