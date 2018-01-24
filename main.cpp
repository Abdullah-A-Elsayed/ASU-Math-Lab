#include <iostream>
#include <string>
#include <fstream>
#include "matrix.h"
#include <map>
#include <algorithm>
#include <ctype.h>
#include <math.h>
#include <sstream>
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
	try{
		//whole = matrix::column_by_column(a,b);
		//whole.print_matrix();
		//string data ="C^3 * sin(1./D)";
		string data ="(1.2 + 3.4 - 5.6)/(2.1*3.2 + 4.6) - 12.1*3.1 + (1.2 + 5.2)^(4/(3.2+5.6))";
		//matrix::partial_Solve2("C^3 * sin(1./D)");
		matrix::Solve(data).print_matrix();
	}
	catch(string e){cout<<e<<endl;}
	/*--------------------------------------end of testing area------------------------------*/
}
else{ // passing file as a command argument
	string fpath = argv[1]; //file name
	//matrix::run(fpath); //phase1 function
	matrix::run_adv(fpath);
}
	//system ("pause");
	return 0;
}
