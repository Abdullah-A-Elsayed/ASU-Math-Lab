#include <iostream>
#include <string>
#include <fstream>
#include "matrix.h"
#include <map>
#include <algorithm>
using namespace std;
int main(int argc, char** argv)
{
	if(argc < 2){
		//matrix::run("example.m");
		//to be changed ...
		cout<<"Enter your commands ..."<<endl;
		int op_index;
		string command,name0,name1,name2;
		matrix input;
		map<const string, matrix> matrices;
		short print_flag = 1;
		try {
			while(1){
				if(print_flag&&command!=""){
					//cout<<name0<<endl;
					matrices[name0].print_matrix();cout<<endl;
				}
				getline(cin,command);	
				if(command[command.length()-1]==';') print_flag = 0 ;
				else print_flag = 1;
				if(command == "") continue;

				name0 = command.substr(0,command.find('=')-1);
				transform(name0.begin(),name0.end(),name0.begin(),::toupper);

				op_index = command.find('[');
				if(op_index != -1){
					matrix::handle_read(matrices,command,name0,op_index);
					continue;
				}

				op_index = command.find('+');
				if(op_index != -1){
					matrix::decode(command,name1,name2,op_index);
					matrices[name0] = matrices[name1].add_matrix(matrices[name2]);
					continue;
				}

				op_index = command.find('-');
				if(op_index != -1){
					matrix::decode(command,name1,name2,op_index);
					matrices[name0] = matrices[name1].sub_matrix(matrices[name2]);
					continue;
				}

				op_index = command.find('*');
				if(op_index != -1){
					matrix::decode(command,name1,name2,op_index);
					matrices[name0] = matrices[name1].mult_matrix(matrices[name2]);
					continue;
				}

				op_index = command.find("./");
				if(op_index != -1){
					matrix::decode(command,name1,name2,op_index+1);//+1 to get correct name2 and name1 is not important
					matrices[name0] = matrices[name2].inverse_matrix();
					continue;
				}

				op_index = command.find('/');
				if(op_index != -1){
					matrix::decode(command,name1,name2,op_index);
					matrices[name0] = matrices[name1].div_matrix(matrices[name2]);
					continue;
				}

				op_index = command.find("'");
				if(op_index != -1){
					command+="extra";
					matrix::decode(command,name1,name2,op_index+1);
					matrices[name0] = matrices[name1].transpose_matrix();
					continue;
				}
			}
		}
		catch(string e){ cout<<e<<endl;}
	}
	else{ // passing file as a parameter
		string fpath = argv[1];
		matrix::run(fpath);
	}
	return 0;
}
