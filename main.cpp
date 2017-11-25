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
				}
				getline(cin,command);	
				if(command[command.length()-1]==';') print_flag = 0 ;
				else print_flag = 1;
				if(command == "") continue;

				name0 = command.substr(0,command.find('=')-1);
				transform(name0.begin(),name0.end(),name0.begin(),::toupper);

				//op_index = command.find('[');
				//if(op_index != -1){
				//	matrix::handle_read(matrices,command,name0,op_index);
				//	continue;
				//}

				//op_index = command.find('+');
				//if(op_index != -1){
				//	matrix::matrix::decode(command,name1,name2,op_index);
				//	matrices[name0] = matrices[name1].add_matrix(matrices[name2]);
				//	continue;
				//}

				//op_index = command.find('-');
				//if(op_index != -1){
				//	matrix::matrix::decode(command,name1,name2,op_index);
				//	matrices[name0] = matrices[name1].sub_matrix(matrices[name2]);
				//	continue;
				//}

				//op_index = command.find('*');
				//if(op_index != -1){
				//	matrix::matrix::decode(command,name1,name2,op_index);
				//	matrices[name0] = matrices[name1].mult_matrix(matrices[name2]);
				//	continue;
				//}

				///*op_index = command.find("./");
				//if(op_index != -1){
				//	matrix::matrix::decode(command,name1,name2,op_index);//+1 to get correct name2 and name1 is not important
				//	matrices[name0] = matrices[name1].div_matrix(matrices[name2]);
				//	continue;
				//}*/op_index = command.find('/');
				//if(op_index != -1){
				//	matrix::matrix::decode(command,name1,name2,op_index);
				//	matrices[name0] = matrices[name1].div_matrix(matrices[name2]);

				//	continue;
				//}

				//op_index = command.find('./');
				//int equal_index = command.find_last_of('=');
				//	string b = command.substr(equal_index+2,op_index-equal_index-3);
				//	for (int i =0; i<b.size(); i++) {
    //        if  ((b[i] >= 'A' && b[i] <= 'Z') ||
    //               (b[i] >= 'a' && b[i] <= 'z')) {
				//if(op_index != -1){
				//	matrix::matrix::decode(command,name1,name2,op_index);
				//	matrices[name0] = matrices[name1].bitwisediv_matrix(matrices[name2]);
				//	continue;}
				//break;}
			 //else{  
				//	double c =stod(b);
				//	if(op_index != -1){
				//		matrix::matrix::decode(command,name1,name2,op_index+1);
				//		matrices[name0] = matrices[name2].bitwisediv2_matrix(c);
				//		continue;}
			
				//	   break;}}
		

				//op_index = command.find("'");
				//if(op_index != -1){
				//	command+="extra";
				//	matrix::matrix::decode(command,name1,name2,op_index+1);
				//	matrices[name0] = matrices[name1].transpose_matrix();
				//	continue;
				//}
				/*******************/
				op_index = command.find('[');
					if(op_index != -1){
						//cout<<command<<endl<<endl;
						
						matrix::handle_read(matrices,command,name0,op_index);
						if (command[command.length()-1]!=';')
						{
							cout<<name0<<": "<<endl;
							matrices[name0].print_matrix();cout<<endl;
						}
						continue;
					}

					op_index = command.find('+');
					if(op_index != -1){
						
						matrix::decode(command,name1,name2,op_index);
						cout<<name0<<": "<<endl;
						matrices[name0] = matrices[name1].add_matrix(matrices[name2]);
						matrices[name0].print_matrix();cout<<endl;
						continue;
					}

					op_index = command.find('-');
					if(op_index != -1){
						
						matrix::decode(command,name1,name2,op_index);
						cout<<name0<<": "<<endl;
						matrices[name0] = matrices[name1].sub_matrix(matrices[name2]);
						matrices[name0].print_matrix();cout<<endl;
						continue;
					}

					op_index = command.find('*');
					if(op_index != -1){
						
						matrix::decode(command,name1,name2,op_index);
						cout<<name0<<": "<<endl;
						matrices[name0] = matrices[name1].mult_matrix(matrices[name2]);
						matrices[name0].print_matrix();cout<<endl;
						continue;
					}

					op_index = command.find("'");
					if(op_index != -1){
						
						command+="extra";
						matrix::decode(command,name1,name2,op_index+1);
						cout<<name0<<": "<<endl;
						matrices[name0] = matrices[name1].transpose_matrix();
						matrices[name0].print_matrix();cout<<endl;
						continue;
					}

					//
					op_index = command.find_last_of('.');
					int bitWise =  command.find('/');
					if(op_index != -1 && bitWise != -1){
					//int equal_index = command.find_last_of('=');
					//string b = command.substr(equal_index+2,op_index-equal_index-3);
					matrix::decode(command,name1,name2,op_index); string b = name1;
					for (int i =0; i<b.length(); i++) {
						if((b[i] >= 'A' && b[i] <= 'Z') || (b[i] >= 'a' && b[i] <= 'z')) {
							cout<<name0<<": "<<endl;
							matrices[name0] = matrices[name1].bitwisediv_matrix(matrices[name2]);
							matrices[name0].print_matrix();cout<<endl;
							break;
						}
						else if(i==b.length()-1){
						    double c = atof(b.c_str());
							cout<<name0<<": "<<endl;
							matrices[name0] = matrices[name2].bitwisediv2_matrix(c);
							matrices[name0].print_matrix();cout<<endl;
					   }
					}
					continue ;
				}
					//

					op_index = command.find("/");
					if(op_index != -1){
						matrix::decode(command,name1,name2,op_index);
						cout<<name0<<": "<<endl;
						matrices[name0] = matrices[name1].div_matrix( matrices[name2]);
						matrices[name0].print_matrix();cout<<endl;
						continue;
						

					}
				/*********************/
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
