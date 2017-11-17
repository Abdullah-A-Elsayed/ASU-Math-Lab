#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
#include <algorithm>
#include <cstdlib>
#include <math.h>
#include "matrix.h"
using namespace std;
//private
	double matrix::cal_determin_sq(int num_rows)
	{
		double s = 1, det = 0;
		matrix b; b.initialize(this->num_rows, this->num_columns);
		int i, j, m, n, c, k = num_rows;
		if (k == 1)
		{
			return (this->values[0][0]);
		}
		else
		{
			det = 0;
			for (c = 0; c<k; c++)
			{
				m = 0;
				n = 0;
				for (i = 0; i<k; i++)
				{
					for (j = 0; j<k; j++)
					{
						b.values[i][j] = 0;
						if (i != 0 && j != c)
						{
							b.values[m][n] = this->values[i][j];
							if (n<(k - 2))
								n++;
							else
							{
								n = 0;
								m++;
							}
						}
					}
				}
				det = det + s * (this->values[0][c] * b.cal_determin_sq((k - 1)));
				s = -1 * s;
			}
		}

		return (det);
	}

	// calculate co-factor func
	matrix matrix::cal_cofactor(int num_rows)
	{
		matrix b; b.initialize(this->num_rows, this->num_columns);
		matrix fac; fac.initialize(this->num_rows, this->num_columns);
		matrix inverse; inverse.initialize(this->num_rows, this->num_columns);
		int f = num_rows;
		int p, q, m, n, i, j;
		for (q = 0; q < f; q++)
		{
			for (p = 0; p < f; p++)
			{
				m = 0;
				n = 0;
				for (i = 0; i < f; i++)
				{
					for (j = 0; j < f; j++)
					{
						if (i != q && j != p)
						{
							b.values[m][n] = this->values[i][j];
							if (n < (f - 2))
								n++;
							else
							{
								n = 0;
								m++;
							}
						}
					}
				}
				fac.values[q][p] = pow(double(-1), p+q) * b.cal_determin_sq((f - 1));
			}
		}

		return (this->transpose(fac, f));
	}

	// find transpose func
	matrix matrix::transpose(matrix& fac, double r)
	{
		int i, j;
		matrix b; b.initialize(this->num_rows, this->num_columns);
		matrix inverse; inverse.initialize(this->num_rows, this->num_columns);
		double d;

		for (i = 0; i<r; i++)
		{
			for (j = 0; j<r; j++)
			{
				b.values[i][j] = fac.values[j][i];
			}
		}
		d = this->cal_determin_sq(r);
		for (i = 0; i<r; i++)
		{
			for (j = 0; j<r; j++)
			{
				inverse.values[i][j] = b.values[i][j] / d;
			}
		}
		return inverse;
	}
    /*-----------------------Gasser end of assisting private functions to calculate matrix inverse --------------*/
    //public:
    matrix::matrix(){ // constructing ...
        this->num_rows = this->num_columns =0;
	}
	
	matrix::matrix(string nums){ // constructing ...
        this->fill_matrix(nums);
    }

    void matrix::initialize(int rows, int cols){ // taking dimensions
        this->num_rows = rows;
        this->num_columns = cols;
         //pushing values with zeros (initialization)
         for(int i=0 ; i<rows ; ++i){ //rows
            vector<double> row;
            for(int j=0 ; j<cols ; ++j){ //columns
                row.push_back(0);
            }
            values.push_back(row);
        }
    }

    void matrix::print_matrix(){ // print matrix for testing
        for(int i=0 ; i< this-> num_rows ; ++i){
            for(int j=0 ; j< this->num_columns ; ++j){
                printf("%g \t",this->values[i][j]);
            }
            cout<<endl;
        }
    }


    // tasks:
    void matrix::fill_matrix (string data){
        // Aly
        // data will be like this "1.1 2 3.5; 9.6 5.2 4.7"
        // these are 2 rows and three columns ('; ' separates rows .. ' ' separates colums)
		// initialize using initialize function provided above then assign values
		if(num_rows){ // resetting
			values.clear();
			num_columns =0;
			num_rows = 0;
		}
		int start=0;
		int end;
		if(data[data.length()-1]!=';') {data = data +";";}

		vector<double> row;
		for(unsigned int i = 0 ; i< data.length(); i++){
			if ((data[i]==' '&&data[i-1]!=';')||(data[i]==';')){
				end=i;
				//ading previous num
				row.push_back(atof(data.substr(start, end-start).c_str()));
				start=i+1;
				if(data[i]==';'){
					this->values.push_back(row);
					row.clear();
					if(data[start]==' ') start++;
				}
			}
		}
		this->num_rows = this->values.size();
		this->num_columns = this->values[0].size();
	}

    matrix matrix::add_matrix( matrix& m){
        //Aya
        // create a result matrix with correct dimensions then initialize it using initialize function provided above
        // result = this + m
        // return result

		//error handling:
		if(this->num_columns!=m.num_columns || this->num_rows != m.num_rows){
			string error = "can't sum 2 matrices with different dimensions, Aborting ...";
			throw(error);
		}

         matrix result;
     result.initialize(this->num_rows,this->num_columns);
     for(int i=0;i<this->num_rows;i++){
      for(int j=0;j<this->num_columns;j++){
     result.values[i][j] = this->values[i][j]+m.values[i][j];
        }
     }

        return result;
    }

    matrix matrix::sub_matrix( matrix& m){
        //Do'aa
        // create a result matrix with correct dimensions then initialize it using initialize function provided above
        // result = this - m
        // return result

		//error handling:
		if(this->num_columns!=m.num_columns || this->num_rows != m.num_rows){
			string error = "can't subtract 2 matrices with different dimensions, Aborting ...";
			throw(error);
		}
              matrix r;
     r.initialize(this->num_rows,this->num_columns);
     for(int i=0;i<this->num_rows;i++){
      for(int j=0;j<this->num_columns;j++){
	  r.values[i][j] = this->values[i][j]- m.values[i][j];
        }
     }

        return r;
    }

    matrix matrix::mult_matrix( matrix& m){
        //Amira
        // create a result matrix with correct dimensions then initialize it using initialize function provided above
        // result = this * m
        // return result
        int a =this-> num_rows;
		int b = m.get_num_columns();
		int c = m.get_num_rows();
        int d =this-> num_columns;
        if(d!=c) throw("can't multiply 2 matrices while 1st cols not equal to 2nd rows");
        
        matrix result;
		result.initialize(a,b);
		for(int i=0; i< a; i++){
			for(int j=0; j<b; j++){
				for (int k =0; k<c; k++){
					result.values[i][j]+=this->values[i][k]*m.values[k][j];
				}
			}
		}

		return result;

    }

    // find inverse matrix
	matrix matrix::inverse_matrix()
	{
		double det_val = this->cal_determin_sq(this->num_rows); // calculate determine value for the matrix

		// intialize the inverse matrix with zeros & check if it is a square matrix
		matrix m; m.initialize(this->num_rows, this->num_columns);
		string error;
		if (this->num_rows != this->num_columns){ error = "No inverse for non-square matrix, calculating inverse is aborted"; throw(error); }

		// check if determine equals zero
		else if (det_val == 0){error = "No inverse for this zero-determine matrix, calculating inverse is aborted" ; throw(error); }

		// strat to get the inverse for the matrix
		else
		{
			m = this->cal_cofactor(this->num_rows);
			return (m);
		}

	}
	
	matrix matrix::transpose_matrix(){
		matrix r;
		r.initialize(this->num_columns,this->num_rows);
		for(int i = 0; i<num_rows ; ++i){
			for(int j = 0; j<num_columns ; ++j){
				r.values[j][i] = this->values[i][j];
			}
		}
		return r;
	}

    matrix matrix::div_matrix( matrix& m){
        //Alaa Ayman
        // create a result matrix with correct dimensions then initialize it using initialize function provided above
        // result = this * m inversed
        // use previous functions
        int a = this ->num_rows;
		int b = m.get_num_columns();

		matrix result;
		result.initialize(a,b);
		matrix x = m.inverse_matrix();
		result = this->mult_matrix(x);

		return result;



    }
	void matrix::handle_read(map<const string, matrix>& matrices,string command,string name0,int op_index){
		int n_deleted = 1;
		if(command[command.length()-1]==';') n_deleted++;
		string values = command.substr(op_index+1,command.length()-op_index-1-n_deleted);
		matrix x(values);
		matrices[name0] = x;
	}
	void matrix::decode(string command,string& name1,string& name2,int op_index){
		int equal_index = command.find_last_of('=');
		name1 = command.substr(equal_index+2,op_index-equal_index-3);
		transform(name1.begin(),name1.end(),name1.begin(),::toupper);

		name2 = command.substr(op_index+2,command.length()-op_index-2);
		transform(name2.begin(),name2.end(),name2.begin(),::toupper);
	}

	void matrix::remove_back_slashes(string& s){
		string u="";
		for(int i=0; i<s.length();++i){
			if(s[i]=='\r') continue;
			u+=s[i];
		}
		s = u;
	}

	void matrix::run(string fpath){
		ifstream file (fpath.c_str());
		if(file){//opened safely
			map<const string, matrix> matrices;
			string command, name0, name1, name2, sub_command="",line;
			int op_index;
			while(getline(file,command)){
				if(command == "") continue;

				name0 = command.substr(0,command.find('=')-1);
				transform(name0.begin(),name0.end(),name0.begin(),::toupper);

				op_index = command.find('[');
				if(op_index != -1){
					if(command.find("]")==-1){//not exist
						while(getline(file,line)){
							remove_back_slashes(command);
							command+=line;
							if(line.find(']')!=-1)break;
						}
					}
					//cout<<command<<endl<<endl;
					remove_back_slashes(command);
					handle_read(matrices,command,name0,op_index);
					continue;
				}

				op_index = command.find('+');
				if(op_index != -1){
					remove_back_slashes(command);
					decode(command,name1,name2,op_index);
					matrices[name0] = matrices[name1].add_matrix(matrices[name2]);
					continue;
				}

				op_index = command.find('-');
				if(op_index != -1){
					remove_back_slashes(command);
					decode(command,name1,name2,op_index);
					matrices[name0] = matrices[name1].sub_matrix(matrices[name2]);
					continue;
				}

				op_index = command.find('*');
				if(op_index != -1){
					remove_back_slashes(command);
					decode(command,name1,name2,op_index);
					matrices[name0] = matrices[name1].mult_matrix(matrices[name2]);
					continue;
				}

				op_index = command.find("./");
				if(op_index != -1){
					remove_back_slashes(command);
					decode(command,name1,name2,op_index+1);//+1 to get correct name2 and name1 is not important
					matrices[name0] = matrices[name2].inverse_matrix();
					continue;
				}

				op_index = command.find('/');
				if(op_index != -1){
					remove_back_slashes(command);
					decode(command,name1,name2,op_index);
					matrices[name0] = matrices[name1].div_matrix(matrices[name2]);
					continue;
				}

				op_index = command.find("'");
				if(op_index != -1){
					remove_back_slashes(command);
					command+="extra";
					decode(command,name1,name2,op_index+1);
					matrices[name0] = matrices[name1].transpose_matrix();
					continue;
				}
			}
			file.close();
			for(map<const string, matrix>::iterator i = matrices.begin(); i!=matrices.end();++i){
				cout<<i->second.get_num_rows()<<"*"<<i->second.get_num_columns();
				cout<<" "<<i->first<<":\n";i->second.print_matrix();cout<<endl;
			}
		}
		else{
			cout<<"error opening file"<<endl;
		}
	}