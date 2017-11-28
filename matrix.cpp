#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
#include <algorithm>
#include <cstdlib>
#include <math.h>
#include "matrix.h"
#include <iomanip>
using namespace std;
double my_abs(double& m ){
	return (m<0)? -m:m;
}
//private
	matrix matrix::inverse_2(){

		//*********************gauss jordan***************************
		int i, j, k, n;
		n = this->num_rows;
		//a should have dimensions 2n+1 * 2n+1
		vector< vector<double> > a(2*n+1);
		//pushing with zeros
		for(i=0;i<2*n+1;++i){
			for(j=0;j<2*n+1;++j){
				a[i].push_back(0);
			}
		}
		// taking from this->values from 1 to n
		for(i=0;i<n;++i){
			for(j=0;j<n;++j){
				a[i+1][j+1]=this->values[i][j];
			}
		}
		// modifing extension to eye
		double d; //?

		for (i = 1; i <= n; i++){

			for (j = 1; j <= 2 * n; j++){

				if (j == (i + n)){

					a[i][j] = 1;
				}
			}
		}

 

		/************** partial pivoting **************/

		for (i = n; i > 1; i--)

		{

			if (a[i - 1][1] < a[i][1])

				for (j = 1; j <= n * 2; j++)

				{

					d = a[i][j];

					a[i][j] = a[i - 1][j];

					a[i - 1][j] = d;

				}

		}

		/*cout << "pivoted output: " << endl;

		for (i = 1; i <= n; i++)

		{

			for (j = 1; j <= n * 2; j++)

				cout << a[i][j] << "    ";

			cout << endl;

		}*/

		/********** reducing to diagonal  matrix ***********/

 

		for (i = 1; i <= n; i++)

		{

			for (j = 1; j <= n * 2; j++)

				if (j != i)

				{
					//if diagonal has zeros?
					//me:******avoiding zeros in i i position
			if(a[i][i]==0){
				int good_row = i;//will be changed
				for(int s=i+1; s<=n ;++s){ //looking for good row
					if(a[s][i]!=0){
						good_row= s;
						break;
					}
				}
				if(good_row==i){string e = "non invertable matrix"; throw(e);};
				//swapping  what is pay for swapping ?? row in good pos is * -1
				double st;
				for(int s=1; s<=2*n;++s){
					st = a[i][s];
					a[i][s] = a[good_row][s];
					a[good_row][s] = -st;
				}
			}
			//***********************************
					d = a[j][i] / a[i][i];

					for (k = 1; k <= n * 2; k++)

						a[j][k] -= a[i][k] * d;

				}

		}

		/************** reducing to unit matrix *************/

		for (i = 1; i <= n; i++)

		{

			d = a[i][i];

			for (j = 1; j <= n * 2; j++)

				a[i][j] = a[i][j] / d;

		}

		matrix r; r.num_columns=n;r.num_rows=n;
	//	cout << "your solutions: " << endl;
		vector<double> row;
		for (i = 1; i <= n; i++)
		{
			
			for (j = n + 1; j <= n * 2; j++)
				row.push_back(a[i][j]);
				//cout << a[i][j] << "    ";
			r.values.push_back(row);
			row.clear();
		}

		return r;
		//************************************************
	}
	double matrix::determinant_2(int n){ //gauss elimination
		int i,j,k,swaps=1;
		//cout.precision(4);        //set precision
		//cout.setf(ios::fixed);
        //declare an array to store the elements of augmented-matrix  
		vector< vector<double> > a;
		a = this->values;
		double det=1; 
		int flag=0; 
		
		for (i=0;i<n;i++)                    //Pivotisation
			for (k=i+1;k<n;k++)
				if (my_abs(a[i][i])<my_abs(a[k][i])){
            		flag++;
            		for (j=0;j<n;j++){
						double temp=a[i][j];
						a[i][j]=a[k][j];
						a[k][j]=temp;
					}
				}
                
		//cout<<"\nThe matrix after Pivotisation is:\n";
		//for (i=0;i<n;i++)            //print the new matrix
		//{
		//	for (j=0;j<n;j++)
		//		cout<<a[i][j]<<setw(16);
		//	cout<<"\n";
		//}   
		for (i=0;i<n-1;i++){            //loop to perform the gauss elimination
			//me:******avoiding zeros in i i position
			if(a[i][i]==0){
				int good_row = i;//will be changed
				for(int s=i+1; s<n ;++s){ //looking for good row
					if(a[s][i]!=0){
						swaps*=-1;
						good_row= s;
						break;
					}
				}
				if(good_row==i) return 0;
				//swapping
				double st;
				for(int s=0; s<n;++s){
					st = a[i][s];
					a[i][s] = a[good_row][s];
					a[good_row][s] = st;
				}
			}
			//***********************************
			for (k=i+1;k<n;k++)
				{
					double t=a[k][i]/a[i][i];  //dividing here <<<
					for (j=0;j<n;j++)
						a[k][j]=a[k][j]-t*a[i][j];    //make the elements below the pivot elements equal to zero or elimnate the variables
				}
		}
    
		//cout<<"\n\nThe matrix after gauss-elimination is as follows:\n";
		//for (i=0;i<n;i++)            //print the new matrix
		//{
		//	for (j=0;j<n;j++)
		//		cout<<a[i][j]<<setw(16);
		//	cout<<"\n";
		//}
	
		for(i=0;i<n;i++){
			det=det*a[i][i];
		}            
		if (flag%2==0){
			det=det;
		}else{
			det=-det;
		}
		det*=swaps;
		//cout<<"\n\n The determinant is: "<<det;
		//if(det != det) return 0;
		return det;
	}


	int matrix::check_zero_dete()
	{
		int zfg=0;
		for(int i =0;i<num_rows;i++)
		{
					for(int j=i+1;j<num_rows;j++)
				{
					if(values[i][0]==values[j][0])
						{
							int counter=1;
							for (int k=1;k<num_columns;k++) 
							{if (values[i][k]==values[j][k]) counter++;else break;}
							if (counter==num_columns) { zfg=1; return zfg;}
						}
				}
		}

	return zfg;
	}
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
				fac.values[q][p] = pow(double(-1), p+q) * b.determinant_2((f - 1));
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
	void matrix::initialize_by1(int rows, int cols){ 
        this->num_rows = rows;
        this->num_columns = cols;
         //pushing values with zeros (initialization)
         for(int i=0 ; i<rows ; ++i){ //rows
            vector<double> row;
            for(int j=0 ; j<cols ; ++j){ //columns
                row.push_back(1);
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
					if(start<data.length()){
						if(data[start]==' ') start++;
					}
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
	matrix matrix:: bitwisediv2_matrix(double c){
		 matrix result;
		 
     result.initialize(this->num_rows,this->num_columns);
     for(int i=0;i<this->num_rows;i++){
      for(int j=0;j<this->num_columns;j++){
		  if(this->values[i][j]==0){
			  string e = "Warning: Divivding by Zero element, aborting \n";
			  throw(e);
		  }
     result.values[i][j] = c / this->values[i][j];
        }
     }

        return result;

	}

	matrix matrix:: bitwisediv_matrix(matrix &m){
	if(this->num_columns!=m.num_columns || this->num_rows != m.num_rows){
			string error = "can't div 2 matrices with different dimensions, Aborting ...";
			throw(error);
		}

     matrix result;
     result.initialize(this->num_rows,this->num_columns);
     for(int i=0;i<this->num_rows;i++){
      for(int j=0;j<this->num_columns;j++){
		  if(m.values[i][j]==0){
			  string e = "Warning: Divivding by Zero element, aborting.. \n";
			  throw(e);
		  }
     result.values[i][j] = this->values[i][j]/m.values[i][j];
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
		string error;
        int a =this-> num_rows;
		int b = m.get_num_columns();
		int c = m.get_num_rows();
        int d =this-> num_columns;
		if(d!=c){ error="can't multiply 2 matrices while 1st cols not equal to 2nd rows"; throw(error);}
        
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
		string error;
		double det_val;
		if (this->num_rows != this->num_columns){
			error = "No inverse for non-square matrix, calculating inverse is aborted"; throw(error);
		}
		det_val = this->determinant_2(this->num_rows);
		if(det_val == 0){
			error = "warning: Dividing by zero, aborting... \n" ; throw(error);
		}
				
		// strat to get the inverse for the matrix
		else
		{
			//matrix m; m.initialize(this->num_rows, this->num_columns);
			//m = this->cal_cofactor(this->num_rows); //very very slow
			//using inverse2
			return this->inverse_2();
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
	matrix matrix ::Rand(int a,int b){
	 matrix result;
		result.initialize(a,b);
		for(int i=0; i< a; i++){
			for(int j=0; j<b; j++){
				result.values[i][j]+=rand();
				
			}
		}

		return result;
	    }
	matrix matrix:: Eye (int a,int b){
	 matrix result;
		result.initialize(a,b);
		for(int i=0; i< a; i++){
			for(int j=0; j<b; j++){
				if(i==j){
				result.values[i][j]=1;
				}
				else {
                result.values[i][j]=0;
		}
		}
		}
		return result;
	
	}
	void matrix::handle_read(map<const string, matrix>& matrices,string command,string name0,int op_index){
		int n_deleted = 1;
		if(command[command.length()-1]==';') n_deleted++;
		string values = command.substr(op_index+1,command.length()-op_index-1-n_deleted);
		matrix x(values);
		matrices[name0] = x;
	}
	/*void matrix::decode(string command,string& name1,string& name2,int op_index){
		int equal_index = command.find_last_of('=');
		name1 = command.substr(equal_index+2,op_index-equal_index-3);
		transform(name1.begin(),name1.end(),name1.begin(),::toupper);

		name2 = command.substr(op_index+2,command.length()-op_index-2);
		transform(name2.begin(),name2.end(),name2.begin(),::toupper);
	}*/
		void matrix::decode(string command,string& name1,string& name2,int op_index){
		int equal_index = command.find_last_of('=');
		name1 = command.substr(equal_index+2,op_index-equal_index-3);
		transform(name1.begin(),name1.end(),name1.begin(),::toupper);

		int name2_begin = op_index + 2;
		if(command[op_index]=='.') name2_begin++;
		name2 = command.substr(name2_begin,command.length()-name2_begin);
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
			int op_index; //holds position of the operation
			while(getline(file,command)){
				if(command == "" || command[0]=='#'|| (command[0]=='/'&&command[1]=='/')) continue;

				name0 = command.substr(0,command.find('=')-1);
				transform(name0.begin(),name0.end(),name0.begin(),::toupper);
				try{
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
						if (command[command.length()-1]!=';')
						{
							cout<<name0<<": "<<endl;
							matrices[name0].print_matrix();cout<<endl;
						}
						continue;
					}

					op_index = command.find('+');
					if(op_index != -1){
						remove_back_slashes(command);
						decode(command,name1,name2,op_index);
						cout<<name0<<": "<<endl;
						matrices[name0] = matrices[name1].add_matrix(matrices[name2]);
						matrices[name0].print_matrix();cout<<endl;
						continue;
					}

					op_index = command.find('-');
					if(op_index != -1){
						remove_back_slashes(command);
						decode(command,name1,name2,op_index);
						cout<<name0<<": "<<endl;
						matrices[name0] = matrices[name1].sub_matrix(matrices[name2]);
						matrices[name0].print_matrix();cout<<endl;
						continue;
					}

					op_index = command.find('*');
					if(op_index != -1){
						remove_back_slashes(command);
						decode(command,name1,name2,op_index);
						cout<<name0<<": "<<endl;
						matrices[name0] = matrices[name1].mult_matrix(matrices[name2]);
						matrices[name0].print_matrix();cout<<endl;
						continue;
					}

					op_index = command.find("'");
					if(op_index != -1){
						remove_back_slashes(command);
						command+="extra";
						decode(command,name1,name2,op_index+1);
						cout<<name0<<": "<<endl;
						matrices[name0] = matrices[name1].transpose_matrix();
						matrices[name0].print_matrix();cout<<endl;
						continue;
					}

					//
					op_index = command.find_last_of('.');
					int bitWise =  command.find('/');
					if(op_index != -1 && bitWise != -1 && bitWise==op_index+1){
					//int equal_index = command.find_last_of('=');
					//string b = command.substr(equal_index+2,op_index-equal_index-3);
					remove_back_slashes(command);
					decode(command,name1,name2,op_index); string b = name1;
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
						remove_back_slashes(command);
						decode(command,name1,name2,op_index);
						cout<<name0<<": "<<endl;
						matrices[name0] = matrices[name1].div_matrix( matrices[name2]);
						matrices[name0].print_matrix();cout<<endl;
						continue;
						

					}
			}
			catch(string e){ cout<<e<<endl;}
		}
			file.close();
			// for(map<const string, matrix>::iterator i = matrices.begin(); i!=matrices.end();++i){
			// 	cout<<i->second.get_num_rows()<<"*"<<i->second.get_num_columns();
			// 	cout<<" "<<i->first<<":\n";i->second.print_matrix();cout<<endl;
			// }
		}
		else{
			cout<<"error opening file"<<endl;
		}
	}
	matrix matrix ::Sin(){
		matrix result ;
		result.initialize(this->num_rows,this->num_columns);
		for(int i=0;i<this->num_rows;i++){
			for(int j=0;j<this->num_columns;j++){
				result.values[i][j]= sin(this->values[i][j]);
			
			}
		}

	return result;
	}

	matrix matrix:: Sqrt(){

		matrix result ;
		result.initialize(this->num_rows,this->num_columns);
		for(int i=0;i<this->num_rows;i++){
			for(int j=0;j<this->num_columns;j++){
				result.values[i][j]=sqrt(this->values[i][j]);
			
			}
		}
		return result;

	}

	