#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
#include <algorithm>
#include <cstdlib>
#include <math.h>
#include "matrix.h"
#include <iomanip>
#include<vector>
#include <wctype.h>
#define pi 3.1416
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
	int matrix::is_identify(double n){
		int m=0;
		for(int i=0;i<10;i++){
			if(n==(2*i+1)*pi/2){
				return m=1;
				break;
			}
			else return m=0;
		}
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

	matrix matrix:: ones(int n,int m){
		matrix result;
		result.initialize(n,m);
		for (int i=0; i<n;i++){
			for(int j=0; j<m;j++){
				result.values[i][j]= 1;
			}
		}
		return result;
	}

	matrix matrix::zeros(int r, int c)
	{
		matrix x;
		x.initialize(r, c);
		return x;
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

		/*matrix matrix ::Sin(){
		matrix result ;
		result.initialize(this->num_rows,this->num_columns);
		for(int i=0;i<this->num_rows;i++){
			for(int j=0;j<this->num_columns;j++){
				result.values[i][j]= sin(this->values[i][j]);
			
			}
		}

	return result;
	}*/

		matrix matrix ::Cos(){
		matrix result ;
		result.initialize(this->num_rows,this->num_columns);
		for(int i=0;i<this->num_rows;i++){
			for(int j=0;j<this->num_columns;j++){
				result.values[i][j]= cos( this->values[i][j]);
			
			}
		}

		return result;
	}
		
		matrix matrix ::Tan(){
		matrix result ;
		result.initialize(this->num_rows,this->num_columns);
		for(int i=0;i<this->num_rows;i++){
			for(int j=0;j<this->num_columns;j++){
	
			if(is_identify(this->values[i][j])==1){
			string error="math error";
			 throw(error);
			}
			else{result.values[i][j]= tan(this->values[i][j]);}
			
			}
		}

	return result;
	}




		matrix matrix :: element_wise_power(double a){
		
		matrix result;
		result.initialize(this->num_rows,this->num_columns);

		for(int i=0; i< this->num_rows; i++){
			for(int j=0; j< this->num_columns; j++){
			
				result.values[i][j]= pow(this->values[i][j], a);

			}
		}
//		result.print_matrix();
		return result;
		
		}


	/*matrix matrix:: Sqrt(){

		matrix result ;
		result.initialize(this->num_rows,this->num_columns);
		for(int i=0;i<this->num_rows;i++){
			for(int j=0;j<this->num_columns;j++){
				result.values[i][j]=sqrt(this->values[i][j]);
			
			}
		}
		return result;

	}*/
	matrix matrix:: Pow(int n){
		matrix result;
		result=matrix::Eye(this->num_rows,this->num_columns); //as EYE is el mohayed el darbee for matrix     ex : 1 0   *   2  3
		                                                                                                     //    0 1       4  6
		for(int i=0;i<n;i++){                                                                                 //ans=   2     3
			result = this->mult_matrix(result);                                                              //        4      6
		}
		return result;
	}
	
	//////////

	string matrix :: getString(){
		string result;
		result.clear();

		char substring[100];
		for(int i=0; i< this->num_rows; i++){
			for(int j=0; j<this->num_columns; j++){
				
				sprintf_s(substring, "%f", this->values[i][j]);
				result+=substring;

				if(j+1<this->num_columns) result+=" ";
				else if(i+1<this->num_rows) result+= "; ";
				else continue;
			}
		}
	//	cout<<result<<endl;
		return result;
	}

	
	matrix matrix:: add_const(double a){
		matrix result;
		result.initialize(this->num_rows, this->num_columns);

		for(int i=0; i<this->num_rows; i++){
			for(int j=0; j<num_columns;j++){
				result.values[i][j]=a+this->values[i][j];
			}
		}
	//	result.print_matrix();
		return result;
	}



	matrix matrix::column_by_column(matrix &a , matrix &b){
		matrix r;
		int n;
         int c;
	 if(a.get_num_rows()!=b.get_num_rows()){
		 string e="mismach number of rows";
		throw(e);

		}
			
		else{c=(a.get_num_columns() + b.get_num_columns());  

			r.initialize(a.get_num_rows(),c);
	
			for(int i=0;i<a.get_num_rows();i++){
				 n=a.get_num_columns();
				for(int j=0;j<c;j++){  
					if(j<n)
                  r.values[i][j]=a.values[i][j];

				else{
					r.values[i][j]=b.values[i][j-n];
				
				}
					
				
				}

			}
		}
		
			return r;
}
	//////////
		matrix matrix::row_by_row(matrix &a , matrix &b){
		matrix r;
		int n;
         int c;
	 if(a.get_num_columns()!=b.get_num_columns()){
		 string e="mismach number of columns";
		throw(e);

		}
			
		else{c=(a.get_num_rows() + b.get_num_rows());  

			r.initialize(c,a.get_num_columns());
	              n=a.get_num_rows();
			for(int i=0;i<c;i++){
		       
				for(int j=0;j<a.get_num_columns();j++){  
					if(i<n)
                  r.values[i][j]=a.values[i][j];

				else{
					r.values[i][j]=b.values[i-n][j];
						}
					
				
				}

			}
		}
		
			return r;
}
	//////////

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
	void matrix::run_old_command(string command, map<const string, matrix>& matrices){//processes given phase1 command(save&print)
		if(command == "" || command[0]=='#'|| (command[0]=='/'&&command[1]=='/')) return;
		int op_index;
		string name0,name1,name2;
		matrix input;
		short print_flag = 1;
		if(command[command.length()-1]==';') print_flag = 0 ;
		else print_flag = 1;
		name0 = command.substr(0,command.find('=')-1);
		transform(name0.begin(),name0.end(),name0.begin(),::toupper);

		op_index = command.find('[');
		if(op_index != -1){
			matrix::handle_read(matrices,command,name0,op_index);
			if (command[command.length()-1]!=';')
			{
				cout<<name0<<"= "<<endl;
				matrices[name0].print_matrix();cout<<endl;
			}
			return;
		}

		op_index = command.find('+');
		if(op_index != -1){				
			matrix::decode(command,name1,name2,op_index);
			cout<<name0<<"= "<<endl;
			matrices[name0] = matrices[name1].add_matrix(matrices[name2]);
			matrices[name0].print_matrix();cout<<endl;
			return;
		}

		op_index = command.find('-');
		if(op_index != -1){	
			matrix::decode(command,name1,name2,op_index);
			cout<<name0<<"= "<<endl;
			matrices[name0] = matrices[name1].sub_matrix(matrices[name2]);
			matrices[name0].print_matrix();cout<<endl;
			return;
		}

		op_index = command.find('*');
		if(op_index != -1){			
			matrix::decode(command,name1,name2,op_index);
			cout<<name0<<"= "<<endl;
			matrices[name0] = matrices[name1].mult_matrix(matrices[name2]);
			matrices[name0].print_matrix();cout<<endl;
			return;
		}

		op_index = command.find("'");
		if(op_index != -1){		
			command+="extra";
			matrix::decode(command,name1,name2,op_index+1);
			cout<<name0<<"= "<<endl;
			matrices[name0] = matrices[name1].transpose_matrix();
			matrices[name0].print_matrix();cout<<endl;
			return;
		}
		//
		op_index = command.find_last_of('.');
		int bitWise =  command.find('/');
		if(op_index != -1 && bitWise != -1 && bitWise==op_index+1){
		//int equal_index = command.find_last_of('=');
		//string b = command.substr(equal_index+2,op_index-equal_index-3);
		matrix::remove_back_slashes(command);
		matrix::decode(command,name1,name2,op_index); string b = name1;
		for (int i =0; i<b.length(); i++) {
			if((b[i] >= 'A' && b[i] <= 'Z') || (b[i] >= 'a' && b[i] <= 'z')) {
				cout<<name0<<"= "<<endl;
				matrices[name0] = matrices[name1].bitwisediv_matrix(matrices[name2]);
				matrices[name0].print_matrix();cout<<endl;
				break;
			}
			else if(i==b.length()-1){
				double c = atof(b.c_str());
				cout<<name0<<"= "<<endl;
				matrices[name0] = matrices[name2].bitwisediv2_matrix(c);
				matrices[name0].print_matrix();cout<<endl;
			}
		}
		return;
	}
		//
		op_index = command.find("/");
		if(op_index != -1){
			matrix::remove_back_slashes(command);
			matrix::decode(command,name1,name2,op_index);
			cout<<name0<<"= "<<endl;
			matrices[name0] = matrices[name1].div_matrix( matrices[name2]);
			matrices[name0].print_matrix();cout<<endl;
			return;
		}
	} 
	void matrix::run(string fpath)
	{
		ifstream file (fpath.c_str());
		if(file){//opened safely
			map<const string, matrix> matrices;
			string command,sub_command="",line;
			int op_index;
			while(getline(file,command)){ //read line by line
				//make sure that command is complete (may be multi lines)
				op_index = command.find('[');
				if(op_index != -1){
					if(command.find("]")==-1){//not exist
						while(getline(file,line)){
							remove_back_slashes(command);
							command+=line;
							if(line.find(']')!=-1)break;
						}
					}
					remove_back_slashes(command);
				}
				//process the command
				try{
					matrix::run_old_command(command, matrices);
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
/* -----------------------------------------Advanced File example------------------------------------------*/

	/* run advanced */

	void matrix::run_adv(string fpath)
	{
		ifstream file(fpath.c_str());
		if (file){//opened safely
		map<const string, matrix> matrices;
		string command, name0, name1, name2, sub_command = "", line;
		int op_index; //holds position of the operation
		while (getline(file, command))
		{
			remove_spaces(command); /* makes the line doesn't start with a space*/
			int prnt_fg = 1; /*this is to rmove the semicolon at the end cuz it breaks if it has*/
			if (command[command.length() - 1] == ';') { command = command.substr(0, command.length() - 1); prnt_fg = 0; }
			if (command == "" || command[0] == '#' || (command[0] == '/'&&command[1] == '/')) continue;

			/* if the command didn't have a name it will be named ans */
			if ((command[0] >= 'A' && command[0] <= 'Z') || (command[0] >= 'a' && command[0] <= 'z'))
			{
				name0 = command.substr(0, command.find('=') - 1); /*this means the name must have a space after it*/
				transform(name0.begin(), name0.end(), name0.begin(), ::toupper);
			}
			else
			{
				name0 = "ans";
			}
			/* end if the command didn't have a name it will be named ans */

			try{
				
				/* detect joined matrix  gasser */

				int chk_mat = command.find_first_of('[');
				if (chk_mat != -1){
					op_index = chk_mat;
						int opn_brac_count = 0, cls_brac_count = 0, end_comd_fg = 0;
						for (int i = 0; i < command.length(); i++)
						{
							if (command[i] == '[') opn_brac_count++;
							if (command[i] == ']') cls_brac_count++;
						}
						if (opn_brac_count == cls_brac_count){ end_comd_fg = 1; }//cout << command << endl; }
						else {
							while (end_comd_fg == 0)
							{
								getline(file, line);
								remove_spaces(line); /* makes the line doesn't start with a space*/
								remove_back_slashes(command);
								command += ' ';		/* adds a space for detiction purposes*/
								command += line;
								prnt_fg = 1;/*this is to rmove the semicolon at the end cuz it breaks if it has*/
								if (command[command.length() - 1] == ';') { command = command.substr(0, command.length() - 1); prnt_fg = 0; }
								opn_brac_count = 0; cls_brac_count = 0;
								for (int i = 0; i < command.length(); i++)
								{
									if (command[i] == '[') opn_brac_count++;
									if (command[i] == ']') cls_brac_count++;
								}
								if (opn_brac_count == cls_brac_count){ end_comd_fg = 1; } //cout << command << endl;}
							}
						}
						remove_back_slashes(command);
						handle_read_adv(matrices, command, name0, op_index);
						if (prnt_fg)
						{
						cout << name0 << ": " << endl;
						matrices[name0].print_matrix();
						}
						continue;
				}

				/* end detect joined matrix*/
				

		/*------------------------------------------------- Adv file Tasks -----------------------------------------------------*/

				/* detect lines [ rand / eye / zeros / ones ] */

				/* End detect lines [ rand / eye / zeros / ones ] */




				/* detect lines x / y & showing matrix from just name */


				/* End detect lines x / y & showing matrix from just name */


		/*------------------------------------------------------- Adv file Tasks ----------------------------------------------------*/


				/* detect solve  A / L */
				else
				{
					if (command.find('=')||command.find('+')||command.find('-')||command.find('/')||command.find('*')||command.find('^')
				 || command.find('s') || command.find('c') || command.find('t') || command.find('('))
					{
						/*there is a problem with the solve function "vector out of range" 
						I tested it with this test case "A = 5.5 + 12 * sin(0.4) + 2.2^4;" */

						matrices[name0] = solve_elemnt(command);
						if (prnt_fg)
						{
							cout << name0 << ": " << endl;
							matrices[name0].print_matrix();
						}
						continue;
					}
					else { continue; }
				}
				/* end detect solve  A / L */

				}

			catch (string e){ cout << e << endl; }
		}

		file.close();

		}
		else{

			cout << "error opening file" << endl;
		}
	}

	/* end run advanced */

	/*handle read adv-gasser*/
	void matrix::handle_read_adv(map<const string, matrix>& matrices, string command, string name0, int opn_index)
	{
		int semi_index; matrix jn_mat[20]; int jn_mat_fg = 0, jn_mat_ord = 0; int c_fg = 0;
		for (int i = 0; i < command.length(); i++)
		{
			/*c - line*/

			if  (command[i] == '['&& command[i + 1] == '[' && 
				((command[i + 2] >= 'A' && command[i + 2] <= 'Z') || (command[i + 2] >= 'a' && command[i + 2] <= 'z')))
			{
				c_fg = 1; int nxt_opn_brac;
				for (int f = i + 2; f < command.length(); f++)
				{
					if (command[f] == '[') nxt_opn_brac = f;
				}
				string mat_ltr = command.substr(i + 2, nxt_opn_brac - i - 3);
				jn_mat[jn_mat_ord] = matrices[mat_ltr];
				jn_mat_ord++;
				continue;
			}

			if (command[i] == ']'&& command[i + 1] == ' ' && c_fg == 1)
			{
				int nxt_cls_brac; //int frst_num;
				for (int f = i + 1; f < command.length(); f++)
				{
					if (command[f] == ']') nxt_cls_brac = f;
				}
				string mat_val = command.substr(i + 2, nxt_cls_brac - i - 2);
				remove_space_after_semis(mat_val);
				cut_mat_solve(mat_val);
				/*for (int f = 0; f < new_comnd.length(); f++)
				{
					if (new_comnd[f] != ' ') { frst_num = f; break; }
				}

				nxt_cls_brac = new_comnd.find(']');
				string final_comnd = new_comnd.substr(frst_num);*/

				matrix z;
				z.fill_matrix_adv(mat_val, matrices);
				jn_mat[jn_mat_ord] = z;
				jn_mat_ord++;
				continue;
			}

			/*End c - line*/

			if (command[i] == ';')
			{ 
				semi_index = i;
				if (command[i + 1] == '[') 
				{
					jn_mat_fg = 1;
					string mat_vals = command.substr(opn_index + 1, semi_index - opn_index);
					remove_space_after_semis(mat_vals);
					cut_mat_solve(mat_vals);
					matrix x;
					x.fill_matrix_adv(mat_vals, matrices);
					jn_mat[jn_mat_ord] =x;
					jn_mat_ord++;
				}
				else if (command[i - 1] == ']')
				{
					jn_mat_fg = 1; int cls_brac_indx;
					for (int i = semi_index + 1; i < command.length(); i++) if (command[i] == ']'){ cls_brac_indx = i; break; }
					string mat_vals = command.substr(semi_index + 1, cls_brac_indx - semi_index - 1);
					remove_space_after_semis(mat_vals);
					cut_mat_solve(mat_vals);
					matrix x;
					x.fill_matrix_adv(mat_vals, matrices);
					jn_mat[jn_mat_ord] = x;
					jn_mat_ord++;
				}
			}

			if (command[i] == '[') opn_index = i;

			if (command[i] == ']' && command[i - 1] != ']')
			{
				string mat_vals = command.substr(opn_index + 1, i - 1 - opn_index);
				remove_space_after_semis(mat_vals);
				cut_mat_solve(mat_vals);
				matrix y;
				y.fill_matrix_adv(mat_vals, matrices);
				jn_mat[jn_mat_ord] = y;
				jn_mat_ord++;
				if (jn_mat_fg == 0 ) matrices[name0] = y;
			}
		}

		//combine joined matrix
		if (jn_mat_fg||c_fg)
		{
			int mx_colms_ord = 0; int mn_colms_ord = 0;
			for (int i = 0; i < jn_mat_ord-1; i++)
			{
				//cout << "this: "<< endl; jn_mat[i].print_matrix();

				if (jn_mat[i].num_columns != jn_mat[i + 1].num_columns && jn_mat[i].num_rows == jn_mat[i + 1].num_rows) //colm by colm
				{
					jn_mat[i] = column_by_column(jn_mat[i], jn_mat[i + 1]);
					jn_mat[i + 1] = jn_mat[i + 2];/*look how to delete a row from array*/
					i = 0; jn_mat_ord--;
				}
				if (jn_mat[i].num_columns == jn_mat[i + 1].num_columns && jn_mat[i].num_rows != jn_mat[i + 1].num_rows) //row by row
				{
					jn_mat[i] = row_by_row(jn_mat[i], jn_mat[i + 1]);
					jn_mat[i + 1] = jn_mat[i + 2]; /*look how to delete a row from array*/
					i = 0; jn_mat_ord--;
				}

			}

			matrices[name0] = jn_mat[0];
		}
	}

	/* end handle read adv-gasser*/

	/*fill mat adv gasser*/

	void matrix::fill_matrix_adv(string data,map<const string, matrix>matrices){

		for (int i = 0; i < data.length(); i++)
		{
			if ((data[i] >= 'A' && data[i] <= 'Z') || (data[i] >= 'a' && data[i] <= 'z'))
			{
				// replace letter with value string
				string letter = data.substr(i, 1);
				string new_str = matrices[letter].getString();
				data.replace(i, 1, new_str);
			}
		}
		//cout << "this: "<< data << endl; // end gasser edit


		// Aly
		// data will be like this "1.1 2 3.5; 9.6 5.2 4.7"
		// these are 2 rows and three columns ('; ' separates rows .. ' ' separates colums)
		// initialize using initialize function provided above then assign values
		if (num_rows){ // resetting
			values.clear();
			num_columns = 0;
			num_rows = 0;
		}
		int start = 0;
		int end;
		if (data[data.length() - 1] != ';') { data = data + ";"; }

		vector<double> row;
		for (unsigned int i = 0; i< data.length(); i++){
			if ((data[i] == ' '&& data[i - 1] != ';') || (data[i] == ';')){
				end = i;
				//ading previous num

				row.push_back(atof(data.substr(start, end - start).c_str()));
				start = i + 1;
				if (data[i] == ';'){
					this->values.push_back(row);
					row.clear();
					if (start<data.length()){
						if (data[start] == ' ') start++;
					}
				}
			}
		}
		this->num_rows = this->values.size();
		this->num_columns = this->values[0].size();
	}

	/*end fill mat adv gasser*/

	/* to make sure line doesn't start with a space*/

	void matrix::remove_spaces(string& s)
	{
		int frst_num;
		for (int f = 0; f < s.length(); f++)
		{
			if (s[f] != ' ') { frst_num = f; break; }
		}
		s = s.substr(frst_num);
	}

	/* to make sure line doesn't start with a space*/


	/* Cut matrix into elements then send it to solve*/
	string matrix::cut_mat_solve(string &mat_val)
	{
		string mat_elemnt;
		//int frst_num = mat_val.find_first_not_of(' ');
		//mat_val = mat_val.substr(frst_num);  //cout << "this::" << mat_val <<"kok"<< endl;
		//int spcchk = mat_val.find(' '); int semichck = mat_val.find(';');

			/* first element */
			int frstspace = mat_val.find_first_of(' ');
			if (frstspace == -1) { frstspace = mat_val.find_first_of(';'); }
			mat_elemnt = mat_val.substr(0, frstspace - 0); solve_elemnt(mat_elemnt);
			mat_val.replace(0, frstspace - 0, mat_elemnt); //cout << "this1::" << mat_val << endl;
			/* first element */

			/* last element */
			int lastspace = mat_val.find_last_of(' ');
			if (lastspace == -1) { lastspace = mat_val.find_last_of(';'); }
			mat_elemnt = mat_val.substr(lastspace + 1);		 solve_elemnt(mat_elemnt);
			mat_val.replace(lastspace + 1, mat_val.length() - lastspace - 1, mat_elemnt); //cout << "this2::" << mat_val << endl;
			/* last element */

			for (int i = 0; i < mat_val.length(); i++)
			{
				int cut1;
				if (mat_val[i] == ' ' || mat_val[i] == ';')
				{
					cut1 = i; 
					if (mat_val[i + 1] == ' ') mat_val = mat_val.erase(i + 1, 1);
						for (int j = cut1 + 1; j < mat_val.length(); j++)
						{
							int cut2;
							if (mat_val[j] == ' ' || mat_val[j] == ';')
							{
								cut2 = j;
								if (mat_val[j + 1] == ' ') mat_val = mat_val.erase(i + 1, 1);
								mat_elemnt = mat_val.substr(cut1 + 1, cut2 - cut1 - 1);
								solve_elemnt(mat_elemnt);
								mat_val.replace(cut1 + 1, cut2 - cut1 - 1, mat_elemnt); //cout << "this132::" << mat_val << endl;
								break;
							}

						}
				}
			}

		return (mat_val);
	}

	string matrix::solve_elemnt(string &mat_elemnt)
	{
		/* check if it's just a matrix name don't change it*/
		for (int i = 0; i < mat_elemnt.length(); i++)
		{
			if (((mat_elemnt[i] >= 'A' && mat_elemnt[i] <= 'Z') || (mat_elemnt[i] >= 'a' && mat_elemnt[i] <= 'z')))
			{
				int rond_brac = mat_elemnt.find('(');
				if (rond_brac == -1)
				{
					return (mat_elemnt);
				}
			}
		}

			//matrix ans;
			// ans = solve_adv(mat_elemnt); still to be done be aly
			//mat_elemnt = ans.getString();
			int x = rand();
			char substring[100];
			sprintf_s(substring, "%d", x);
			mat_elemnt = substring;
			return (mat_elemnt);
	}

	string matrix::remove_space_after_semis(string &mat_vals)
	{
		int frstnum=0;
		for (int i = 0; i < mat_vals.length(); i++)
		{
			if (mat_vals[i] == ';'&&mat_vals[i + 1] == ' ')
			{
				for (int k = i+1; k < mat_vals.length(); k++)
				{
					if (mat_vals[k] != ' ') {frstnum = k; break;}
				}

				mat_vals = mat_vals.erase(i + 1, frstnum - i-1);
			}
		}

		return mat_vals;
	}

	/* End - Cut matrix into elements then send it to solve*/

	/* Solve advanced by alY
	
	matrix matrix::solve_adv(string value)
	{
		this a funtion to send to solve if it's a number and return value as 1x1 matrix 
		or if it's matrix opperation it does it.
	}
	
	*/

	/* -----------------------------------------End Advanced File example------------------------------------------*/



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

    void matrix :: call(vector<string>&arr2,vector<double>&fix_arr1,int index,double result){
	int sec_element=index+1;
	        fix_arr1.erase(fix_arr1.begin() + sec_element);
			fix_arr1.erase(fix_arr1.begin() + index);
	       fix_arr1.insert(fix_arr1.begin() + index,result);
		   arr2.erase( arr2.begin() + index );
		   }
	

	matrix  matrix::Solve(string data){
	vector<string>arr1;
	vector<string>arr2;
    vector<double>fix_arr1;
	vector<string>::iterator it;
	int am=0;
	string first_element="";
	for(unsigned int i=0;i<data.length();i++){
		if (data[i]=='^'||data[i]=='*'||data[i]=='/'||data[i]=='+'||data[i]=='-'){
		first_element=data.substr(0,i);	
		arr1.push_back(first_element);
		 am=i;
		 string of_op=data.substr(i,1);
			arr2.push_back(of_op);
		 data=data.erase(0,am+1);
		 i=0;
		}}
	arr1.push_back(data);
	
	//for(unsigned int j=0;j<arr2.size();j++){
	//cout<<arr2[j]<<endl;}
	//cout<<nof_op;
	/*******fix arr1 ***********************************/
	double res_tri;
   for(unsigned int j=0;j<arr1.size();j++)
   {
	string a=arr1[j];//put every data in string a 
	if(a[0] != ' '){
	if((a[0] >= 'A' && a[0] <= 'Z') || (a[0] >= 'a' && a[0] <='z'))
	  {//check on first char 
		string ins=a.substr(4,a.find(')')-4);
		double inside=stod(ins); 
		if(a[0]=='s')
		{res_tri=sin(inside);
		fix_arr1.push_back(res_tri);
		}
		if(a[0]=='c')
		{ res_tri=cos(inside);
		fix_arr1.push_back(res_tri);
		}
		if(a[0]=='t')
		{ res_tri=tan(inside);
		fix_arr1.push_back(res_tri);
		}
      }
	else{//number
		double num=stod(a);fix_arr1.push_back(num);
        }
	}
	else
	{    if((a[1] >= 'A' && a[1] <= 'Z') || (a[1] >= 'a' && a[1] <='z'))
	  {//check on second char 
		string ins=a.substr(5,a.find(')')-5);
		double inside=stod(ins);
		if(a[1]=='s')
		{res_tri=sin(inside);
		fix_arr1.push_back(res_tri);
		}
		if(a[1]=='c')
		{ res_tri=cos(inside);
		fix_arr1.push_back(res_tri);
		}
		if(a[1]=='t')
		{ res_tri=tan(inside);
		fix_arr1.push_back(res_tri);
		}
      }
	else{//number
		double num=stod(a);fix_arr1.push_back(num);
        }
	}
	}
//for(unsigned int k=0;k<fix_arr1.size();k++){
	//cout<<fix_arr1[k]<<endl;}
//cout<<endl;
/************************operations*******************/
   while(arr2.size()>0){
	   if (find(arr2.begin(), arr2.end(), "^") != arr2.end() )
	   {
		   it=find(arr2.begin(), arr2.end(), "^");
		   int pos = distance(arr2.begin(), it);
		   double part_result=pow(fix_arr1[pos],fix_arr1[pos+1]);
		   call(arr2,fix_arr1,pos,part_result);

	   }
	  
    else if(find(arr2.begin(), arr2.end(), "*") != arr2.end() )
	  { it=find(arr2.begin(), arr2.end(), "*");
	   int pos = distance(arr2.begin(), it);
		   double part_result=fix_arr1[pos]*fix_arr1[pos+1];
		   //fix_arr1.push_back(part_result);
		    call(arr2,fix_arr1,pos,part_result);
	  }
	else if(find(arr2.begin(), arr2.end(), "/") != arr2.end() )
	  {
		   it=find(arr2.begin(), arr2.end(), "/");
		   int pos = distance(arr2.begin(), it);
		   double part_result=fix_arr1[pos]/fix_arr1[pos+1];
		  // fix_arr1.push_back(part_result);
		   call(arr2,fix_arr1,pos,part_result);
	  }
	 else if(find(arr2.begin(), arr2.end(), "+") != arr2.end() )
	 {
				it=find(arr2.begin(), arr2.end(), "+");
		        int pos = distance(arr2.begin(), it);
		        double part_result=fix_arr1[pos]+fix_arr1[pos+1];   
				call(arr2,fix_arr1,pos,part_result);
	 }
	 else{ if(find(arr2.begin(), arr2.end(), "-") != arr2.end() )
	 {
						 it=find(arr2.begin(), arr2.end(), "-");
		        int pos = distance(arr2.begin(), it);
		        double part_result=fix_arr1[pos]-fix_arr1[pos+1];
				call(arr2,fix_arr1,pos,part_result);
		          
	 }}
	   
   }
   for(unsigned int z=0;z<fix_arr1.size();z++){
	cout<<fix_arr1[z]<<endl;}
 double re=fix_arr1[0];
   matrix result;
     result.initialize(1,1);
     for(int i=0;i<1;i++){
      for(int j=0;j<1;j++){
     result.values[i][j] = re;
        }
     }

        return result;
 
}
	