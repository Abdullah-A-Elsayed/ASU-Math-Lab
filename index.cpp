#include <iostream>
#include <vector>
#include <string>
using namespace std;

class matrix {
    //setup:
private:
    int num_rows;
    int num_columns;
    vector< vector<float> > values; //2d dynamic array of float using vector class
public:
    matrix(){ // constructing ...
        this->num_rows = this->num_columns =0;
    }

    void initialize(int rows, int cols){ // taking dimensions
        this->num_rows = rows;
        this->num_columns = cols;
         //pushing values with zeros (initialization)
         for(int i=0 ; i<rows ; ++i){ //rows
            vector<float> row;
            for(int j=0 ; j<cols ; ++j){ //columns
                row.push_back(0);
            }
            values.push_back(row);
        }
    }

    int get_num_rows (){
        return this->num_rows;
    }

    int get_num_columns (){
        return this->num_columns;
    }

    void print_matrix(){ // for testing
        for(int i=0 ; i< this-> num_rows ; ++i){
            for(int j=0 ; j< this->num_columns ; ++j){
                cout<<this->values[i][j]<<"    ";
            }
            cout<<endl;
        }
    }


    // tasks:
    void fill_matrix (string data){
        // Aly
        // data will be like this "1.1 2 3.5;9.6 5.2 4.7"
        // these are 2 rows and three columns (';' separates rows .. ' ' separates colums)
        // initialize using initialize function provided above then assign values
	int start=0;
	int end;
	for(int i = 0 ; i<data.length ; i++){
	if (data[i]==' '||data[i]==';'){	end=i;
						row.push_back(data.substr(start,end));
						start=i+1;
						if(data[i]==';')
							{
							values.push_back(row);
							row.clear();
							}
				
					}
	
        }

	
    }

    matrix add_matrix( matrix m){
        //Aya
        // create a result matrix with correct dimensions then initialize it using initialize function provided above
        // result = this + m
        // return result
    }

    matrix sub_matrix( matrix m){
        //Do'aa
        // create a result matrix with correct dimensions then initialize it using initialize function provided above 
        // result = this - m
        // return result
    }

    matrix mult_matrix( matrix m){
        //Amira
        // create a result matrix with correct dimensions then initialize it using initialize function provided above
        // result = this * m
        // return result
    }

    matrix inverse_matrix(){
        //gasser
        //create a result matrix with correct dimensions then initialize it using initialize function provided above
        // result = this inverse
        // return result
        // good luck
    }

    matrix div_matrix( matrix m){
        //Alaa Ayman
        // create a result matrix with correct dimensions then initialize it using initialize function provided above
        // result = this * m inversed
        // use previous functions
    }
};
int main(){
    
    return 0;
}
