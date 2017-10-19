#include <iostream>
#include <vector>
using namespace std;

class matrix {
private:
    int num_rows;
    int num_columns;
    vector< vector<int> > values; //2d dynamic array using vector class
public:
    matrix(int rows, int cols){// constructor
        this->num_rows = rows;
        this->num_columns = cols;
         //pushing values with zeros (initialization)
         for(int i=0 ; i<rows ; ++i){
            vector<int> row;
            for(int j=0 ; j<cols ; ++j){
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
                cout<<values[i][j]<<"    ";
            }
            cout<<endl;
        }
    }

};
int main(){
    matrix u(3,5);
    u.print_matrix();
    return 0;
}