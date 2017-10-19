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
    }
    int get_num_rows (){
        return this->num_rows;
    }
    int get_num_columns (){
        return this->num_columns;
    }

};
int main(){
    matrix u(3,5);
    cout<<u.get_num_columns()<<endl; //check
    return 0;
}