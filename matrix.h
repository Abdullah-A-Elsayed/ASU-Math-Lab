#ifndef MAT_H

#define MAT_H
#include <string>
#include <map>
#include <vector>
using namespace std;
class matrix {
    //setup:
private:
    int num_rows;
    int num_columns;
    vector< vector<double> > values; //2d dynamic array of double using vector class
		// calculate determine func
	double cal_determin_sq(int num_rows);
	// calculate co-factor func
	matrix cal_cofactor(int num_rows);
	// find transpose func
    matrix transpose(matrix& fac, double r);
    static void handle_read(map<const string, matrix>& matrices,string command,string name0,int op_index);
    static void decode(string command,string& name1,string& name2,int op_index);
public:
    matrix();
    matrix(string values);
    void initialize(int rows, int cols);
    int get_num_rows (){ return this->num_rows;}
    int get_num_columns (){  return this->num_columns; }
    void print_matrix();
    // tasks:
    void fill_matrix (string data1);
    matrix add_matrix( matrix& m);
    matrix sub_matrix( matrix& m);
    matrix mult_matrix( matrix& m);
    matrix inverse_matrix();
    matrix transpose_matrix();
    matrix div_matrix( matrix& m);
    static void run(string fpath);
};

#endif
