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
	double cal_determin_sq(int num_rows); //bad
	matrix inverse_2(); //fast
	// calculate co-factor func
	matrix cal_cofactor(int num_rows); //needed to calculate inverse ,bad
	// find transpose func
    matrix transpose(matrix& fac, double r); //bad
	int check_zero_dete(); //bad
	int is_identify(double n); // check if n is odd multiple of pi/2
public:
    matrix();
    matrix(string values);
    void initialize(int rows, int cols); 
	
    int get_num_rows (){ return this->num_rows;}
    int get_num_columns (){  return this->num_columns; }
    void print_matrix();
    // tasks:
	double determinant_2(int n); // good fast
    void fill_matrix (string data);
    matrix add_matrix( matrix& m);
    matrix sub_matrix( matrix& m);
    matrix mult_matrix( matrix& m);
    matrix inverse_matrix();  //you can use it as it calls the good one
    matrix transpose_matrix();
    matrix div_matrix( matrix& m);
	matrix bitwisediv_matrix(matrix &m); //solved issue
	matrix bitwisediv2_matrix(double c);//double c diveded by each value: c/val

	/*--------------------------------------phase 1 read file work----------------------------------------------*/
    static void run(string fpath);
	//only for phase1: reads the file and calls run_old_command

	static void run_old_command(string command, map<const string, matrix>& matrices);
	//processes given phase1 command(save&print)
	/*ex if command is: a = [1 2 3] then a will be saved to the map and printed to the screen
		and then called again with b = [3 5 6] same thing will happen,
		then called with c = a + b, c will be saved in map, and will be printed (all phase 1 operations supported)
	*/

    static void handle_read(map<const string, matrix>& matrices,string command,string name0,int op_index);
	//for phase 1
	/* to handle a line like this A = [2.2 7.3 4.8 2.4; 2.3 6.5 8.9 1.2;] -> only updates the map */
   
	static void decode(string command,string& name1,string& name2,int op_index);
	/*ex: C = A + B -> updates name1 & name2 given position of '+' operator (op_index)*/
   
	static void remove_back_slashes(string& s);//takes string and removes backSlashes from it
	/*--------------------------------------end of phase 1 read file work----------------------------------------------*/


	/*------------------------------------phase2 team1 work ---------------------------------------------------------*/
	static vector<int> get_braces_data(string data);
	/*
		get first good () positions:
		if string is ((7)) returns [1,3]
		if string is ()+() returns [3,4]
		if string is 1+2*4 returns [0] //only one element means no braces
		--note
		it ignores braces of log(),sin(),sqrt()...etc
	*/

	static void call(vector<string>&arr2,vector<double>&fix_arr1,int index,double result);
	//call is in partial_Solve function go to impelementation to see more comments

	static string partial_Solve(string data);//data is 5.5 + 12 * sin(0.4) + 2.2^4

	static matrix  Solve(string data);//data is (5.5 + 12) * (sin(0.4) + 2.2^4)

	static matrix column_by_column (matrix& a, matrix& b);//AYA
	/* ex: if a=  3    3     b= 1                c=a b  ->c= 3  3  1  2    // 2*2 2*1->2*3
	                                                           
	              4    4        6                            4  4  6  6
    */
	static matrix row_by_row(matrix& a,matrix& b); //AYA
	                  /* ex: if a=  3    3     b= 1  6           c=a   ->c= 3  3     // 2*2 1*2->3*2
					                4    4                         b        4  4
	                                                                        1  6
	              
    */                                                     





	matrix Sin();//AYA
	/* sin(B)
ans =

   0.932039   0.745705   0.818560
   0.963558   0.675463  -0.058374
  -0.993691   0.963558   0.998543
*/
	matrix Log();//AYA
    matrix Cos();//AYA
    matrix Tan();//AYA
	matrix Sqrt();//AYA
	/* sqrt(B)
ans =

   1.0954   1.5166   5.7964
   1.1402   1.5492   1.7889
   2.1448   1.1402   2.7928*/
	matrix Pow(int n);//AYA
	/*  B^2
ans =

   158.984    51.958   309.748
    19.400    12.910    76.318
    43.090    23.840   219.554
*/

static	matrix Rand(int a,int b);//AMERA
	/*rand(4,4)
	ans =

   0.506569   0.918578   0.739724   0.105990
   0.426771   0.614705   0.318234   0.448953
   0.510929   0.436886   0.610394   0.190750
   0.713867   0.698260   0.013628   0.215899


	*/

	static matrix Eye (int a,int b);//AMERA
	/*eye(4, 4)

Diagonal Matrix

   1   0   0   0
   0   1   0   0
   0   0   1   0
   0   0   0   1  */
	

	static matrix ones (int a,int b);//ALAA
	/*ones(3, 6)
	ans =
   1   1   1   1   1   1
   1   1   1   1   1   1
   1   1   1   1   1   1
	*/

	static matrix zeros(int r, int c); /* return a matrix of rXc of zeros */
	
	  string getString();
	  /* 
	  returns a string with the values of
	  matrix elements
	  */

	  matrix add_const(double a); //works for positive and negative
	  /*
	  adds double constant and
	  every element in the matrix then
	  returns the result as a matrix
	  */

	  matrix mult_const(double a); //works for mltiplying and dividing
		/*if you want to divide pass 1/a instead of a*/

	  matrix element_wise_power(double a);
	  /*
	  Every element in the matrix to
	  power of (double constant)
	  */

	  matrix strassen(matrix& u);
	  /*
	  strassen algorithm for multiplication is used
	  in the power function to optimize the code
	  */
	/* --------------------------------end of phase2 team1 work----------------------------------------------------------*/



	/* --------------------------------phase2 team2 work----------------------------------------------------------*/
	void fill_matrix_adv(string data,map<const string, matrix>& matrices);
	//only translates names to numbers then calls fill matrix to update this->values , num_rows and num_columns
	static void run_adv(string fpath);
	static void handle_read_adv(map<const string, matrix>& matrices, string command, string name0, int op_index);
	static void remove_spaces(string& s);
	static string remove_space_after_semis(string &mat_vals);
	static string cut_mat_solve(string &mat_val);
	static string solve_elemnt(string &mat_elemnt);
	static bool Isnt_num(string f);
	static matrix partial_Solve2 (string data);
	 //void call2(vector<string>&arr2,vector<matrix>&fix_arr1,int index,double result)

	/* --------------------------------end of phase2 team2 work----------------------------------------------------------*/
};

#endif