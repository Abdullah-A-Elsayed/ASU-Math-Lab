#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
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

    void print_matrix(){ // print matrix for testing
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
		vector<float> row;
		for(int i = 0 ; i< data.length(); i++){
			if (data[i]==' '||data[i]==';'){
				end=i;
				row.push_back(atof(data.substr(start,end).c_str()));
				start=i+1;
				if(data[i]==';'){
					this->values.push_back(row);
					row.clear();
				}
			}
		}
		this->num_rows = this->values.size();
		this->num_columns = this->values[0].size();
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

	/*------------------------------------------------START GASSER inverse-mat--------------------------------------------------------------*/
	
	// fill matrix with test values
	void fill_mat_test()
	{
		this->values[0][0] = 3; this->values[0][1] = 0; this->values[0][2] = 2;
		this->values[1][0] = 2; this->values[1][1] = 0; this->values[1][2] = -2;
		this->values[2][0] = 0; this->values[2][1] = 1; this->values[2][2] = 1;
	}

	// calculate determine func
	float cal_determin_sq(int num_rows)
	{
		float s = 1, det = 0; 
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
	matrix cal_cofactor(int num_rows)
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
				fac.values[q][p] = pow(-1, q + p) * b.cal_determin_sq((f - 1));
			}
		}
		
		return (this->transpose(fac, f));
	}

	// find transpose func
	matrix transpose(matrix fac, float r)
	{
		int i, j;
		matrix b; b.initialize(this->num_rows, this->num_columns);
		matrix inverse; inverse.initialize(this->num_rows, this->num_columns);
		float d;

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

    // find inverse matrix
	matrix inverse_matrix()
	{
		float det_val = this->cal_determin_sq(this->num_rows); // calculate determine value for the matrix

		// intialize the inverse matrix with zeros & check if it is a square matrix
		matrix m; m.initialize(this->num_rows, this->num_columns); 
		if (this->num_rows != this->num_columns){ cout << "No inverse for non-square matrix" << endl; return (m); }

		// check if determine equals zero
		else if (det_val == 0){ cout << "No inverse for this zero-determine matrix" << endl; return (m); }

		// strat to get the inverse for the matrix
		else
		{
			m = this->cal_cofactor(this->num_rows);
			return (m);
		}

    }

	/*------------------------------------------------END GASSER inverse-mat--------------------------------------------------------------*/

    matrix div_matrix( matrix m){
        //Alaa Ayman
        // create a result matrix with correct dimensions then initialize it using initialize function provided above
        // result = this * m inversed
        // use previous functions
    }
};

int main()
{
	/* 
		// code test for inverse

		matrix x; matrix y;
		x.initialize(3, 3);
		x.fill_mat_test();
		x.print_matrix();
		cout << endl;
		y= x.inverse_matrix();
		cout << "Inverse Matrix\n" << endl;
		y.print_matrix();
		cout << endl;
	*/

    return 0;
}
