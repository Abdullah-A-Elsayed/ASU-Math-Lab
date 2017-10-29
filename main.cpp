#include <iostream>
#include <string>
#include "matrix.h"
using namespace std;
int main()
{
	// slides example
	matrix A;
    A.fill_matrix("1.4 2.2 3.2; 4.4 5.4 6.4; 3.3 4.2 2.2");
    cout<<"A is:"<<endl;
    A.print_matrix();
    cout<<endl;

    matrix B;
    B.fill_matrix("1.5 4.1 5.4; 3.1 4.2 1.2; 3.2 4.3 2.2");
    cout<<"B is:"<<endl;
    B.print_matrix();
    cout<<endl;

	try{
		matrix C = A.add_matrix(B);
		cout<<"summing: C"<<endl;
        C.print_matrix();
        cout<<endl;

		matrix D = A.sub_matrix(B);
		cout<<"subtracting: D"<<endl;
        D.print_matrix();
        cout<<endl;

		matrix E = A.mult_matrix(B);
		cout<<"multiplying: E"<<endl;
        E.print_matrix();
        cout<<endl;

        matrix F = A.div_matrix(B);
		cout<<"dividing: F"<<endl;
        F.print_matrix();
        cout<<endl;

        matrix G = A.inverse_matrix();
		cout<<"A' : G"<<endl;
        G.print_matrix();
        cout<<endl;
	}
	catch(string err){ cout<<err;}
    return 0;
}