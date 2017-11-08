#include "QuantumGates.h"
#include<process.h>
#include <iostream>
#include<string>
#include <Windows.h>
#include<time.h>
#include<stdio.h>
const double PI = 3.141593;
using std::cout;
using std::endl;
using std::string;
using Gate::QState;

int main()
{
	Gate::QState psi,matrix;
    //SHOR'S ALGORITHM
    string s("00000000000100");
    psi = Gate::Init_State(s);
    double t1 = clock();
    for(int i=0;i!=100;i++)
    Gate::Hadamard(psi, 1);
    double t2 = clock();
    cout << t2 - t1 << endl;
    /*   Gate::Hadamard(psi, 2);
    Gate::Hadamard(psi, 3);
    Gate::CR(psi, 3,4, 3.1415926);
    Gate::CR(psi, 3,5, 3.1415926);
    Gate::CR(psi, 3,6, 3.1415926);
    Gate::CR(psi, 3,7, 3.1415926);
    Gate::cswap(psi, 3,5,6);
    Gate::cswap(psi, 3, 4, 6);
    Gate::cswap(psi, 3, 4, 7);
    Gate::cswap(psi, 2, 4, 6);
    Gate::cswap(psi, 2, 5, 7);
    Gate::CNOT(psi,1,3); Gate::CNOT(psi, 3, 1); Gate::CNOT(psi, 1, 3);
    Gate::Hadamard(psi, 3);
    Gate::CR(psi, 3, 2, -3.1415926*1.5);
    Gate::Hadamard(psi, 2);
    Gate::CR(psi, 3, 1, -3.1415926*0.25);
    Gate::CR(psi, 3, 1, -3.1415926*1.5);
    Gate::Hadamard(psi, 1);
    */
	//double P;
    /*
	psi.push_back(Gate::COMPLEX(0));
	psi.push_back(Gate::COMPLEX(0));
	psi.push_back(Gate::COMPLEX(0));
	psi.push_back(Gate::COMPLEX(0));
	psi.push_back(Gate::COMPLEX(0));
	psi.push_back(Gate::COMPLEX(1));
	psi.push_back(Gate::COMPLEX(0));
	psi.push_back(Gate::COMPLEX(0));
    */
    //Gate::Rand_gen();
	//Gate::Hadamard(psi, 3,0.5);
   // Gate::RX(psi, 2,3.1415926, 0.5);
   // srand(time(NULL));
   // int rand_seed = Gate::Rand_seed_Init();
   // for (int i = 1; i != 10; i++)
        //cout <<Gate:: Rand_gen()<< endl;
        //Gate::Rand_gen(); Gate::Rand_gen();
	//Gate::Hadamard(psi, 3,0.5);
   // Gate::cswap(psi, 1,2,3 );
	//P=Gate::QMeasure(psi,2);
	//Gate::iSWAP(psi, 1,2);
	//Gate::COMPLEX a(3,4);
    /*****************************************************************
    test Init_State
    *****************************************************************/
    /*
    Gate::QState psi_1;
    string s("11001");
    psi_1 = Gate::Init_State(s);*/
    //long RandNum::
	//for (QState::iterator iter = psi.begin(); iter != psi.end(); iter++)
	//	cout << *iter << endl;
        
    /*****************************************************************
    test Usqg(QState&, size_t, QState&)
    *****************************************************************/
    //mat_reverse
    /*
    Gate::COMPLEX mat[4][4] = { 1,0,0,0,
                                0,1,0,0,
                                0,0,0,1,
                                0,0,1,0 };
    Gate::COMPLEX mat_rev[4][4];
    Gate::Mat_Reverse(mat, mat_rev);
    for (size_t i = 0; i != 4; i++)
    {
        for(size_t j=0;j!=4;j++)
            cout << mat_rev[i][j] << "\t";
        cout << "\n";
    }
    //finish
    */
	system("pause");
	return 0;

}