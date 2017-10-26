#include "QuantumGates.h"

#include <iostream>

#include <Windows.h>
const double PI = 3.141593;
using std::cout;
using std::endl;
using Gate::QState;

int main()
{
	Gate::QState psi;
	double P;
	psi.push_back(Gate::COMPLEX(0));
	psi.push_back(Gate::COMPLEX(0));
	psi.push_back(Gate::COMPLEX(0));
	psi.push_back(Gate::COMPLEX(0));
	psi.push_back(Gate::COMPLEX(1/sqrt(2)));
	psi.push_back(Gate::COMPLEX(0));
	psi.push_back(Gate::COMPLEX(1/sqrt(2)));
	psi.push_back(Gate::COMPLEX(0));
	//Gate::Hadamard(psi, 3);
	//Gate::Hadamard(psi, 3);
	//P=Gate::QMeasure(psi,2);
	Gate::SetQ(psi, 3, 0, 1);
	//Gate::COMPLEX a(3,4);
	for (QState::iterator iter = psi.begin(); iter != psi.end(); iter++)
		cout << *iter << endl;

	system("pause");

	return 0;

}