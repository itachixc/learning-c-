#pragma once
/*

class QuantumGates
{
public:
	QuantumGates();
	~QuantumGates();
};

*/


#include<vector>
#include<complex>

namespace Gate
{

	using std::vector;
	using std::complex;
	using std::string;
	typedef vector<complex<double>> QState;                                
	typedef complex<double> COMPLEX;

	void Hadamard(QState&, size_t);
	void RX(QState&, size_t,double);
	void RY(QState&, size_t, double);
	void RZ(QState&, size_t, double);
	void CNOT(QState&, size_t,size_t);
	void CR(QState&, size_t, size_t,double);
	void iSWAP(QState&, size_t, size_t);
	void sqiSWAP(QState&, size_t, size_t);
	double QMeasure(QState&, size_t);                   //返回测量的比特处于0的概率
	void SetQ(QState&, size_t, double,int);
	string::iterator fgetw(string::iterator);
}
