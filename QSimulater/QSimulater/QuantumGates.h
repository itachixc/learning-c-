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
    void Hadamard(QState& psi, size_t qn, double error_rate = 0);
	void RX(QState&, size_t,double, double error_rate = 0);
	void RY(QState&, size_t, double, double error_rate = 0);
	void RZ(QState&, size_t, double, double error_rate = 0);
	void CNOT(QState&, size_t,size_t, double error_rate = 0);
	void CR(QState&, size_t, size_t,double, double error_rate = 0);
	void iSWAP(QState&, size_t, size_t, double error_rate = 0);
	void sqiSWAP(QState&, size_t, size_t, double error_rate = 0);
    void cswap(QState&, size_t, size_t, size_t, double error_rate = 0);
	double GetProb(QState&, size_t);                   //返回测量的比特处于0的概率
	void QMeasure(QState&, size_t, double);
    void Usqg(QState&, size_t, QState&, double error_rate = 0);                //unitary single qubit gate 
	string::iterator fgetw(string::iterator);
    QState Init_State(string&);
    void Mat_Reverse(COMPLEX (*)[4],COMPLEX (*)[4]);    
    double Rand_gen();
}
