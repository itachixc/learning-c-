/***********************************************************************
Copyright:
Author:Xue Cheng
Date:2017-10-26
Description:definition of quantum logic gates
************************************************************************/
#include<stdio.h>
#include<iostream>
#include<complex>
#include<vector>
#include <cmath>
#include "QuantumGates.h"
#include<time.h>
#include<process.h>
using std::cout;
using std::endl;

/*************************************************
  Name:	Rand_gen()
  Description:  16807 random number generator
  Argin:    None
  Argout:	None
  return:	random number in the region of [0,1]
*************************************************/
double Gate::Rand_gen()
{   
    double x, y;
    int i, a = 16807, m = 2147483647, q = 127773, r = 2836;             //difine constant number in 16807 generator                                                                  
    time_t rawtime;   
    struct tm * timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);                           
    static int I0 = timeinfo->tm_year + 70 * (timeinfo->tm_mon          //generate seed of the generator based on time
                    + 1 + 12 * (timeinfo->tm_mday 
                    + 31 * (timeinfo->tm_hour + 23 * (timeinfo->tm_min 
                    + 59 * timeinfo->tm_sec))));
    static int I ;
    if (a*(I0%q) - r*(I0 / q) >= 0)
        I = a*(I0%q) - r*(I0 / q);
    else
        I = a*(I0%q) - r*(I0 / q) + m;
    I0 = I;
    return (double)I / m;      
   
}
/*************************************************
  Name:	Hadamard
  Description:Hadamard gate,the matrix is:[1/sqrt(2),1/sqrt(2);1/sqrt(2),-1/sqrt(2)]
  Argin:    psi 	reference of vector<complex<double>> which contains quantum state information.
            qn	qubit number that the Hadamard gate operates on.
            error_rate  the errorrate of the gate
  Argout:	psi     quantum state after the gate
  return:	None
*************************************************/
void Gate::Hadamard(QState& psi, size_t qn,double error_rate)
{    
    double b = Rand_gen();                                              //Rand_gen() generates a random number of [0,1]    
   // cout << b << endl;
    if (b > error_rate)
    {
        size_t step = pow(2, qn - 1);
        COMPLEX alpha, beta;
        for (size_t i = 0; i< psi.size(); i += step * 2)
        {
            for (size_t j = i; j<i + step; j++)
            {
                alpha = psi[j];
                beta = psi[j + step];
                psi[j] = (alpha + beta)*0.70710678;                     // in j,the goal qubit is in |0>
                psi[j + step] = (alpha - beta)*0.70710678;			    //in j+step,the goal qubit is in |1>
            }
        }
    }
}
/*************************************************
  Name:	RX
  Description:RX gate,quantum state rotates θ by x axis.The matric is:
              [cos(θ/2),-i*sin(θ/2);i*sin(θ/2),cos(θ/2)]
  Argin:	psi	    reference of vector<complex<double>> which contains quantum state information.
            qn	    qubit number that the Hadamard gate operates on.
            theta	rotation angle
            error_rate  the errorrate of the gate
  Argout:	psi     quantum state after the gate
  return:	None
*************************************************/
void Gate::RX(QState& psi, size_t qn, double theta, double error_rate)
{
    double b = Rand_gen();                                              //Rand_gen() generates a random number of [0,1]    
    cout << b << endl;
    if (b > error_rate)
    {
        COMPLEX alpha, beta;
        size_t step = pow(2, qn - 1);
        double Ct = cos(theta / 2);
        double St = sin(theta / 2);
        COMPLEX compl(0, 1);

        for (size_t i = 0; i != psi.size(); i += step * 2)
            for (size_t j = i; j != i + step; j++)
            {
                alpha = psi[j];
                beta = psi[j + step];
                psi[j] = alpha * Ct - compl*St*beta;
                psi[j + step] = beta * Ct - compl*St*alpha;
            }
    }  
}
/*************************************************
  Name:	RY
  Description:RY gate,quantum state rotates θ by y axis.The matric is
              [cos(θ/2),-sin(θ/2);sin(θ/2),cos(θ/2)]
  Argin:	psi	    reference of vector<complex<double>> which contains quantum state information.
            qn	qubit number that the Hadamard gate operates on.
            theta	rotation angle
            error_rate  the errorrate of the gate
  Argout:	psi     quantum state after the gate
  return:	None
*************************************************/
void Gate::RY(QState& psi, size_t qn, double theta,double error_rate)
{
    double b = Rand_gen();                                              //Rand_gen() generates a random number of [0,1]    
    if (b > error_rate)
    {
        size_t step = pow(2, qn - 1);
        double Ct = cos(theta / 2);
        double St = sin(theta / 2);
        COMPLEX alpha, beta;
        for (size_t i = 0; i != psi.size(); i += step * 2)
            for (size_t j = i; j != i + step; j++)
            {
                alpha = psi[j];
                beta = psi[j + step];
                psi[j] = alpha * Ct - St*beta;
                psi[j + step] = beta * Ct + St*alpha;
            }
    }   
}
/*************************************************
  Name:	RZ
  Description:RZ gate,quantum state rotates θ by z axis.The matric is 
              [1 0;0 exp(i*θ)]
  Argin:	psi	    reference of vector<complex<double>> which contains quantum state information.
            qn	    qubit number that the Hadamard gate operates on.
            theta	rotation angle
            error_rate  the errorrate of the gate
  Argout:	psi     quantum state after the gate
  return:	None
*************************************************/
void Gate::RZ(QState& psi, size_t qn, double theta, double error_rate)
{
    double b = Rand_gen();                                              //Rand_gen() generates a random number of [0,1]    
    if (b > error_rate)
    {
        size_t step = pow(2, qn - 1);
        double Ct = cos(theta);
        double St = sin(theta);
        COMPLEX alpha, beta;
        COMPLEX compl(Ct, -St);
        for (size_t i = 0; i != psi.size(); i += step * 2)
            for (size_t j = i; j != i + step; j++)
            {
                alpha = psi[j];
                beta = psi[j + step];
                psi[j] = alpha;
                psi[j + step] = compl*psi[j + step];
            }
    }    
}
/*************************************************
  Name:	CNOT
  Description:CNOT gate,when control qubit is |0>,goal qubit does flip,
              when control qubit is |1>,goal qubit flips.the matric is:
              [1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0]
  Argin:	psi 	reference of vector<complex<double>> which contains quantum state information.
            qn_1	control qubit number
            qn_2	goal qubit number
            error_rate  the errorrate of the gate
  Argout:	psi     quantum state after the gate
  return:	None
*************************************************/
void Gate::CNOT(QState& psi, size_t qn_1, size_t qn_2, double error_rate)
{
    double b = Rand_gen();                                              //Rand_gen() generates a random number of [0,1]    
    if (b > error_rate)
    {
        int step1 = pow(2, qn_1 - 1), step2 = pow(2, qn_2 - 1);
        bool mark = 1;
        COMPLEX alpha, beta;
        if (qn_1>qn_2) {                                                // control qubit numer is higher
            for (size_t i = 0; i<psi.size(); i = i + step1) {
                mark = !mark;
                for (size_t j = i; j<i + step1; j = j + 2 * step2)
                    for (size_t k = j; k<j + step2; k++) {
                        if (mark) {
                            alpha = psi[k];
                            beta = psi[k + step2];
                            psi[k] = beta;                              // k:control |1>,goal |0>
                            psi[k + step2] = alpha;					    // k+step:control |1>,goal |1>


                        }

                    }
            }
        }
        else {                                                          // control qubit numer is lower
            for (size_t i = 0; i<psi.size(); i = i + 2 * step2)
                for (size_t j = i; j<i + step2; j = j + step1) {
                    mark = !mark;
                    for (size_t k = j; k<j + step1; k++) {
                        if (mark) {
                            alpha = psi[k];
                            beta = psi[k + step2];
                            psi[k] = beta;
                            psi[k + step2] = alpha;
                        }

                    }
                }
        }
    }       
}
/*************************************************
  Name:	    CR
  Description:	CR gate,when control qubit is |0>,goal qubit does not rotate,
                when control qubit is |1>,goal qubit rotate θ by z axis.the matric is:
                [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 exp(i*θ)]
  Argin:	psi 	reference of vector<complex<double>> which contains quantum state information.
            qn_1	control qubit number
            qn_2	goal qubit number
            theta	roration angle
            error_rate  the errorrate of the gate
  Argout:	psi     quantum state after the gate
  return:	None
*************************************************/
void Gate::CR(QState& psi, size_t qn_1, size_t qn_2,double theta, double error_rate)
{
    double b = Rand_gen();                                              //Rand_gen() generates a random number of [0,1]    
    if (b > error_rate)
    {
        int step1 = pow(2, qn_1 - 1), step2 = pow(2, qn_2 - 1);
        bool mark = 1;
        double Ct = cos(theta);
        double St = sin(theta);
        COMPLEX alpha, beta;
        COMPLEX compl(Ct, St);
        if (qn_1>qn_2) {                                                // control qubit numer is higher
            for (size_t i = 0; i<psi.size(); i = i + step1) {
                mark = !mark;
                for (size_t j = i; j<i + step1; j = j + 2 * step2)
                    for (size_t k = j; k<j + step2; k++) {
                        if (mark) {
                            alpha = psi[k];
                            beta = psi[k + step2];
                            psi[k] = alpha;                             // k:control |1>,goal |0>
                            psi[k + step2] = beta*compl;				// k+step:control |1>,goal |1>
                        }

                    }
            }
        }
        else {                                                          // control qubit numer is lower
            for (size_t i = 0; i<psi.size(); i = i + 2 * step2)
                for (size_t j = i; j<i + step2; j = j + step1) {
                    mark = !mark;
                    for (size_t k = j; k<j + step1; k++) {
                        if (mark) {
                            alpha = psi[k];
                            beta = psi[k + step2];
                            psi[k] = alpha;
                            psi[k + step2] = beta*compl;
                        }

                    }
                }
        }
    }    
}
/*************************************************
  Name:  	iSWAP
  Description:	iSWAP gate,both qubits swap and rotate π by z-axis,the matric is:
                [1 0 0 0;0 0 -i 0;0 -i 0 0;0 0 0 1]
  Argin:	psi	    reference of vector<complex<double>> which contains quantum state information.
            qn_1	first qubit number
            qn_2	second qubit number
            error_rate  the errorrate of the gate
  Argout:	psi     quantum state after the gate
  return:	None
*************************************************/
void Gate::iSWAP(QState& psi, size_t qn_1, size_t qn_2, double error_rate)
{
    double b = Rand_gen();                                              //Rand_gen() generates a random number of [0,1]    
    if (b > error_rate)
    {
        int temp, step1 = pow(2, qn_1 - 1), step2 = pow(2, qn_2 - 1);
        if (qn_1 < qn_2) {                                              // iSWAP(qn_1,qn_2) is agree with 
                                                                        //iSWAP(qn_2,qn_1)

            temp = step1;
            step1 = step2;
            step2 = temp;
        }
        temp = step1 - step2;
        bool mark = 1;
        COMPLEX compl(0, -1.0);
        COMPLEX alpha, beta;
        for (size_t i = 0; i<psi.size(); i = i + 2 * step1) {
            for (size_t j = i + step2; j<i + step1; j = j + 2 * step2)
                for (size_t k = j; k<j + step2; k++) {
                    alpha = psi[k];
                    beta = psi[k + temp];
                    psi[k] = compl*beta;                                // k:|01>
                    psi[k + temp] = compl*alpha;						// k+temp:|10>			
                }
        }
    }
    
}

void Gate::cswap(QState& psi, size_t qn_1, size_t qn_2, size_t qn_3, double error_rate)
{
    double b = Rand_gen();                                              //Rand_gen() generates a random number of [0,1]   
    COMPLEX alpha, beta;
    if (b > error_rate)
    {
        int temp, step1 = pow(2, qn_1 - 1), step2 = pow(2, qn_2 - 1), step3 = pow(2, qn_3 - 1);
        for (size_t i = step1; i < psi.size(); i = i + 2 * step1)
            for (size_t j = i; j < i + step1; ++j) {
                if (j % (2 * step2) >= step2 && j % (2 * step3) < step3)
                {
                    alpha = psi[j]; psi[j] = psi[j - step2 + step3]; psi[j - step2 + step3] = alpha;
                }                
            }                                      

    }
}
/*************************************************
  Name:  	sqiSWAP
  Description:	sqiSWAP gate,both qubits swap and rotate π by z-axis,the matrix is:
                [1 0 0 0;0 1/sqrt(2) -i/sqrt(2) 0;0 -i/sqrt(2) 1/sqrt(2) 0;0 0 0 1]
  Argin:	psi	reference of vector<complex<double>> which contains quantum state information.
            qn_1	first qubit number
            qn_2	second qubit number
            error_rate  the errorrate of the gate
  Argout:	psi     quantum state after the gate  
  return:	None
*************************************************/
void Gate::sqiSWAP(QState& psi, size_t qn_1, size_t qn_2, double error_rate)
{
    double b = Rand_gen();                                              //Rand_gen() generates a random number of [0,1]    
    if (b > error_rate)
    {
        int temp, step1 = pow(2, qn_1 - 1), step2 = pow(2, qn_2 - 1);
        if (qn_1 < qn_2) {                                              // sqiSWAP(qn_1,qn_2) is agree with 
                                                                        //sqiSWAP(qn_2,qn_1)

            temp = step1;
            step1 = step2;
            step2 = temp;
        }
        temp = step1 - step2;
        bool mark = 1;
        COMPLEX compl(0, -1.0);
        COMPLEX alpha, beta;
        for (size_t i = 0; i < psi.size(); i = i + 2 * step1) {
            for (size_t j = i + step2; j < i + step1; j = j + 2 * step2)
                for (size_t k = j; k < j + step2; k++) {
                    alpha = psi[k];
                    beta = psi[k + temp];
                    psi[k] = alpha / sqrt(2) + compl*beta / sqrt(2);        // k:|01>
                    psi[k + temp] = compl*alpha / sqrt(2) + beta / sqrt(2);	// k+temp:|10>			
                }
        }
    }
}

/*************************************************
  Name:  	GetProb
  Description:	get the probability of the target qubit in state |0>
  Argin:	psi	    reference of vector<complex<double>> which contains quantum state information.
            qn	    qubit number of the measurement
  Argout:	None
  return:	the probability of the measuremnt qubit in state |0>
*************************************************/
double Gate :: GetProb(QState& psi, size_t qn)
{
    double P0(0);         
    size_t step = pow(2, qn - 1);
    for (size_t i = 0; i< psi.size(); i += step * 2)
    {
        for (size_t j = i; j<i + step; j++)
        {
            P0 += abs(psi[j])*abs(psi[j]);           
        }
    }
    return P0;
}


/*************************************************
Name:  	GetProb
Description:	get the probability of the target qubit in state |0>
Argin:	psi	    reference of vector<complex<double>> which contains quantum state information.
qn	    qubit number of the measurement
Argout:	None
return:	the probability of the measuremnt qubit in state |0>
*************************************************/
void Gate::QMeasure(QState& psi, size_t qn,double P0)          
{
    int REG = 0, step = pow(2, qn - 1);
    double b = Rand_gen();
    if (b > P0)
        REG = 1;
    if (REG == 0)
    {		 
        P0 = 1 / sqrt(P0);
        for (int i = 0; i<psi.size(); i = i + 2 * step)
            for (int j = i; j<i + step; j++) {
                psi[j] *= P0;
                psi[j + step] = 0;
            }
    }
    else
    {
        P0 = 1 / sqrt(1 - P0);
        for (int i = 0; i<psi.size(); i = i + 2 * step)
            for (int j = i; j<i + step; j++) {
                psi[j] =0;
                psi[j + step] *=P0;
            }
    }
    
}
Gate::QState Gate::Init_State(string& s)         //Init
{
    size_t n = size(s);
    size_t length = pow(2, n);
    Gate::QState psi_init(length, 0);
    size_t origin(0);
    for (size_t i = 0; i != n;i++) 
        if (s[i] == '1') {
            origin += pow(2, n - i - 1);
            psi_init[origin] = 1;
        }            
    return psi_init;  
}
void Gate::Usqg(QState& psi, size_t qn, QState& matrix, double error_rate )
{
    double b = Rand_gen();                                              //Rand_gen() generates a random number of [0,1]    
    if (b > error_rate)
    {
        COMPLEX det = matrix[0] * matrix[3] - matrix[1] * matrix[2];
        if (abs(abs(det) - 1) > 0.000001)
            return;
        size_t step = pow(2, qn - 1);
        COMPLEX alpha, beta;
        for (size_t i = 0; i< psi.size(); i += step * 2)
        {
            for (size_t j = i; j<i + step; j++)
            {
                alpha = psi[j];
                beta = psi[j + step];
                psi[j] = alpha*matrix[0] + beta*matrix[1];                        // in j,the goal qubit is in |0>
                psi[j + step] = alpha*matrix[2] + beta*matrix[3];;					//in j+step,the goal qubit is in |1>
            }
        }
    }
}
void  Gate :: Mat_Reverse(COMPLEX (*mat)[4],COMPLEX (*mat_rev)[4])
{
    COMPLEX  temp;
   /* mat_rev = { { mat[0][0],mat[0][2],mat[0][1],mat[0][3] },
                { mat[2][0],mat[2][2],mat[2][1],mat[2][3] },
                { mat[1][0],mat[1][2],mat[1][1],mat[1][3] },
                { mat[3][0],mat[3][2],mat[3][1],mat[3][3] }, };*/
    for (size_t j = 0; j != 4; j++) {                            //第二行第三行交换
        *(*mat_rev  + j) = *(*mat  + j);
        *(*(mat_rev + 1) + j) = *(*(mat + 2) + j);
        *(*(mat_rev + 2) + j) = *(*(mat + 1) + j);
        *(*(mat_rev + 3) + j) = *(*(mat + 3) + j);
    }
        
    for (size_t j = 0; j != 4; j++) {                            //第二列第三列交换
        temp = *(*(mat_rev + j) + 1);
        *(*(mat_rev + j) + 1) = *(*(mat_rev + j) + 2);
        *(*(mat_rev + j) + 2) = temp;        
    }
   // return mar_rev;

}














