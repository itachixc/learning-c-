/***********************************************************************
Copyright:
Author:Xue Cheng
Date:2017-10-26
Description:difinition of quantum logic gates
************************************************************************/

#include<iostream>
#include<complex>
#include<vector>
#include <cmath>
#include "QuantumGates.h"




/*************************************************
  Name:	Hadamard
  Description:Hadamard gate,the matrix is:[1/sqrt(2),1/sqrt(2);1/sqrt(2),-1/sqrt(2)]
  Argin:	argin0	reference of vector<complex<double>> which contains quantum state information.
			argin1	qubit number that the Hadamard gate operates on.
  Argout:	None
  return:	None
*************************************************/
void Gate::Hadamard(QState& psi, size_t qn)
{

	size_t step = pow(2, qn - 1);
	COMPLEX alpha, beta;
	for (size_t i = 0; i< psi.size(); i += step*2)
	{
		for (size_t j=i;j<i+step;j++)                                         
		{
			alpha = psi[j];
			beta = psi[j+step];
			psi[j] = (alpha + beta)*0.70710678;                         // in j,the goal qubit is in |0>
			psi[j+step] = (alpha - beta)*0.70710678;					//in j+step,the goal qubit is in |1>
		}
	}

}
/*************************************************
  Name:	RX
  Description:RX gate,quantum state rotates θ by x axis.
  Argin:	argin0	reference of vector<complex<double>> which contains quantum state information.
            argin1	qubit number that the Hadamard gate operates on.
			argin2	rotation angle
  Argout:	None
  return:	None
*************************************************/
void Gate::RX(QState& psi, size_t qn, double theta)
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
void Gate::RY(QState& psi, size_t qn, double theta)
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
			psi[j + step] = beta * Ct +St*alpha;
		}
}
void Gate::RZ(QState& psi, size_t qn, double theta)
{
	
	size_t step = pow(2, qn - 1);
	double Ct = cos(theta);
	double St = sin(theta);
	COMPLEX alpha, beta;
	COMPLEX compl(Ct,-St);
	for (size_t i = 0; i != psi.size(); i += step * 2)
		for (size_t j = i; j != i + step; j++)
		{
			alpha = psi[j];
			beta = psi[j + step];
			psi[j] = alpha;
			psi[j + step] = compl*psi[j + step];
		}
}//rotate by z axis
void Gate::CNOT(QState& psi, size_t qn_1, size_t qn_2)
{
	
	int step1 = pow(2, qn_1 - 1), step2 = pow(2, qn_2 - 1);
	bool mark = 1;
	COMPLEX alpha, beta;
	if (qn_1>qn_2) {
		for (size_t i = 0; i<psi.size(); i = i + step1) {
			mark = !mark;     //判断i态中控制位为0还是1
			for (size_t j = i; j<i + step1; j = j + 2 * step2)
				for (size_t k = j; k<j + step2; k++) {         //k都是从目标位为0开始
					if (mark) {
						alpha = psi[k];
						beta = psi[k + step2];
						psi[k] = beta;
						psi[k + step2] = alpha;
						
						
					}//控制位为1
					
				}
		}
	}
	else {                            //k1<k2
		for (size_t i = 0; i<psi.size(); i = i + 2 * step2)
			for (size_t j = i; j<i + step2; j = j + step1) {
				mark = !mark;
				for (size_t k = j; k<j + step1; k++) {
					if (mark) {            //控制比特为0
						alpha = psi[k];
						beta = psi[k + step2];
						psi[k] = beta;
						psi[k + step2] = alpha;
					}
					
				}
			}
	}//k1<k2情况         */
	
}
void Gate::CR(QState& psi, size_t qn_1, size_t qn_2,double theta)
{

	int step1 = pow(2, qn_1 - 1), step2 = pow(2, qn_2 - 1);
	bool mark = 1;
	double Ct = cos(theta);
	double St = sin(theta);
	COMPLEX alpha, beta;
	COMPLEX compl(Ct, St);
	if (qn_1>qn_2) {
		for (size_t i = 0; i<psi.size(); i = i + step1) {
			mark = !mark;     //判断i态中控制位为0还是1
			for (size_t j = i; j<i + step1; j = j + 2 * step2)
				for (size_t k = j; k<j + step2; k++) {         //k都是从目标位为0开始
					if (mark) {
						alpha = psi[k];
						beta = psi[k + step2];
						psi[k] = alpha;
						psi[k + step2] = beta*compl;
					}//控制位为1

				}
		}
	}
	else {                            //k1<k2
		for (size_t i = 0; i<psi.size(); i = i + 2 * step2)
			for (size_t j = i; j<i + step2; j = j + step1) {
				mark = !mark;
				for (size_t k = j; k<j + step1; k++) {
					if (mark) {            //控制比特为0
						alpha = psi[k];
						beta = psi[k + step2];
						psi[k] = alpha;
						psi[k + step2] = beta*compl;
					}

				}
			}
	}//k1<k2情况         */

}
void Gate::iSWAP(QState& psi, size_t qn_1, size_t qn_2)           //[1 0 0 0;0 0 -i 0;0 -i 0 0;0 0 0 1]
{

	int step1 = pow(2, qn_1 - 1), step2 = pow(2, qn_2 - 1);
	bool mark = 1;
	COMPLEX compl(0, -1.0);
	COMPLEX alpha, beta;
	if (qn_1>qn_2) {
		for (size_t i = 0; i<psi.size(); i = i + step1) {
			mark = !mark;     //判断i态中控制位为0还是1
			for (size_t j = i; j<i + step1; j = j + 2 * step2)
				for (size_t k = j; k<j + step2; k++) {         //k都是从目标位为0开始
					if (mark) {			
						alpha = psi[k];
						beta = psi[k - step2];
						psi[k] = compl*beta;
						psi[k - step2] = compl*alpha;
					}//控制位为1					

				}
		}
	}
	else {                            //k1<k2
		for (size_t i = 0; i<psi.size(); i = i + 2 * step2)
			for (size_t j = i+step1; j<i + step2; j = j + 2*step1) {
				//mark = !mark;
				for (size_t k = j; k<j + step1; k++) {
					//if (!mark) {            //控制比特为1
						alpha = psi[k];
						beta = psi[k+step2-step1];
						psi[k] = compl*beta;
						psi[k + step2-step1] = compl*alpha;
				//	}
					

				}
			}
	}//k1<k2情况         */

}
void Gate::sqiSWAP(QState& psi, size_t qn_1, size_t qn_2)         //unfinished
{

	int step1 = pow(2, qn_1 - 1), step2 = pow(2, qn_2 - 1);
	bool mark = 1;
	COMPLEX alpha, beta;
	COMPLEX compl(0, -1.0);
	if (qn_1>qn_2) {
		for (size_t i = 0; i<psi.size(); i = i + step1) {
			mark = !mark;     //判断i态中控制位为0还是1
			for (size_t j = i+step2; j<i + step1; j = j + 2 * step2)
				for (size_t k = j; k<j + step2; k++) {         //k都是从目标位为0开始
					if (!mark) {          //mark==1,控制位为1
						alpha = psi[k];
						beta = psi[k +step2];
						psi[k] = alpha*0.70710678 + compl*beta*0.70710678;
						psi[k+step2]= beta*0.70710678 + compl*alpha*0.70710678;
					}//控制位为1					
				}
		}
	}
	else {                            //k1<k2
		for (size_t i = step2; i<psi.size(); i = i + 2 * step2)    //从目标为1开始
			for (size_t j = i; j<i + step2; j = j + step1) {
				mark = !mark;
				for (size_t k = j; k<j + step1; k++) {
					if (!mark) {            //控制比特为0
						alpha = psi[k];
						beta = psi[k + step2];
						psi[k] = alpha*0.70710678 + compl*beta*0.70710678;
						psi[k + step2] = beta*0.70710678 + compl*alpha*0.70710678;
					}					
				}
			}
	}//k1<k2情况         */

}
double Gate :: QMeasure(QState& psi, size_t qn)
{
	double P0(0);         //处于0态的概率
	size_t step = pow(2, qn - 1);
	for (size_t i = 0; i< psi.size(); i += step * 2)
	{
		for (size_t j = i; j<i + step; j++)
		{
			P0 += abs(psi[j])*abs(psi[j]);
			//psi[j] = (alpha + beta)*0.70710678;
			//psi[j + step] = (alpha - beta)*0.70710678;
			//psi[i].
		}
	}
	return P0;
}
void Gate::SetQ(QState& psi, size_t qn,double P0,int REG)          //PO为处于0态的概率
{
	int step = pow(2, qn - 1);
	//double Prob;
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
















