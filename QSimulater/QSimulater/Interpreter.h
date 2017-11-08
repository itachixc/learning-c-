#pragma once
#include"QuantumGates.h"
#include<vector>
#include<string>
using namespace std;
using Gate::QState;
typedef vector<string>  QCode;
typedef vector<int>     CReg;
class Interpreter
{
public:
    Interpreter();
    ~Interpreter();
};
    /*
    void Init(QCode);
    (QState::iterator)Find_repeat_point();
    void InterpretNext(QState&, string::iterator);
    QCode Code;
    const int N=30;
    CReg c;


    
};
*/
/*
namespace Interpreter
{
    Init_Interpret 
}
*/

