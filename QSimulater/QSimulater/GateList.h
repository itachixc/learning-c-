#pragma once
#include<vector>
using std::vector;
using std::string;
typedef vector<string> gate_list;
class GateList
{
public:
    GateList();
    ~GateList();
    gate_list gl;
    int Interpret(Gate::QState&,gate_list&);

};

