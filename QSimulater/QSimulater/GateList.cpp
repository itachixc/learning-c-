#include "GateList.h"
#include<vector>
#include "QuantumGates.h"
#include<string>
#include<iostream>
using std::cout;
using std::endl;
GateList::GateList()
{
    head = NULL;
    length = 0;
}


GateList::~GateList()
{
    Gate_seq *temp;
    for (int i = 0; i < length; i++)
    {
        temp = head;
        head = head->next;
        delete temp;
    }
}
int GateList::Length()        //得到链表长度
{
    return length;
}
void GateList::InsertHead(Gate_node val)
{
    Insert(val, 0);
}
void GateList::Insert(Gate_node val,int pos)       //插入节点
{
    if (pos < 0)
    {
        cout << "pos must be positive" << endl;
        return;
    }
    int index = 1;
    Gate_seq *temp = head;
    Gate_seq *gate = new Gate_seq(val);       //分配空间
    if (pos == 0)
    {
        gate->next = temp;
        head = gate;
        length++;
        return;
    }
    while (temp != NULL && index < pos)
    {
        temp = temp->next;
        index++;
    }
    if (temp == NULL)
    {
        cout << "Insert failed" << endl;
        return;
    }
    gate->next = temp->next;
    temp->next = gate;
    length++;
}
void GateList::Remove(Gate_node val)
{
    int pos = find(val);
    if (pos == -1)
    {
        cout << "Delete failed" << endl;
        return;
    }
    if (pos == 1)
    {
        head = head->next;
        return;
    }
    int index = 2;
    Gate_seq *temp = head;
    while (index < pos)
        temp = temp->next;
    temp->next = temp->next->next;
    length--;
}
int GateList::Find(Gate_node val)
{
    /*
    Gate_seq *temp = head;
    int index = 1;
    while (temp != NULL)
    {
        if(temp->val.)
    }
    */
    return 1;
}
void GateList::Reverse()         //反转链表
{
    if (head == NULL)
        return;
    Gate_seq *curgate = head, *nextgate = head->next,*temp;
    while (nextgate != NULL)
    {
        temp = nextgate->next;
        nextgate->next = curgate;
        curgate = nextgate;
        nextgate = temp;
    }
    head->next = NULL;
    head = curgate;
}
void GateList::Print()
{
    if (head == NULL)
    {
        cout << "GateList is empty" << endl;
        return;
    }
    Gate_seq *temp = head;
    while (temp != NULL)
    {
        cout << temp->val.gatename << "," << temp->val.para_1 << endl;
    }
}
int GateList::Interpret(Gate::QState& psi,gate_list& gl) {
    for (gate_list::iterator i = gl.begin; i != gl.end; ++i)
        switch((*i))

}