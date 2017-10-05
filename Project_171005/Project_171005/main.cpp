#include "TicketMachine.h"
#include <iostream>
int main()
{
	TicketMachine tm;
	tm.insertMoney(100);
	tm.showBalance();
	tm.showPrompt();
	system("pause");
	return 0;
}