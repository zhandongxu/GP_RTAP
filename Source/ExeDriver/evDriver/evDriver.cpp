#include <iostream>
#include <tchar.h>
#include "..\..\DllProj\evNet\ev.h"
//#include <head.h>


string evdir = GetSourcePath();
// Two-OD small network
string twoodevNet(int n)
{
	ENET* net = new ENET(evdir + "\\Networks\\evNet\\2odNet\\2odNet");
	net->BuildTAPASNet();
	net->SetConv(1e-10);
	net->SetMaxIter(5000);
	net->SetMaxIterTime(25);
	//net->SetODLengthLimit(35);
	net->IniRandomDistByLowUpValue(300, 450, DistanceDist::lognormal);
	net->SetNumofClass(n);
	net->Solver(ENET::fd_multiple_gp);
	ostringstream a("");
	a << n << "," << TNM_FloatFormat(net->linkVector[0]->volume, 10, 4) << "," << net->linkVector[1]->volume << "," << net->linkVector[2]->volume << "," << net->linkVector[3]->volume;
	return a.str();
}
// Winnipeg network
void WinnipegNetFlow(float up_bound)
{
	ENET* net = new ENET(evdir + "\\Networks\\evNet\\winni\\winnipeg");
	net->BuildTAPASNet();
	net->SetConv(1e-6);
	net->SetMaxIter(2000);
	net->SetMaxIterTime(1000);
	//net->SetODLengthLimit(23);
	//net->Set_distance_Dist(10, 20, DistanceDist::uni);
	//net->InfinityDistance(1000000, 10000000, DistanceDist::uni);
	net->IniRandomDistByFacor(1,up_bound, DistanceDist::uni);
	net->Solver(ENET::fd_continuous_gp);
	string IterConvName = net->networkName + "Link flow pattern_contiConstr_GP" + to_string(up_bound) + ".conv";
	ofstream outfile;
	if (!TNM_OpenOutFile(outfile, IterConvName))
	{
		cout << "\n\tFail to prepare report: Cannot open .conv file to write convergence information!" << endl;
		return;
	}
	cout << "\tWriting iterative information into " << IterConvName << " for covergence info!" << endl;
	outfile << "Link ID,Link flow" << endl;
	for (int i = 0; i < net->numOfLink; i++)
	{
		outfile << i + 1 << "," << TNM_FloatFormat(net->linkVector[i]->volume, 10, 4) << endl;
	}
	outfile.close();
}

void WinnipegNetComput(float up_bound)
{
	ENET* net = new ENET(evdir + "\\Networks\\evNet\\winni\\winnipeg");
	net->BuildTAPASNet();
	net->SetConv(1e-10);
	net->SetMaxIter(2000);
	net->SetMaxIterTime(1000);
	net->IniRandomDistByFacor(1, up_bound, DistanceDist::uni);
	net->Solver(ENET::fd_continuous_gp); //fd_continuous_fw: FW algorithm
	net->ReportIter();
}

void WinnipegNetSingle()
{
	ENET* net = new ENET(evdir + "\\Networks\\evNet\\winni\\winnipeg");
	net->BuildTAPASNet();
	net->SetConv(1e-10);
	net->SetMaxIter(2000);
	net->SetMaxIterTime(1000);
	net->SetODLengthLimit(10000);
	net->Solver(ENET::fd_single_gp);
	net->ReportIter();
}

void winnipeg_multi(int n)
{
	ENET* net = new ENET(evdir + "\\Networks\\evNet\\winni\\winnipeg");
	net->BuildTAPASNet();
	net->SetConv(1e-10);
	net->SetMaxIter(5000);
	net->SetMaxIterTime(25);
	net->IniRandomDistByFacor(1, 3, DistanceDist::lognormal);
	net->SetNumofClass(n);
	net->Solver(ENET::fd_multiple_gp);
	net->ReportIter();
}


void main(int argc, char** argv)
{
	//winnipeg_multi(1); // Discrete model 
	/* WinnipegNetComput(): Computational performance of Continuous model
	 * WinnipegNetFlow(): Flow pattern of Continuous mdoel
	*/
	cout << "Run finisheds!!" << endl;
	int l = getchar();
}




