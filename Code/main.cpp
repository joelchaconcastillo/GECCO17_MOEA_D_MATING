#include <iostream>
#include <omp.h>
#include <ctime>
#include "MOEAD_Diversity.hpp"

using namespace std;
string DiccionarioProblemas(int Id);
void PrepareDirectory(string Dest, string Instance);
int strInstancetoId(string Instance);
void SetConfiguration(int argc, char*argv[], double &ProbCross, double &ProbMut, string &Path, int &SizePool, int &TotalGenerations, int &Sed, int &PeriodReport, string &Label, string &Instance, string &Crossover, string &Mutation);
int Str_to_Id_Crossover(string Crossover);
int Str_to_Id_Mutation(string Mutation);
void PrintHelp();
int main(int argc, char *argv[])
{
	if(argc<2)
         {
	    
	    cout << "Unknown Argument.."<<endl;
	    PrintHelp();
	    exit(0);
	 }

	int Sed = 1;
        double ProbCross= 0.9, ProbMut = 1.0/24.0;
        int SizePool = 100, TotalGenerations =250, PeriodReport=10;
	string Path="./", Instance="WFG1", Label ="0000", Crossover="SBX", Mutation="Polynomial";

	
	//Se realiza la lectura de la configuración que pertenece a los parámetros...
	SetConfiguration(argc, argv, ProbCross, ProbMut, Path, SizePool, TotalGenerations, Sed, PeriodReport, Label, Instance, Crossover, Mutation);
 	PrepareDirectory(Path, Instance);
	int IdCrossover = Str_to_Id_Crossover(Crossover);	
	int IdMutation = Str_to_Id_Mutation(Mutation);
        int i = strInstancetoId(Instance);
	clock_t TimeInit = clock();
	
	string ProblemName = DiccionarioProblemas(i);
	cout << " Iniciando...."<<ProblemName<<endl;
	Benchmark ObjBenchmark(i);
	MOEAD ObjMOEAD(SizePool, ProbCross, ProbMut,  ObjBenchmark, TotalGenerations, ProblemName, Sed, Path, Label, PeriodReport, IdCrossover, IdMutation);
	ObjMOEAD.Init_MOEAD();
	cout << "Problema "+ProblemName+" terminado.. TiempoTotal: " << (clock()-TimeInit)/ (double)  CLOCKS_PER_SEC << endl  ;
    
	return 0;
}

void PrintHelp()
{
	cout << "Instructions:"<<endl;
	cout << "--Instance NAMEINSTANCE (WFG1)"<<endl;
	cout << "--Sed SED (299)" <<endl;
	cout << "--Label OPERATOR_SBX, is a postfix by filename"<<endl;
	cout << "--Pcross 0.9, is the Crossover probability" <<endl;
	cout << "--Pmut 0.3, is the Mutation Probability " << endl;
	cout << "--Path ./RESULT, is the directory where will save results"<<endl;
	cout << "--SizePool 100, is the number of individual by generation"<<endl;
	cout << "--Generations 25000, is the number of generations"<<endl;
	cout << "--Crossover SBX, is the Crossover operator to use"<<endl;
	cout << "--Mutation Polynomial ,is the mutation operator to use"<<endl;
	cout << "--PeriodReport 10, each N generations will save informaton in the \"--Path\" files"<<endl;
	cout <<"---------Crossover Operators Availables-------------"<<endl;
	cout << "SBX, BLX, HBX, Fuzzy, SBXHybrid"<<endl;
	cout <<"---------Mutation Operators Availables--------------"<<endl;
	cout << "Polynomial, NormalDistribution, Random"<<endl;
}
void PrepareDirectory(string Dest, string Instance)
{
	string term = "mkdir -p ";
	term += Dest;
	term += "/"+Instance;
	system(term.c_str());
}
int Str_to_Id_Crossover(string Crossover)
{
	if(Crossover == "SBX")
		return _SBX;
	if(Crossover == "BLX")
		return _BLX;
	if(Crossover == "HBX")
		return _HBX;
	if(Crossover == "Fuzzy")
		return _FUZZY;
	if(Crossover == "SBXHybrid")
		return _SBXH;
	else
	{
		cout << "Unknown Crossover..."<<endl;
		exit(0);
	}
}
int Str_to_Id_Mutation(string Mutation)
{
	if(Mutation == "Polynomial")
		return _POLYNOMIAL;
	else  if(Mutation == "NormalDistribution")
		return _NORMAL;
	else if(Mutation == "Random")
		return _RANDOM;
	else
	{
		cout << "Unknown Mutation...."<<endl;
		exit(0);
	}
}
string DiccionarioProblemas(int Id)
{
	switch(Id)
	{
		case Type_WFG1:
			return "WFG1";
		break;
		case Type_WFG2:
			return "WFG2";
		break;
		case Type_WFG3:
			return "WFG3";
		break;
		case Type_WFG4:
			return "WFG4";
		break;
		case Type_WFG5:
			return "WFG5";
		break;
		case Type_WFG6:
			return "WFG6";
		break;
		case Type_WFG7:
			return "WFG7";
		break;
		case Type_WFG8:
			return "WFG8";
		break;
		case Type_WFG9:
			return "WFG9";
		break;
	}
}
int strInstancetoId(string Instance)
{
	if(Instance == "WFG1" )
		return Type_WFG1;
	else if(Instance == "WFG2")
		return Type_WFG2;	
	else if(Instance == "WFG3")
		return Type_WFG3;	
	else if(Instance == "WFG4")
		return Type_WFG4;	
	else if(Instance == "WFG5")
		return Type_WFG5;	
	else if(Instance == "WFG6")
		return Type_WFG6;
	else if(Instance == "WFG7")
		return Type_WFG7;		
	else if(Instance == "WFG8")
		return Type_WFG8;	
	else if(Instance == "WFG9")
		return Type_WFG9;	
	else
	{
		cout << "Unknown instance... " << endl;
		exit(0);
	}
}
void SetConfiguration(int argc, char*argv[], double &ProbCross, double &ProbMut, string &Path, int &SizePool, int &TotalGenerations, int &Sed, int &PeriodReport, string &Label, string &Instance, string &Crossover, string &Mutation)
{
	for(int i = 1; i < argc ; i++)
    	{
		string Terminal(argv[i]);
		if( Terminal == "--Instance")
			Instance = string(argv[++i]);
		else if(Terminal == "--Sed")
			Sed = atoi(argv[++i]);
		else if(Terminal == "--Label")
			Label = string(argv[++i]);
		else if(Terminal == "--Pcross")
			ProbCross = atof(argv[++i]);
		else if(Terminal == "--PMut")
			ProbMut = atof(argv[++i]);
		else if(Terminal == "--Path")
			Path = string(argv[++i]);
		else if(Terminal =="--SizePool")
			SizePool = atoi(argv[++i]);
		else if(Terminal == "--Generations")
			TotalGenerations = atoi(argv[++i]);
		else if(Terminal == "--Crossover")	
			Crossover = string(argv[++i]);
		else if(Terminal == "--Mutation")
			Mutation = string(argv[++i]);
		else if(Terminal =="--PeriodReport")
			PeriodReport = atoi(argv[++i]);

		else if(Terminal == "--help" || Terminal == "--h")
			PrintHelp();
		else
		{
			cout << Terminal<<endl;
			cout << "Unknown Argument...";
			exit(0);
		}
	    }
}
