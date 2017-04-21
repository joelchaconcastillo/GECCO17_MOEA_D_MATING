/**
    Autor: Joel Chacón Castillo
    Fecha: 01/10/2016
    Descripción:
	This is the principal script representing the MOEA/D-EVSD
**/
#include <fstream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <cmath>
#include <string>
#include <algorithm>
#include "Benchmark.hpp"
#include "Individual.hpp"
#include "MOEAD_Diversity.hpp"

MOEAD::MOEAD(int SizePool, double ProbCross, int NBitsMut, Benchmark &ObjBenchmark, int NumberGenerations, string ProblemName, int Semilla, string Ruta, string Label, int PeriodReport, int IdCrossover, int IdMutation )
{
    srand(Semilla);
    this->T = 3;
    this->SizePool = SizePool;
    this->NumberObjectives = ObjBenchmark.getNObjectives();
    this->Dimension = ObjBenchmark.getDimension();
    this->ProbCross = ProbCross;
    this->NBitsMut = NBitsMut;
    this->ObjBenchmark = &ObjBenchmark;
    this->NumberGenerations = NumberGenerations;
    this->CurrentGeneration=0;
    this->IdCrossover = IdCrossover;
    this->IdMutation = IdMutation;
    this->PeriodReport = PeriodReport;
    this->VariableRepresentation = BINARY_ENCODE;
    this->FilenameGenotype = Ruta+"/"+ProblemName+"/Genotype_"+Label+"_"+to_string(Semilla)+".txt";
    this->FilenameFenotype =  Ruta+"/"+ProblemName+"/Fenotype_"+Label+"_"+to_string(Semilla)+".txt";
    this->FilenameFitness = Ruta+"/"+ProblemName+"/Fitness_"+Label+"_"+to_string(Semilla)+".txt";
    this->FileNameSummary = Ruta+"/"+ProblemName+"/SummaryFront_"+Label+"_"+to_string(Semilla)+".txt";
    this->ProblemName = ProblemName;
	//Reemplazar el archivo en caso de que exista
	ofstream SummaryGen, SummaryFen, SummaryFit;
	SummaryGen.open (this->FilenameGenotype);
	SummaryFen.open (this->FilenameFenotype);
	SummaryFit.open (this->FilenameFitness);
	SummaryGen << "Total_Generations\t" << this->NumberGenerations<< "\tNumber_Objectives\t"<< this->NumberObjectives <<endl;
	SummaryFen << "Total_Generations\t" << this->NumberGenerations<< "\tDimension\t"<< this->Dimension <<endl;
	SummaryFit << "Total_Generations\t" << this->NumberGenerations<< "\tNumber_Objectives\t"<< this->NumberObjectives <<endl;


}
/**
    Constructor del MOEAD para problemas continuos
**/
MOEAD::MOEAD(int SizePool, double ProbCross, double ProbMut ,Benchmark &ObjBenchmark, int NumberGenerations, string ProblemName, int Semilla, string Ruta, string Label, int PeriodReport, int IdCrossover, int IdMutation )
{
    srand(Semilla);

    this->T= 0.01*SizePool;	
    this->WAG =200; 
    this->ExplotationStage = 0.9;
    this->SizePool = SizePool;
    this->NumberObjectives = ObjBenchmark.getNObjectives();
    this->Dimension = ObjBenchmark.getDimension();
    this->ProbCross = ProbCross;
    this->ProbMutation = ProbMut;
    this->ObjBenchmark = &ObjBenchmark;
    this->NumberGenerations = NumberGenerations;
    this->CurrentGeneration=0;
    this->IdCrossover = IdCrossover;
    this->IdMutation = IdMutation;
    this->PeriodReport = PeriodReport;
    this->VariableRepresentation = REAL_ENCODE;
    this->FilenameGenotype = Ruta+"/"+ProblemName+"/Genotype_"+Label+"_"+to_string(Semilla)+".txt";
    this->FilenameFenotype =  Ruta+"/"+ProblemName+"/Fenotype_"+Label+"_"+to_string(Semilla)+".txt";
    this->FilenameFitness = Ruta+"/"+ProblemName+"/Fitness_"+Label+"_"+to_string(Semilla)+".txt";
    this->FileNameSummary = Ruta+"/"+ProblemName+"/SummaryFront_"+Label+"_"+to_string(Semilla)+".txt";
    this->ProblemName = ProblemName;
	//Eliminar el archivo en caso de que existan
	ofstream SummaryGen, SummaryFen, SummaryFit;

	SummaryFen.open (this->FilenameFenotype);
	SummaryFit.open (this->FilenameFitness);

;
}
/**
  Destructor de la clase
**/
MOEAD::~MOEAD()
{
}
/**
    En esta función se ejecuta el bloque principal del algoritmo
*/
void MOEAD::Init_MOEAD()
{
    //STEP 1) Initialization....
	//1.2 Compute Euclidean distances between any two weight vectors...
	InitUniformWeight();	
	ComputeDirectionVectors();

   	IndividualChecked = vector<bool>(SizePool, false); 
	InitNeighborhood();
	//1.3 Initialize population
        InitializePool(Pool_P);
	//1.4 Initialize z = (z_1, ..., z_m)
	InitIdealPoint();
	//Improve();
        while(CriterionStop())
        {

	 this->Neightborhood.clear();
	    InitNeighborhood();
	  vector<Individual> tmp = Pool_P;	
	   //2) Update....
	   for(int k = 0; k < this->SizePool; k++)
	   {
		  int i =k;// rand()%SizePool;
		   //2.1) Reproduction... selection, mutation polynomial, crossover SBX
		   Individual X = tmp[i]; //Initialize with same configuration...
		   Reproduction(X, i);
		
		   //2.2) Improvement...
		   if( !(CurrentGeneration % 10) && CurrentGeneration > NumberGenerations * ExplotationStage)
		   LocalSearch(X, i); //Disabled
		   if( CurrentGeneration > NumberGenerations * ExplotationStage) this->T = this->SizePool*0.1;
		   //2.3) Update of z
		   UpdateReference(X);
		   //2.4) Update Neightboring Solutions.
		   UpdateNeighBors(X, i);
	
	   }

	  if( !(this->CurrentGeneration % PeriodReport))
	   {
		ExportIndividualsFile(this->Pool_P, this->FilenameGenotype, this->FilenameFenotype, this->FilenameFitness, this->CurrentGeneration);
	   }
	    
        }
        End();
}

double MOEAD::DistanceVariable(Individual &A, Individual &B)
{
   double Dist = 0;
   for(int i = 0; i < A.getDimension(); i++)
	Dist +=(A.getVariable(i)-B.getVariable(i))*(A.getVariable(i)-B.getVariable(i))  ;
   return sqrt(Dist);
}
/**
	Here is applied the variation of mating streategy
	Note:
		In this code the Delta value is decremented from 1.0 to 0.5, and represent the probability 
		of select global individuals at difference that was indicated in the paper from GECCO17 where is 
		the probability of select a neighbour, then in this code Delta is decremented while in the paper Delta
		is increased. Anyway the process and results are the same.
*/
void MOEAD::Reproduction(Individual &ind, int IndexProblem) // Here add especific mating restriction....
{
	vector<int> Mating;
	double Delta = 1.0 - (double)CurrentGeneration / NumberGenerations;
	int SizeMating =SizePool*0.1;//SizePool;// max(0.2*SizePool, SizePool*( ((double) ExplotationStage*NumberGenerations-CurrentGeneration)/NumberGenerations  ));
	for(int i = 0; i < SizeMating; i++)
	{
		int r = rand()%((int)Neightborhood[IndexProblem].size());
		r = Neightborhood[IndexProblem][r];
		int rg = rand()%SizePool;	
		int indexparent; 
		if(getRandom(0.0,1.0) < Delta)
	  	 indexparent = rg;
		else 
		 indexparent = r;	
		Mating.push_back(indexparent);
	}
	double MaxDist = -10000;
	int p1, p2;
	if( CurrentGeneration < NumberGenerations*ExplotationStage )
	//get the two individuals so far....
	for(int i = 0; i < Mating.size(); i++)
	{
	   int a =  Mating[i];
	   for(int j = i+1; j < Mating.size(); j++)
	   {
		int b = Mating[j];
		if( DistanceVariable(Pool_P[a], Pool_P[b]) > MaxDist )	
		{
		   MaxDist = DistanceVariable(Pool_P[a], Pool_P[b]);
		   p1 = a;
		   p2 = b;
		}
	   }
	}
	else
	{
	   p1 = Mating[rand()%Mating.size()];
	   p2 = Mating[rand()%Mating.size()];
	}
		
	Individual Trash = ind;//Initialize with same configuration...
	SBXIndividualHybrid(Pool_P[p1], Pool_P[p2], ind, Trash ,this->ProbCross);
	if(getRandom(0,1.0) >=0.5)
	ind = Trash;
	int IndexDistributionMutation = 20;
	PolynomialMutationIndividual(ind, IndexDistributionMutation, this->ProbMutation);
	ind.EvalIndividual();
}
void MOEAD::LocalSearch(Individual &ind, int k)
{
	double CR = ProbCross;	
	double F = 0.5;
	for(int i = 0; i < 20; i++)
	{
		
		int p1 = Neightborhood[k][rand()%Neightborhood[k].size()];
		int p2 = Neightborhood[k][rand()%Neightborhood[k].size()];
		Individual Nuevo = ind;
	   int j = rand()%ind.getDimension();
	   for(int d = 0; d < ind.getDimension(); d++)
	   {
		if(getRandom(0,1.0) < CR || d == j)
	 	Nuevo.TruncateVariable(d, ind.getFenotype(d) + F *( Pool_P[p1].getFenotype(d) - Pool_P[p2].getFenotype(d) ));
		else
	 	Nuevo.TruncateVariable(d, ind.getFenotype(d));
	   }
	   Nuevo.EvalIndividual();
	   double F1 =  getFitness(Nuevo, Lambda[k] );
	   double F2 =  getFitness(ind, Lambda[k] );
	
	   if(  F1 < F2 )
	   {
		ind = Nuevo;
	   }
	}
	
}
void MOEAD::UpdateExternalPopulation()
{
	//Eliminate Individual dominated by current 
	for(int j = 0; j < Pool_P.size(); j++)
	{
	    for(int i = EP.size()-1; i>=0; i--)
		if( Pool_P[j].Dominate(EP[i]))
		   EP.erase(EP.begin()+i);	
	}  
	for(int j = 0; j < Pool_P.size(); j++)
	{
	   bool flag = true;
	   for(int i = 0; i < EP.size(); i++)
	   {
		if(EP[i].Dominate(Pool_P[j])  )
		{
		   flag = false;
		  break; 
		}
	   }
	  if(flag) EP.push_back(Pool_P[j]);
	}
}

void MOEAD::UpdateNeighBors(Individual &ind, int IdProblem)
{
   double suiteable = 10000000; 
   //For each individual of the Neightborhood	
   for(int i = 0; i < this->Neightborhood[IdProblem].size(); i++)
   {
	int k = Neightborhood[IdProblem][i];
	double f1 = getFitness(Pool_P[k], Lambda[k]);
	double f2 = getFitness(ind, Lambda[k]);
	vector<double>v1 = Pool_P[k].getAllFenotype(), v2 = ind.getAllFenotype();
	if(f2 < f1 )
	{
	   this->Pool_P[k] = ind;
	   f1 = f2;
	}
   }
}

void MOEAD::UpdateReference(Individual &ind)
{
   for(int i = 0; i < this->NumberObjectives; i++)
	this->Z_ref[i] = min(this->Z_ref[i], ind.getObjectiveValue(i));
}
double MOEAD::DistVector(vector<double> &A, vector<double> &B)
{
	double Suma = 0;	
	for(int i = 0; i < A.size(); i++)
	{
	   Suma += (A[i]-B[i])*(A[i]-B[i]); //Avoid use pow....
	}
   return sqrt(Suma);
}
void MOEAD::FastSort(vector<double> &X, vector<int> &Index, int T)
{
   for(int i = 0; i < T; i++)
   {
	for(int j = i + 1; j < this->SizePool; j++)
	{
	   if(X[i] > X[j])
	   {
	      //Swap values...	
		iter_swap(X.begin()+i, X.begin()+j);
	      //Swap index...
		iter_swap(Index.begin()+i, Index.begin()+j);
	   }
	}
   }
}
/**
   With two objectives this complexity could be O(N log(N) ) but now is O(n^2)
**/
void MOEAD::InitNeighborhood()
{
   this->Neightborhood.resize(this->SizePool, vector<int>(T, 0));	
    vector<double> X(this->SizePool);
    vector<int> Index(this->SizePool);
   for(int i = 0; i < this->SizePool; i++)
   {
	//Get all euclidean distances of lamnda...
	for(int j = 0; j < this->SizePool; j++)
	{
	   X[j] = DistVector(this->DirectionVectors[i], this->DirectionVectors[j]);
	   Index[j] = j;
	}
	//Sort and get the first T index...
	FastSort(X, Index, this->T);	
	for(int k = 0; k < this->T; k++)
	this->Neightborhood[i][k] = Index[k];
   }
}
double MOEAD::DistanceObjective(Individual &A, Individual &B)
{
   double Dist = 0;
   for(int i = 0; i < this->NumberObjectives; i++)
	Dist +=(A.getObjectiveValue(i)-B.getObjectiveValue(i))*(A.getObjectiveValue(i)-B.getObjectiveValue(i))  ;
   return sqrt(Dist);
}
/**
   Tchebycheff Approach
**/
double MOEAD::getFitness(Individual &ind, vector<double> &SubLambda)
{
   double Max = -1.0e+30;
   for(int i = 0; i < this->NumberObjectives; i++)
   {
	double Difference = fabs(ind.getObjectiveValue(i) - this->Z_ref[i]  );
	double Eval;
	if(SubLambda[i] == 0)
		Eval = 0.0001 * Difference;
	else 
		Eval = Difference * SubLambda[i];
	Max = max(Max, Eval);
   }
   return Max; 
}
void MOEAD::Transformation_WS(vector<double> &W, vector<double> &L)
{
   double Sum = 0;
   double SumZero = 0;
   double prod = 1;
   double epsilon = 0.0001;
   for(int i = 0; i < (int)W.size(); i++)
   {
		Sum += (1.0/W[i]);
		SumZero += (1.0/ (W[i]+epsilon));
		prod *=W[i];
   }
   for(int i = 0; i < (int)W.size(); i++)
   {
	if(!prod)
	L[i] = (1.0/(W[i]+epsilon)) / SumZero;
	else
	L[i] = (1.0/W[i]) / Sum;
   }
}
void MOEAD::ComputeDirectionVectors()
{
	this->DirectionVectors.resize(this->SizePool, vector<double>(2,0));
	for(int i =0 ; i < this->SizePool; i++)
	{
	   Transformation_WS(Lambda[i], DirectionVectors[i]);
	}
}
void MOEAD::InitUniformWeight()
{
	if( this->NumberObjectives == 2)
	{
		this->Lambda.resize(this->SizePool, vector<double>(2,0));
		for(int i = 0; i < this->SizePool; i++)
		{
		   double a = ((double)i/ (this->SizePool-1));
		   this->Lambda[i][0] = a;
		   this->Lambda[i][1] = 1.0-a;	
		}
	}
	else{
		cout << "No soport more than 2 objectives...";
		exit(0);
	}
}
void MOEAD::InitIdealPoint()
{
   Z_ref.resize(this->NumberObjectives, 1.0e+30);
}
/**
    Agregar información al final del archivo...
**/
void MOEAD::InitializePool(vector<Individual> &Pool_P)
{
    Pool_P.resize(this->SizePool);
    for(int i = 0; i < this->SizePool; i++)
    {
        Pool_P[i].InitializeIndividual(this->ObjBenchmark); //Inside evaluate each objective function...
    }
}
double MOEAD::getRandom(double Min, double Max)
{
    double f = (double)rand() / RAND_MAX;
    return Min + f * (Max - Min);
}
bool MOEAD::CriterionStop()
{
    if(this->CurrentGeneration < this->NumberGenerations)
    {
        this->CurrentGeneration++;
        return true;
    }
    else
    {
        return false;
    }
}
void MOEAD::End()
{
  ofstream save;
  save.open (FileNameSummary);
  for(int i = 0; i < (int)this->Pool_P.size(); i++)
  {
    for(int m = 0; m < this->NumberObjectives; m++)
      save << this->Pool_P[i].getObjectiveValue(m)<< " ";
    save << endl;
  }
}
