#ifndef MOEAD_HPP_INCLUDED
#define MOEAD_HPP_INCLUDED
#include <vector>
#include <iostream>
#include <queue>
#include "Individual.hpp"
#include "EAOperators.hpp"
#include "EvolutiveAlgorithm.hpp"
/////Identificador de cada operador de cruza
#define _SBX 101
#define _BLX 102
#define _HBX 103
#define _FUZZY 104
#define _SBXH 105
/////identificador de cada operador de mutación....
#define _POLYNOMIAL 201
#define _NORMAL 202
#define _RANDOM 203
using namespace std;
class MOEAD : public EvolutiveAlgorithm
{
    public:

	MOEAD(int SizePool, double ProbCross, int NBitsMut, Benchmark &ObjBenchmark, int NumberGenerations, string ProblemName, int Semilla, string Ruta, string Label, int PeriodReport, int IdCrossover, int IdMutation  );

	MOEAD(int SizePool, double ProbCross, double ProbMut ,Benchmark &ObjBenchmark, int NumberGenerations, string ProblemName, int Semilla, string Ruta, string Label, int PeriodReport, int IdCrossover, int IdMutation  );

       MOEAD();
        ~MOEAD();
        void Init_MOEAD();
        /**
            Esta estructura define el criterio para ordenar
            un conjunto de individuos en base al valor de
            la función objetivo indicada por el índice proporcionado.
        **/
        struct LocalIndividual
        {
            LocalIndividual(int IndexObjective){this->IndexObjective = IndexObjective;}
            bool operator()(Individual A, Individual B) { return A.getObjectives().SpaceObjectives[IndexObjective].Fitness < B.getObjectives().SpaceObjectives[IndexObjective].Fitness;  }
            //bool operator()(int i, int j){return i < j;}
            int IndexObjective;
        };
        void PlotInterfaceRSpaceObjective();
        void PlotInterfaceRSpaceVariables();
	void PlotInterfaceRSpaceObjectiveEP();
        //void getSolutions(vector<vector<double>> &SetFront, vector<vector<double>> &SetParetoOptimal, vector<int> &IndexBoundsObjectives);
        //void WriteFilePool();
    private:
        Benchmark *ObjBenchmark;
        int SizePool, NumberObjectives, Dimension, TotalGenerations, CurrentGeneration, T;
        double ProbCross, ProbMutation;
	double ExplotationStage;
        string FileNameSummary, FileNameFront, FilenameGenotype, FilenameFitness, FilenameFenotype, ProblemName;
        int NBitsMut, IdCrossover, IdMutation, PeriodReport;
	int WAG; //Period replace lambda vectors...
       	vector<double> Z_ref; //Reference vector...
	vector<vector<double>> Lambda,DirectionVectors;
	vector<double> DCNLambda;
	vector<bool> IndividualChecked;
	typedef pair<double, int> Elt;
	vector<vector<int>> Neightborhood ;
        vector<Individual> Pool_P, EP;
        vector< vector<Individual> > Fronts;
	void InitUniformWeight();
	void InitNeighborhood();
	void InitIdealPoint();
        void InitializePool(vector<Individual> &Pool_P);
	void ComputeDirectionVectors();
	void Transformation_WS(vector<double> &W, vector<double> &L);
	void UpdateReference(Individual &ind);
	void UpdateNeighBors(Individual &ind, int IdProblem);
	void Reproduction(Individual &ind, int IndexProblem);
	void LocalSearch(Individual &ind, int k);
	void UpdateExternalPopulation();
	double getFitness(Individual &ind, vector<double> &SubLambda);
	double DistVector(vector<double> &A, vector<double> &B);
	double DistanceObjective(Individual &A, Individual &B);
	double DistanceVariable(Individual &A, Individual &B);
	void FastSort(vector<double> &X, vector<int> &Index, int T);
        bool CriterionStop();
        double getRandom(double Min, double Max);
        void End();
};
//vector< Individual > operator+(vector<Individual> PoolA, vector<Individual> PoolB );
//ostream & operator << ( ostream &out,vector<Individual> Pool);
//ostream & operator << ( ostream &out,vector<double> data);

#endif // MOEAD_HPP_INCLUDED
