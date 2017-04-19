/**
	Author: Joel Chacón Castillo
**/
#ifndef EVOLUTIVEALGORITHM_HPP_INCLUDED
#define EVOLUTIVEALGORITHM_HPP_INCLUDED
#include <string>
#include "Individual.hpp"
#include "Benchmark.hpp"
#include "EAOperators.hpp"	
using namespace std;
class EvolutiveAlgorithm 
{
	public: 
		EvolutiveAlgorithm();	
		virtual void Run();
		void ExportIndividualsFile(vector<Individual> &Pool, string Genotype, string Fenotype, string Fitness, int CurrenGeneration);
	protected:
		int SizePopulation, NumberGenerations, NBitsMut, CurrentGeneration=0, VariableRepresentation = REAL_ENCODE;
		string FilenameGenotype, FilenameFenotype, FilenameFitness;
		double ProbCross, ProbMutation; 
		Benchmark *ObjBenchmark;
		vector<Individual> Population, OffSpring;

		virtual void Init();	
		virtual void Reproduction();
		virtual void Generation();
		virtual void InitializePopulation();
		virtual void Selection();
		virtual void Recombination();
		virtual void Mutation();
		virtual void Evaluation();
		virtual void Improvement();
		virtual void End();
		};
vector< Individual > operator+(vector<Individual> PoolA, vector<Individual> PoolB );
ostream & operator << ( ostream &out,vector<Individual> Pool);
ostream & operator << ( ostream &out,vector<vector<double>> Pool);
ostream & operator << ( ostream &out,vector<double> data);
ostream & operator << ( ostream &out,vector<int> data);


#endif
