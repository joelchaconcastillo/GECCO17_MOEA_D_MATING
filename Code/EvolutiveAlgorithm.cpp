/**
	Author: Joel Chacón Castillo
**/
#include <fstream>
#include <iostream>
#include <iomanip>
#include "Benchmark.hpp"
#include "EvolutiveAlgorithm.hpp"
#include "EAOperators.hpp"
EvolutiveAlgorithm::EvolutiveAlgorithm()
{
	srand(time(NULL));
}

void EvolutiveAlgorithm::Reproduction()
{

}
void EvolutiveAlgorithm::Generation()
{
	Selection();
	Recombination();
	Mutation();
	Improvement();
	Evaluation();
}
void EvolutiveAlgorithm::Init()
{
	InitializePopulation();

}
void EvolutiveAlgorithm::Run()
{
	Init();
	for(CurrentGeneration = 0; CurrentGeneration  < this->NumberGenerations; CurrentGeneration++)
	{
		Generation();
	}
	End();
}
/**
	Description: Inicializar a la población.
**/
void EvolutiveAlgorithm::InitializePopulation()
{
	this->Population.resize(this->SizePopulation);
	this->OffSpring.resize(this->SizePopulation);
	for(int i = 0; i < this->SizePopulation; i++)
	{
		this->Population[i].setVariableRepresentation(BINARY_ENCODE);
		this->Population[i].InitializeIndividual(this->ObjBenchmark);
	}
}
void EvolutiveAlgorithm::Selection()
{
	TournamentSelection(this->Population, this->OffSpring);
	this->Population = this->OffSpring;
}
void EvolutiveAlgorithm::Recombination()
{
	UniformCross( this->Population ,this->OffSpring, this->ProbCross);
	this->Population = this->OffSpring;
}
void EvolutiveAlgorithm::Mutation()
{
	BitMutation(this->Population, this->NBitsMut);
}
void EvolutiveAlgorithm::Evaluation()
{
        for(int i = 0; i < this->SizePopulation; i++)
        {
                this->Population[i].EvalIndividual();
        }
}
/**
	Se implementan heurísticas y/o búsquedas locales
**/
void EvolutiveAlgorithm::Improvement()
{

}

void EvolutiveAlgorithm::End()
{

}
void EvolutiveAlgorithm::ExportIndividualsFile(vector<Individual> &Pool, string Genotype, string Fenotype, string Fitness, int CuGeneration)
{
    ofstream SummaryGen, SummaryFen, SummaryFit;
	if(VariableRepresentation == BINARY_ENCODE)
    SummaryGen.open (Genotype, std::ofstream::out | std::ofstream::app );
    SummaryFen.open (Fenotype, std::ofstream::out | std::ofstream::app);
    SummaryFit.open (Fitness, std::ofstream::out | std::ofstream::app);

//    SummaryFit << "Size_Pool\t" <<Pool.size() <<" \tCurrent_Generation\t" << CuGeneration << "\t"  << endl;
//		if(VariableRepresentation == BINARY_ENCODE)
//    SummaryGen << "Size_Pool\t"<< Pool.size() << "\tCurrent_Generation\t" << CuGeneration << "\t"  << endl;
//    SummaryFen << "Size_Pool\t" <<Pool.size() << "\tCurrent_Generation\t" << CuGeneration << "\t" << endl;
    for(int i =0 ; i < Pool.size(); i++)
    {
        for(int j = 0; j < Pool[i].getNumberObjectives(); j++)
        {
		SummaryFit<<fixed << setprecision(30) << Pool[i].getObjectiveValue(j)<< "\t";
        }
	for(int d = 0; d < Pool[i].getDimension() ; d++)
	{
		if(VariableRepresentation == BINARY_ENCODE)
		{
			for(int b = 0; b < Pool[i].getNBits(); b++)
			{
			   SummaryGen << Pool[i].DecisionVariables[d][b] << "\t" ;
			}
		}
		SummaryFen << fixed << setprecision(30) <<Pool[i].getFenotype(d) << "\t";
	}
        SummaryFit<<endl;
				if(VariableRepresentation == BINARY_ENCODE)
					SummaryGen<<endl;
	SummaryFen<<endl;
    }
	SummaryFit.close();
if(VariableRepresentation == BINARY_ENCODE)
	SummaryGen.close();
	SummaryFen.close();
}
/**
	Sobrecarga de los operadores..
**/

ostream & operator << ( ostream &out,vector<Individual> Pool)
{
    for(int i = 0; i < Pool.size(); i++)
    {
        for(int j = 0; j < Pool[i].getObjectives().getNumberObjectives() ; j++)
         out << Pool[i].getObjectives().SpaceObjectives[j].Fitness<<" ";
         out << endl;
    }
    return out;
}
ostream & operator << ( ostream &out,vector<double> data)
{
    for(int i = 0; i < data.size(); i++)
    {
        out << data[i] << endl;
    }
    return out;
}
ostream & operator << ( ostream &out,vector<vector<double>> data)
{
    for(int i = 0; i < data.size(); i++)
    {
        for(int j = 0; j < data[i].size() ; j++)
         out << data[i][j]<<" ";
         out << endl;
    }
    return out;
}
ostream & operator << ( ostream &out,vector<int> data)
{
    for(int i = 0; i < data.size(); i++)
    {
         out << data[i]<<" ";
    }
    return out;
}
vector<Individual> operator+(vector<Individual> PoolA, vector<Individual> PoolB )
{
    vector<Individual> C;
    for(int i=0 ; i < PoolA.size(); i++ )
        C.push_back(PoolA[i]);
    for(int i = 0; i < PoolB.size(); i++)
        C.push_back(PoolB[i]);
    return C;
}
