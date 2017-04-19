#include <iostream>
#include <cstdlib>
#include <random>
#include <cmath>
#include "Benchmark.hpp"
#include "Objectives.hpp"
#include "Individual.hpp"
using namespace std;
Individual::Individual()
{
}
void Individual::InitializeIndividual(Benchmark *ObjBenchmark)
{

	this->VariableRepresentation = REAL_ENCODE;


        this->nbits = ObjBenchmark->getNBits();
        this->NumberObjectives = ObjBenchmark->getNObjectives();
        this->Dimension = ObjBenchmark->getDimension();

	ObjObjectives.setBenchmark(ObjBenchmark);
        Bounds = ObjBenchmark->Bounds;
       	if(this->VariableRepresentation == BINARY_ENCODE)
	{
	      	this->DecisionVariables.resize( this->Dimension, vector<bool>(nbits));
      		this->setRandomBinary();
        }
	else
	{
    //this->DecisionVariables.resize( this->Dimension);
		for(int i = 0; i < this->Dimension; i++)
		{
			//Obtener el valor en base a los límites..
			double Min = this->Bounds[i][0], Max = this->Bounds[i][1];
      			this->ObjObjectives.DecisionVariables[i] =  DoubleRandom(Min, Max);
    }
	}

	this->EvalIndividual();
}
void Individual::EvalIndividual()
{
	if(this->VariableRepresentation == BINARY_ENCODE)
	{
		this->DecodeIndividual(ObjObjectives.DecisionVariables);
	}
		ObjObjectives.Eval();
}
void Individual::DecodeIndividual(vector<double> &Genotype)
{
        for(int i = 0; i < this->Dimension; i++)
        {
		          Genotype[i] = (getBase10(this->DecisionVariables[i], Bounds[i][0], Bounds[i][1]));
	      }
}
vector<bool> Individual::random_bool( double p )
{
    vector<bool> bits(this->nbits);
    std::random_device rd;
    std::mt19937 gen( rd());
    std::bernoulli_distribution d( p);
    /**
        si el número de bits son muchos,
        entonces existe la posibilidad de un desboradmiento
        por el tipo de dato, por lo que sólo se randomiza
        una sección de cromosoma, como un centrómero, o
        los extremos adaptadores de 3 prima o 5 prima.
    **/
    for( int n = 0; n < bits.size(); ++n) {
        bits[bits.size()-1-n] = d(gen);
    }
    return bits;
}
void Individual::setRandomBinary()
{
        for(int j = 0; j < this->Dimension; j++)
        {
                DecisionVariables[j] = random_bool();
                /*for(int i = 0;  i < DecisionVariables[j].size(); i++)
                cout << DecisionVariables[j][i]<<" ";
                cout << endl;*/
        }
}
int Individual::BinarytoInt(vector<bool> &Individuo)
{
    int Integer=0;
    for(int i =0; i < Individuo.size(); i++)
    {
        if(Individuo[i])
            Integer+=pow(2, Individuo.size()-1-i);
    }
    return Integer;
}
double Individual::getBase10(vector<bool> &Individuo, double Min , double Max)
{
    double Delta = (double)(Max - Min) / ( pow(2,this->nbits) -1);
    double Xbase10 = (double) BinarytoInt(Individuo);
    Xbase10 =(Xbase10 *Delta) +Min;
    return Xbase10;
}
/**
	Description: Este método regresa el valor de aptitud, aptitud elevada es mejor por lo que se traducen las funciones objetivo MIN o MAX en una ponderación correspondiente a una aptitud del individuo.
**/
double Individual::getEvaluationMethod()
{
	//En caso de que sea una sola función objetivo
	if( this->NumberObjectives == 1 )
	{
		if( this->ObjObjectives.SpaceObjectives[0].Type == MAXIMIZE)
		{
			return this->ObjObjectives.SpaceObjectives[0].Fitness;
		}
		else
		{
			return -this->ObjObjectives.SpaceObjectives[0].Fitness;
		}
	}
}
void Individual::TruncateVariable(int Dim, double Value)
{

	double min = getMinimum(Dim), max =getMaximum(Dim);
	double psi = max - min;

	if(Value > max)
	{
		//double a = (Value-max)/psi;
		//this->setFenotype(Dim, max-(psi*(a-floor(a))));
		this->setFenotype(Dim, max);
	} else if(Value < min) {
		//double a = (min - Value)/psi;
		//this->setFenotype(Dim, min + psi*( a-floor(a))); 	
		this->setFenotype(Dim,min);
	} else {
		this->setFenotype(Dim, Value);
	}
}
double Individual::DoubleRandom(double Min, double Max)
{
	double f = (double) rand()/RAND_MAX;
	return Min + f*(Max - Min);
}
