/**
	Author: Joel Chacón Castillo
**/
#ifndef EAOperators_HPP_INCLUDED
#define EAOperators_HPP_INCLUDED
#include "Individual.hpp"
using namespace std;
//class EAOperators
//{
//	public:
//		EAOperators();
//	protected:

	double NextRand();
	/**
		Operadores relacionados con la selección
	**/
	void TournamentSelection(vector<Individual> &Population, vector<Individual> &OffSpring);
	void RouletteWhell();
	void SUS();//Stochastic Universal Sampling
	/**
		Operadores relacionados con la cruza..
	**/
	//Binary Representation
	void OnePointCrossOver();
	void UniformCross(vector<Individual> &Population ,vector<Individual> &OffSpring, double ProbCross);
	void UniformCrossIndividual(Individual &A, Individual &B, Individual &OffSpring, double ProbCross);
	void NPointCrossOver();
	//Real Representation
	void LinearCrossover();
	void HBX(vector<Individual> &Population, vector<Individual> &OffSpring, double PCross);
	void HBXIndividual(Individual &Parent1, Individual &Parent2, Individual &Parent3, Individual &OffSpring1, double PCross );
	void FuzzyCrossOver(vector<Individual> &Population, vector<Individual> &OffSpring, double PCross);
	void FuzzyIndividual(Individual &Parent1, Individual &Parent2, Individual &OffSpring1, Individual &OffSpring2, double PCross );
	void BLX(vector<Individual> &Parent, vector<Individual> &OffSpring, double PCross); //Blend Crossover
	void BLXIndividual(Individual &Parent1, Individual &Parent2, Individual &OffSpring1, Individual &OffSpring2, double PCross );
	void SBXHybrid(vector<Individual> &Population ,vector<Individual> &OffSpring, double PCross);
	void SBXIndividualHybrid(Individual &Parent1, Individual &Parent2, Individual &OffSpring1, Individual &OffSpring2, double PCross );
	void SBX(vector<Individual> &Population ,vector<Individual> &OffSpring, double PCross); //Simulated Binary Crossover
	void SBXIndividual(Individual &Parent1, Individual &Parent2, Individual &Offspring1, Individual &Offspring2, double PCross );
	 void UNDX(); //Unimodal Normally Distributed Crossover
	  void SPX(); //Simplex Crossover
	 void UnfairAverageCrossover();
	
	/**
		Operadores relacionados con la mutación..
	**/
	//Binary Representation
	 void BitMutation(vector<Individual> &Population, int NBitsMut);
	 void BitMutationIndividual(Individual &A, int NBitsMut);
	//Real Representation
	 void RandomMutation(vector<Individual> &Population);
	 void NonUniformMutation();
	 void NormallyDistribuitedMutation(vector<Individual> &Population);
	 void PolynomialMutation(vector<Individual> &Population, double p_mut);
	 void PolynomialMutationIndividual(Individual &Ind, double eta, double p_mut);
//};

#endif
