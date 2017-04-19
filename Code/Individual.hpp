#ifndef INDIVIDUAL_HPP_INCLUDED
#define INDIVIDUAL_HPP_INCLUDED
#define BINARY_ENCODE 111
#define REAL_ENCODE 112
#include "Objectives.hpp"
using namespace std;
class Individual{
    public:
        Individual();
        void InitializeIndividual(Benchmark *ObjBenchmark);
	void TruncateVariable(int Dim, double Value);
        /**
            Se obtiene el objeto objectives el cual consta de la información
            de todas las funciones objetivo.
        **/

        /**
            Get
        **/
        inline Objectives& getObjectives(){ return ObjObjectives; }
        inline int getNumberObjectives(){return this->NumberObjectives;}
        inline int getDimension(){return this->Dimension;}
      	inline int getVariableRepresentation(){return this->VariableRepresentation;}
      	double getEvaluationMethod();
      	inline double getObjectiveValue(int N){ return this->ObjObjectives.SpaceObjectives[N].Fitness;}
      	inline double getFenotype(int Dimension){ return this->ObjObjectives.DecisionVariables[Dimension];}
	inline vector<double> getAllFenotype(){return this->ObjObjectives.DecisionVariables;}
      	inline bool Dominate(Individual &Ind){ return this->ObjObjectives.Dominate(Ind.ObjObjectives);}
      	inline int getNBits(){ return this->nbits;}
      	inline double getVariable(int d){ return this->ObjObjectives.DecisionVariables[d];}
      	inline double getMaximum(int d){ return this->Bounds[d][1]; }
      	inline double getMinimum(int d){ return this->Bounds[d][0];}
        /**
            Set
        **/
	inline void setVariableRepresentation(int VariableRepresentation){this->VariableRepresentation = VariableRepresentation;}
	inline void setFenotype(int Dimension, double Value){ this->ObjObjectives.DecisionVariables[Dimension] = Value; }
        void EvalIndividual();
        void DecodeIndividual(vector<double> & Genotype);

        /**
		Variable de desición para el caso encod-binary
	**/
	vector<vector<bool>> DecisionVariables;
    private:
        Objectives ObjObjectives;
        int Dimension, NumberObjectives, nbits;
	int VariableRepresentation;


        double getBase10(vector<bool> &Individuo, double Min , double Max);
        int BinarytoInt(vector<bool> &Individuo);
        vector<vector<double>> Bounds;
        void setRandomBinary();
	double DoubleRandom(double Min, double Max);
        //Generar boleanos aleatorios con una distribución bernoulli
        vector<bool>random_bool( double p  = 0.5);

};


#endif // INDIVIDUAL_HPP_INCLUDED
