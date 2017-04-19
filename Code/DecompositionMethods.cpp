#include <vector>
#include <cmath>
#include "DecompositionMethods.hpp"
double ScalarTchebycheffApproach(vector<double> &Weight, vector<double> &IdealVector ,Individual &Ind)
{
	double Max_TValue = -1e30;
	for(int m = 0; m < IdealVector.size(); m++)
	{
		double Difference = fabs( Ind.getObjectiveValue(m) - IdealVector[m] );
		double Eval =0;
		if(Weight[m] == 0)
			Eval = 0.000001*Difference;
		else
			Eval = Weight[m]*Difference;
		if(Eval > Max_TValue) Max_TValue = Eval;
	}
	return Max_TValue;
}
