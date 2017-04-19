#include <iostream>
#include <algorithm>
#include <random>
#include <cmath>
#include "Individual.hpp"
#include "EAOperators.hpp"
#define EPS 1.0e-14


//EAOperators::EAOperators()
//{
//
//}
double NextRand()
{
	return (double)rand()/RAND_MAX;
}
/**
        Operadores relacionados con la selección
**/
void TournamentSelection(vector<Individual> &Population, vector<Individual> &OffSpring)
{
	for(int i = 0; i < Population.size(); i++)
	{
		int Index1 = rand()%Population.size();	
		int Index2 = rand()%Population.size();
		//Se escoge al individuo con mayor aptitud..	
		if( Population[Index1].getEvaluationMethod() > Population[Index2].getEvaluationMethod()  )
		{
		 	OffSpring[i] = Population[Index1];		
		}
		else
		{
			OffSpring[i] = Population[Index2];
		}

	}
}
void RouletteWhell()
{

}
void SUS()
{

}
/**
        Operadores relacionados con la cruza..
**/
//Binary Representation
void OnePointCrossOver()
{

}
void UniformCross(vector<Individual> &Population ,vector<Individual> &OffSpring, double ProbCross)
{
	for(int i = 0; i < Population.size(); i++)
	{	
		int Index1 = rand()%Population.size();	
		int Index2 = rand()%Population.size();	
		UniformCrossIndividual(Population[Index1], Population[Index2], OffSpring[i], ProbCross );	
	}
}
void UniformCrossIndividual(Individual &A, Individual &B, Individual &OffSpring, double ProbCross)
{
		for(int D = 0; D < A.getDimension(); D++)
		{
			for(int bit = 0; bit < A.DecisionVariables[D].size() ; bit++)
			{
				double RandomNumber = rand()/(RAND_MAX);
				if(RandomNumber < ProbCross )
					OffSpring.DecisionVariables[D][bit] = A.DecisionVariables[D][bit]; 		
				else
					OffSpring.DecisionVariables[D][bit] = B.DecisionVariables[D][bit]; 		
			}
		}
}
void NPointCrossOver()
{

}
//Real Representation
void LinearCrossover()
{

}
/**
	HBX: Tomado del artículo "A Novel Adaptative Hybrnid Crossover Operator for Multiobjective Evolutionary" por Qingling Zhu
**/
void HBX(vector<Individual> &Population, vector<Individual> &OffSpring, double PCross)
{
    for(int i = 0; i < Population.size(); i++)
    {
	int x1 = rand()%Population.size();
	int x2 = rand()%Population.size();
	int x3 = rand()%Population.size();
        HBXIndividual(Population[x1], Population[x2], Population[x3], OffSpring[i], PCross);
    }
}
void HBXIndividual(Individual &Parent1, Individual &Parent2, Individual &Parent3, Individual &OffSpring1, double PCross )
{
    double eta_c = 10;
    double rand;
    double y1, y2, yl, yu;
    double c1, c2;
    double alpha, beta, betaq;
    double P1 = 0.5, P2 = 0.6;
    double F=0.5;
    if ( NextRand() <= PCross)
    {
        for (int i=0; i< Parent1.getDimension(); i++)
        {
            if (NextRand()<= P1 )
            {
                if (fabs(Parent1.getVariable(i)-Parent2.getVariable(i)) > EPS)
                {
                    if (Parent1.getVariable(i) < Parent2.getVariable(i))
                    {
                        y1 = Parent1.getVariable(i);
                        y2 = Parent2.getVariable(i);
                    }
                    else
                    {
                        y1 = Parent2.getVariable(i);
                        y2 = Parent1.getVariable(i);
                    }
		    yl = Parent1.getMinimum(i);
		    yu = Parent1.getMaximum(i);

                    rand = NextRand();
                    beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                    alpha = 2.0 - pow(beta,-(eta_c+1.0));
                    if (rand <= (1.0/alpha))
                    {
                        betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
                    }
                    else
                    {
                        betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                    }
                    c1 = 0.5*((y1+y2)-betaq*(y2-y1));
                    beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                    alpha = 2.0 - pow(beta,-(eta_c+1.0));
                    if (rand <= (1.0/alpha))
                    {
                        betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
                    }
                    else
                    {
                        betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                    }
                    c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                    if (c1<yl)
                        c1=yl;
                    if (c2<yl)
                        c2=yl;
                    if (c1>yu)
                        c1=yu;
                    if (c2>yu)
                        c2=yu;

		    double vi = Parent1.getVariable(i) + F*( Parent2.getVariable(i)-Parent3.getVariable(i));
		    if( NextRand() < P2)
		    {
			    if (NextRand()<=0.5)
			    {
				OffSpring1.TruncateVariable(i, c2);
			    }
			    else
			    {
				OffSpring1.TruncateVariable(i, c1);
			    }
		   }
		   else
				OffSpring1.TruncateVariable(i, vi);
                }
                else
                {
		    OffSpring1.TruncateVariable(i, Parent1.getVariable(i));
                }
            }
            else
            {
		OffSpring1.TruncateVariable(i, Parent2.getVariable(i));
            }
        }
    }
    else
    {
	OffSpring1 = Parent1;
    }
}
void BLX(vector<Individual> &Population, vector<Individual> &OffSpring, double PCross)
{
    for(int i = 0; i < Population.size(); i+=2)
    {
        BLXIndividual(Population[i], Population[i+1], OffSpring[i], OffSpring[i+1], PCross);
    }
}
void BLXIndividual(Individual &Parent1, Individual &Parent2, Individual &OffSpring1, Individual &OffSpring2, double PCross )
{
	if(NextRand() > PCross)
	{
		OffSpring1 = Parent1;
		OffSpring2 = Parent2;
		return;
	}
    double gamma, alpha=0.5;
    for(int d=0; d < Parent1.getDimension(); d++)
    {
        double Px1 = min(Parent1.getVariable(d),Parent2.getVariable(d));
        double Px2 = max(Parent1.getVariable(d),Parent2.getVariable(d));
        gamma = (1.0+2.0*alpha)*NextRand()-alpha;
        double xi_1 = (1.0-gamma)*Px1 + gamma*Px2;  
        OffSpring1.TruncateVariable(d, xi_1);
         
        gamma = (1.0+2.0*alpha)*NextRand()-alpha;
        double xi_2 = (1.0-gamma)*Parent1.getVariable(d) + gamma*Parent2.getVariable(d);  
        OffSpring2.TruncateVariable(d, xi_2);
    }
}
void FuzzyCrossOver(vector<Individual> &Population, vector<Individual> &OffSpring, double PCross)
{
	for(int i = 0; i < Population.size(); i+=2)
	{
		FuzzyIndividual(Population[i], Population[i+1], OffSpring[i], OffSpring[i+1], PCross);
	}
}
void FuzzyIndividual(Individual &Parent1, Individual &Parent2, Individual &OffSpring1, Individual &OffSpring2, double PCross )
{

	if(NextRand() > PCross)
	{
		OffSpring1 = Parent1;
		OffSpring2 = Parent2;
		return;
	}
	double D=0.5;
	for(int d=0; d < Parent1.getDimension(); d++)
	{
		double C1 = min(Parent1.getVariable(d),Parent2.getVariable(d)); //Modal Value parent1
		double C2 = max(Parent1.getVariable(d),Parent2.getVariable(d)); //Modal Value parent2
		register double Difference = fabs(C2-C1);
		double a1 =  C1-D*Difference;
		double b1 =  C1+D*Difference;
		double F1 = (C1-a1)/(b1-a1);
		double u1 = NextRand();
		double a2 =  C2-D*Difference;
		double b2 =  C2+D*Difference;
		double F2 = (C2-a2)/(b2-a2);
		double u2 = NextRand();

		double xi_1, xi_2;
		if(u1<=F1)
			xi_1 = a1+sqrt(u1*(b1-a1)*(C1-a1));
		else
			xi_1 = b1-sqrt((1.0-u1)*(b1-a1)*(b1-C1));

		if(u2<=F2)
			xi_2 = a2+sqrt(u2*(b2-a2)*(C2-a2));
		else
			xi_2 = b2-sqrt((1.0-u2)*(b2-a2)*(b2-C2));

		OffSpring1.TruncateVariable(d, xi_1);
		OffSpring2.TruncateVariable(d, xi_2);
	}
}
void SBXHybrid(vector<Individual> &Population ,vector<Individual> &OffSpring, double PCross) //Simulated Binary Crossover
{

	for(int i = 0; i < Population.size(); i+=2)
	{
			SBXIndividualHybrid(Population[i], Population[i+1], OffSpring[i], OffSpring[i+1], PCross);
	}
}
void SBXIndividualHybrid(Individual &Parent1, Individual &Parent2, Individual &OffSpring1, Individual &OffSpring2, double PCross )
{
//Indice de distribución...
    double eta_c = 10;
    double rand;
    double y1, y2, yl, yu;
    double c1, c2;
    double alpha, beta, betaq;
    if ( NextRand() <= PCross)
    {
        for (int i=0; i< Parent1.getDimension(); i++)
        {
            if (NextRand()<=0.5 )
            {
                if (fabs(Parent1.getVariable(i)-Parent2.getVariable(i)) > EPS)
                {
                    if (Parent1.getVariable(i) < Parent2.getVariable(i))
                    {
                        y1 = Parent1.getVariable(i);
                        y2 = Parent2.getVariable(i);
                    }
                    else
                    {
                        y1 = Parent2.getVariable(i);
                        y2 = Parent1.getVariable(i);
                    }
		    yl = Parent1.getMinimum(i);
		    yu = Parent1.getMaximum(i);

                    rand = NextRand();
                    beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                    alpha = 2.0 - pow(beta,-(eta_c+1.0));
                    if (rand <= (1.0/alpha))
                    {
                        betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
                    }
                    else
                    {
                        betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                    }
                    c1 = 0.5*((y1+y2)-betaq*(y2-y1));
                    beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                    alpha = 2.0 - pow(beta,-(eta_c+1.0));
                    if (rand <= (1.0/alpha))
                    {
                        betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
                    }
                    else
                    {
                        betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                    }
                    c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                    if (c1<yl)
                        c1=yl;
                    if (c2<yl)
                        c2=yl;
                    if (c1>yu)
                        c1=yu;
                    if (c2>yu)
                        c2=yu;
                    if (NextRand()<=0.5)
                    {
			OffSpring1.TruncateVariable(i, c2);
			OffSpring2.TruncateVariable(i, c1);
                    }
                    else
                    {
			OffSpring1.TruncateVariable(i, c1);
			OffSpring2.TruncateVariable(i, c2);
                    }
                }
                else
                {
		    OffSpring1.TruncateVariable(i, Parent1.getVariable(i));
		    OffSpring2.TruncateVariable(i, Parent2.getVariable(i));
                }
            }
            else
            {
		OffSpring1.TruncateVariable(i, Parent1.getVariable(i));
		OffSpring2.TruncateVariable(i, Parent2.getVariable(i));
            }
        }
    }
    else
    {
	OffSpring1 = Parent1;
	OffSpring2  = Parent2;
    }
    return;
}
void SBX(vector<Individual> &Population ,vector<Individual> &OffSpring, double PCross) //Simulated Binary Crossover
{

	for(int i = 0; i < Population.size(); i+=2)
	{
			SBXIndividual(Population[i], Population[i+1], OffSpring[i], OffSpring[i+1], PCross);
	}
}
/**
	Esta implementación considera los límites en la función de distribución...
**/	
void SBXIndividual(Individual &Parent1, Individual &Parent2, Individual &OffSpring1, Individual &OffSpring2, double PCross )
{
	double eta=10;
	double Rand;
        double y1, y2, yl, yu;
        double c1, c2;
        double alpha, beta, betaq;
		if(NextRand() > PCross)
		{
			OffSpring1 = Parent1;
			OffSpring2 = Parent2;
			return;
		}
		for(int d = 0; d < Parent1.getDimension(); d++)
		{
			if (fabs(Parent1.getVariable(d) - Parent2.getVariable(d)) > EPS) 
			{
			   //if( NextRand() <= 0.5)
			   //{
			    if ( Parent1.getVariable(d) < Parent2.getVariable(d))
			    {
				y1 = Parent1.getVariable(d);
				y2 = Parent2.getVariable(d);
			    }
			    else
			    {
				y1 = Parent2.getVariable(d);
				y2 = Parent1.getVariable(d);
			    }
			    yl = Parent1.getMinimum(d);
			    yu = Parent1.getMaximum(d);

			    Rand = NextRand();
			    beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1));
			    alpha = 2.0 - pow( beta, -(eta + 1.0));
			    if (Rand <= (1.0/alpha))
			    {
				betaq = pow ( (Rand * alpha), (1.0 / (eta + 1.0)));
			    }
			    else
			    {
				betaq = pow ( (1.0 / (2.0 - Rand * alpha)), (1.0 / (eta+1.0)));
			    }
			    c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));

			    beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
			    alpha = 2.0 - pow( beta, -(eta + 1.0));
			    if (Rand <= (1.0/alpha))
			    {
				betaq = pow ( (Rand * alpha), (1.0 / (eta + 1.0)));
			    }
			    else
			    {
				betaq = pow ( (1.0 / (2.0 - Rand * alpha)), (1.0 / (eta + 1.0)));
			    }
			    c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));
			    //if(NextRand() <= 0.5 )
			    {
			    	OffSpring1.TruncateVariable(d, c1);  
			    	OffSpring2.TruncateVariable(d, c2);  
			    }
			//     else
			//	{
			//		OffSpring1.TruncateVariable(d, c2);  
			//		OffSpring2.TruncateVariable(d, c1);  

			//	}
			//}
			//else
			//{
			//	OffSpring1.TruncateVariable(d, Parent1.getVariable(d));  
			//	OffSpring2.TruncateVariable(d, Parent2.getVariable(d));  


			//}
			}
	}
}
void UNDX()
{

}
void SPX()
{

}
void UnfairAverageCrossover()
{

}

/**
        Operadores relacionados con la mutación..
**/
//Binary Representation
void BitMutation(vector<Individual> &Population, int NBitsMut)
{
	for(int i = 0; i < Population.size(); i++)
	{
	
	}
}
void BitMutationIndividual(Individual &A, int NBitsMut)
{
		for(int D = 0; D < A.getDimension() ; D++)
		{
			int ChangedBits = 0;
                 	while( ChangedBits < NBitsMut)
                 	{
				int LocusBit = 	rand()% A.DecisionVariables[D].size(); 
			        if( A.DecisionVariables[D][LocusBit])
					A.DecisionVariables[D][LocusBit] = false;
				else
					A.DecisionVariables[D][LocusBit] = true;
			        
				ChangedBits++;
		 	}
		 	
		 }
	
}
//Real Representation
void RandomMutation(vector<Individual> &Population)
{
	double Delta = 0.4;
	for(int i = 0; i < Population.size(); i++)
	{
		for(int d = 0; d < Population[i].getDimension(); d++)
		{
			double Value = Population[i].getFenotype(d);
			double Rand = (double)rand()/RAND_MAX;
			Population[i].setFenotype(d, Value + (Rand -0.5)*Delta);

		}
	}
}
void NonUniformMutation()
{

}
void NormallyDistribuitedMutation(vector<Individual> &Population)
{
	random_device rd;
	mt19937 gen(rd());
	normal_distribution<double> Normal(0.0, 1.0);
	for(int i = 0; i < Population.size(); i++)
	{
		for(int d = 0; d < Population[i].getDimension(); d++)
		{
			double Value = Population[i].getFenotype(d);
			double Rand = Normal(gen);
			Population[i].setFenotype(d, Value + Rand);

		}
	}
}
void PolynomialMutation(vector<Individual> &Population, double p_mut)
{
	for(int i = 0; i < Population.size(); i++)
	{
		PolynomialMutationIndividual(Population[i], 50, p_mut);
	}
}


void PolynomialMutationIndividual(Individual &Ind, double eta, double p_mut )
{
	double rnd, delta1, delta2, mut_pow, deltaq, delta_max;
	double y, yl, yu, val, xy;
	mut_pow = 1.0/(eta+1.0);
	for(int j = 0; j < Ind.getDimension(); j++){
		if (((double)rand()/RAND_MAX) <= p_mut){
			y = Ind.getVariable(j);
			yl = Ind.getMinimum(j);
			yu = Ind.getMaximum(j);
			delta1 = (y-yl) / (yu-yl);
			delta2 = (yu-y) / (yu-yl);
			if( (y-yl) > (yu-y))
				delta_max = delta2;
			else
				delta_max = delta1;
			double rnd = (double)rand()/RAND_MAX;
			if(rnd <= 0.5){
				xy = 1.0-delta1;
				val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta+1.0)));
				deltaq =  pow(val,mut_pow) - 1.0;
			} else {
				xy = 1.0-delta2;
				val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta+1.0)));
				deltaq = 1.0 - (pow(val,mut_pow));
			}
			y = y + deltaq*(yu-yl); 
		
		        //////////////
			//double d = (double)rand()/RAND_MAX;
			//if( d < 0.3)
			Ind.TruncateVariable(j, y);	
			//else
			//{
			//	double v = (double)rand()/RAND_MAX;
			//	v = (yu - yl)*v;
			//	Ind.TruncateVariable(j, v/2.0);	
			//}
		}
	}
}

