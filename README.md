This is the code from MOEA/D EVSD, the used language is C++,

Instructions:
*Inside the dirname "Code" execute "make", then a binary file will be generated as "MOEAD_EVSD"

Usage binary:
Instructions:
--Instance NAMEINSTANCE (WFG1)
--Sed 1 -->SED (299)
--Label OPERATOR_SBX --> is a postfix by filename
--Pcross 0.9 --> is the Crossover probability
--Pmut 0.3  --> is the Mutation Probability 
--Path ./RESULT --> is the directory where will save results
--SizePool 100 --> is the number of individual by generation
--Generations 25000  --> is the number of generations
--Crossover SBX --> is the Crossover operator to use
--Mutation Polynomial --> is the mutation operator to use
--PeriodReport 10 --> each N generations will save informaton in the "--Path" files

---------Crossover Operators Availables-------------
SBX, BLX, HBX, Fuzzy, SBXHybrid
---------Mutation Operators Availables--------------
Polynomial, NormalDistribution, Random

Example execution:

./MOEAD_EVSD --Instance WFG8 --Generations 10000 --PeriodReport 2500 --Path .
