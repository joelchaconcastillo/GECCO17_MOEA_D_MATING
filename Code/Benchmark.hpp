#ifndef BENCHMARK_HPP_INCLUDED
#define BENCHMARK_HPP_INCLUDED
#include <vector>
#define MAXIMIZE 2
#define MINIMIZE 3
#define SCH1 10
#define SCH2 11
#define FON 12
#define KUR 13
#define POL 14
#define VNT 15
#define ZDT1 16
#define ZDT2 17
#define ZDT3 18
#define ZDT4 19
#define ZDT5 20
#define ZDT6 21
#define BNH 22
#define OSY 23
#define SRN 24
#define TNK 25
#define SINGLE 26
#define Type_WFG1 27
#define Type_WFG2 28
#define Type_WFG3 29
#define Type_WFG4 30
#define Type_WFG5 31
#define Type_WFG6 32
#define Type_WFG7 33
#define Type_WFG8 34
#define Type_WFG9 35

using namespace std;
class Benchmark{

    private:
            int NObjectives, TypeProblem, Dimension, DefaultDimension;
            vector<int> TypeDuality;
            void SetAllMinimize();
            void SetAllMaximize();
            void SetAllBounds(vector<vector<double>> &MatrixBound,double Min, double Max);
            int SizeChrom;
    public:
        vector<vector<double>> Bounds, OptimalBounds;
        Benchmark(int TypeProblem);
        Benchmark();
        void InitConfig();
        /**
            get
        **/
        inline int getNObjectives(){return this->NObjectives;}
        inline int getDimension(){return this->Dimension;}
        inline int getNBits(){return this->SizeChrom;}
        inline vector<int> getTypeDuality(){return this->TypeDuality;}

        void Eval(vector<double> &X, vector<double> &obj);
	void Single(vector<double> &X,vector<double> &obj);
        void sch1(vector<double> &X, vector<double> &obj);
        void sch2(vector<double> &X, vector<double> &obj);
        void fon(vector<double> &X, vector<double> &obj);
        void kur(vector<double> &X, vector<double> &obj);
        void pol(vector<double> &X, vector<double> &obj);
        void vnt(vector<double> &X, vector<double> &obj);
        void zdt1(vector<double> &X, vector<double> &obj);
        void zdt2(vector<double> &X, vector<double> &obj);
        void zdt3(vector<double> &X, vector<double> &obj);
        void zdt4(vector<double> &X, vector<double> &obj);
        void zdt5(vector<double> &X, vector<double> &obj);
        void zdt6(vector<double> &X, vector<double> &obj);
        void bnh(vector<double> &X, vector<double> &obj);
        void osy(vector<double> &X, vector<double> &obj);
        void srn(vector<double> &X, vector<double> &obj);
        void tnk(vector<double> &X, vector<double> &obj);
	void wfg1(vector<double> &X, vector<double> &obj);
	void wfg2(vector<double> &X, vector<double> &obj);
	void wfg3(vector<double> &X, vector<double> &obj);
	void wfg4(vector<double> &X, vector<double> &obj);
	void wfg5(vector<double> &X, vector<double> &obj);
	void wfg6(vector<double> &X, vector<double> &obj);
	void wfg7(vector<double> &X, vector<double> &obj);
	void wfg8(vector<double> &X, vector<double> &obj);
	void wfg9(vector<double> &X, vector<double> &obj);


};


#endif // BENCHMARK_HPP_INCLUDED
