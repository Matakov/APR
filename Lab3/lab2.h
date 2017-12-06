#include<iostream>
#include<iomanip>
#include<sstream>
#include<functional>
#include<fstream>
#include<algorithm>
#include<vector>
#include<stdlib.h>
#include<math.h>
#include <map>
#include <time.h>

class AbstractFunction 
{
	protected:
		int numberOfCalls = 0;
		AbstractFunction();
		void restartCounting();
		void increase();
	public:
		virtual ~AbstractFunction();
		//virtual void info() = 0;
		int getNumbers();
		virtual void restartCount()=0;
		virtual double function(std::vector<double> lista){};
		virtual double function(double a, double b){};
		virtual double function(double a){};
		
};

void unimodalni(double , double , double& , double& , AbstractFunction& );
void unimodalniMulti(double , std::vector<double> , double& , double& , AbstractFunction& , std::vector<double> );
double Zlatni_rezMulti(double , std::vector<double> ,double , AbstractFunction& , std::vector<double> );

std::vector<std::vector<double>> tockeSimpleksa(std::vector<double> ,double );
void getCentroid(std::vector<std::vector<double>> ,std::vector<double>& , double );
int getMaximumIndex(std::vector<std::vector<double>> , AbstractFunction& ,std::map<std::vector<double>, double>& );
std::vector<double> getMinimum(std::vector<std::vector<double>> , AbstractFunction& );
int getMinimumIndex(std::vector<std::vector<double>> , AbstractFunction& ,std::map<std::vector<double>, double>& );
void ekspanzija(std::vector<double> , std::vector<double> ,double ,double ,std::vector<double>& );
void kontrakcija(std::vector<double> , std::vector<double> ,double ,std::vector<double>& );
void pomakiPremaL(std::vector<std::vector<double>>& ,int );
bool kriterijZaustavljanja(std::vector<std::vector<double>> , std::vector<double> ,double eps,AbstractFunction& ,std::map<std::vector<double>, double>& );

std::vector<double> NelderMead(std::vector<double> ,double , double , double , double , double , AbstractFunction& );
