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
		AbstractFunction(){};
		void restartCounting(){numberOfCalls=0;}
		void increase(){numberOfCalls++;}
		double r=0;
	public:
		virtual ~AbstractFunction(){}
		//virtual void info() = 0;
		int getNumbers(){return numberOfCalls;}
		virtual void restartCount()=0;
		virtual double function(std::vector<double> lista){};
		virtual double function(double a, double b){};
		virtual double function(double a){};
		virtual void setR(double r){};
		virtual double getR(){};
		
};

std::vector<std::string> split(const std::string& , char );

void openFile(std::string ,double& , double& ,double& ,double& ,double& ,double& ,double&, double& ,double&);

class function1: public AbstractFunction 
{
	public:
		function1();
		double function(double a, double b);
		double function(std::vector<double> lista);
		void restartCount();
};

class function3: public AbstractFunction 
{
	public:
		function3();
		double function(std::vector<double> lista);
		void restartCount();
};

class function6: public AbstractFunction 
{
	public:
		function6();
		double function(std::vector<double> lista);
		void restartCount();
};

class function7: public AbstractFunction 
{
	public:
		function7();
		double function(std::vector<double> lista);
		void restartCount();
};
