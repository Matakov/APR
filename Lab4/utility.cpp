#include "utility.h"

//Rosenbrockova 'banana' funkcija
class function1: public AbstractFunction 
{
	public:
		function1():AbstractFunction(){}
		double function(double a, double b)
		{
			increase();
			return 100*pow(b-pow(a,2),2)+pow((1-a),2);
		}
		double function(std::vector<double> lista)
		{	
			if (lista.size()>2)
			{
				std::cout<<"Input vectors are not good!"<<std::endl;
				throw(std::invalid_argument( "Input vectors are not good!" ));
							
			}
			increase();
			return 100*pow(lista[1]-pow(lista[0],2),2)+pow((1-lista[0]),2);
		}
		void restartCount()
		{
			restartCounting();
		}
};

class function3: public AbstractFunction 
{
	public:
		function3():AbstractFunction(){}
		double function(std::vector<double> lista)
		{	
			double output=0;
			for(int i=0;i<lista.size();i++)
			{
				output+=(lista[i]-i)*(lista[i]-i);
			}
			increase();
			return output;
		}
		void restartCount()
		{
			restartCounting();
		}
};

//Schaffer's function
class function6: public AbstractFunction 
{
	public:
		function6():AbstractFunction(){}
		double function(std::vector<double> lista)
		{
			increase();
			int numDim = lista.size();
			double sum = 0;
			for(int i=0;i<numDim;i++)
			{
				sum +=pow((lista[i]),2);
			}
			return 0.5+(pow(sqrt(sum),2)-0.5)/pow((1+0.001*sum),2);
		}
		void restartCount()
		{
			restartCounting();
		}
};


//Almost Schaffer's function
class function7: public AbstractFunction 
{
	public:
		function7():AbstractFunction(){}
		double function(std::vector<double> lista)
		{
			increase();
			int numDim = lista.size();
			double sum = 0;
			for(int i=0;i<numDim;i++)
			{
				sum +=pow((lista[i]),2);
			}
			return pow(sum,0.25)*(1+pow(sin(50*pow(sum,2)),2));
		}
		void restartCount()
		{
			restartCounting();
		}
};

/*
String splitting function
*/
std::vector<std::string> split(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> tokens;
    while (getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}

void openFile(std::string name,double& brPopulacije, double& preciznost,double& vjerojatnost,double& brEval, double& borderLeft, double& borderRight, double& velicinaVektora)
{
	std::ifstream myfile;
	std::string line;
	myfile.open(name);
	if(myfile.is_open())
	{
		while(getline (myfile,line))
		{
			std::vector<std::string> container = split(line,' ');
			if (container[0] == "Veličina" && container[1] == "populacije:")
			{
				brPopulacije=atof(container[2].c_str());			
			}
			if (container[0] == "Broj" && container[1] == "bitova:")
			{
				preciznost=atof(container[2].c_str());			
			}
			if (container[0] == "Vjerojatnost" && container[1] == "mutacije:")
			{
				vjerojatnost=atof(container[2].c_str());
			}
			if (container[0] == "Broj" && container[1] == "evaluacija:")
			{
				brEval=atof(container[2].c_str());			
			}
			if (container[0] == "Početni" && container[1] == "interval:")
			{
				borderLeft=atof(container[2].c_str());
				borderRight=atof(container[3].c_str());			
			}
			if (container[0] == "Veličina" && container[1] == "vektora:")
			{
				velicinaVektora=atof(container[2].c_str());			
			}
			
		}
	}
}