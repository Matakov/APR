#include<iostream>
#include<iomanip>
#include<sstream>
#include<functional>
#include<fstream>
#include<algorithm>
#include<vector>
#include<stdlib.h>
#include<math.h>
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

//Abstraktna klasa koja broji koliko je puta pozvana koja funkcija
class AbstractFunction 
{
	protected:
		int numberOfCalls = 0;
		AbstractFunction(){};
		void restartCounting(){numberOfCalls=0;}
		void increase(){numberOfCalls++;}
	public:
		virtual ~AbstractFunction(){}
		//virtual void info() = 0;
		int getNumbers(){return numberOfCalls;}
		virtual void restartCount()=0;
		virtual double function(std::vector<double> lista){};
		virtual double function(double a, double b){};
		virtual double function(double a){};
		
};

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

class function2: public AbstractFunction 
{
	public:
		function2():AbstractFunction(){}
		double function(double a, double b)
		{
			increase();
			return pow(a-4,2)+ 4*pow(b-2,2);
		}
		double function(std::vector<double> lista)
		{	
			if (lista.size()>2)
			{
				std::cout<<"Input vectors are not good!"<<std::endl;
				throw(std::invalid_argument( "Input vectors are not good!" ));
							
			}
			increase();
			return pow(lista[0]-4,2)+ 4*pow(lista[1]-2,2);
		}
		void restartCount()
		{
			restartCounting();
		}
};

class function3: public AbstractFunction 
{
	private:
		double b;
	public:
		function3():AbstractFunction(){}
		function3(double input):b(input),AbstractFunction(){}
		double function(std::vector<double> lista)
		{
			increase();
			int numDim = lista.size();
			double sum = 0;
			for(int i=0;i<numDim;i++)
			{
				sum +=pow((lista[i]-i),2);
			}
			return sum;
		}
		//jednodimenzijska funkcija sa pomakom
		double function(double a)
		{
			increase();
			return pow(a-b,2);
		}
		void restartCount()
		{
			restartCounting();
		}
};

//Jakobovićeva funkcija
class function4: public AbstractFunction 
{
	public:
		function4():AbstractFunction(){}
		double function(std::vector<double> lista)
		{	
			if (lista.size()>2)
			{
				std::cout<<"Input vectors are not good!"<<std::endl;
				throw(std::invalid_argument( "Input vectors are not good!" ));
							
			}
			increase();
			return abs((lista[0]-lista[1])*(lista[0]+lista[1]))+sqrt(pow(lista[0],2)+pow(lista[1],2));
		}
		double function(double a, double b)
		{	
			increase();
			return abs((a-b)*(a+b))+sqrt(pow(a,2)+pow(b,2));
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


/*
Postupak trazenja unimodalnog intervala

Ulazne velicine:
- tocka: pocetna tocka pretrazivanja
- h: pomak pretrazivanja
- f: ciljna funkcija

Izlazne vrijednosti:
- unimodalni interval [l, r]
*/
void unimodalni(double h, double tocka, double& l, double &r, AbstractFunction& Class)
{
	l = tocka - h, r = tocka + h; 
	double m = tocka;
	double fl, fm, fr;
	int step = 1;

	fm = Class.function(tocka);
	fl = Class.function(l);
	fr = Class.function(r);

	if(fm < fr && fm < fl)	return;
	else if(fm > fr)
		do
		{	l = m;
			m = r;
			fm = fr;
			r = tocka + h * (step *= 2);
			fr = Class.function(r);
		} while(fm > fr);
	else 
		do
		{	r = m;
			m = l;
			fm = fl;
			l = tocka - h * (step *= 2);
			fl = Class.function(l);
		} while(fm > fl);
}

/*
Algoritam zlatnog reza

ulazne velicine:
- a, b: pocetne granice unimodalnog intervala
- e: preciznost

- t: pocetna tocka
- h: korak

- point: da li je dan unimodalni ili pocetna tocka
*/
/*
double Zlatni_rez(bool point,double h, double t,double e, std::function<double (double)> f)
{	
	double a,b;
	if(point)
	{
		//izracunaj prvo unimodalni interval
		unimodalni(h,t,a,b,f);
	}
	else
	{
		a = h;
		b = t;
	}
	double k = 0.5*(sqrt(5)-1);
	double c = b - k * (b - a);
	double d = a + k * (b - a);
	double fc = f(c);
	double fd = f(d);
	while((b - a) > e)
	{
		if(fc < fd) {
			b = d;
			d = c;
			c = b - k * (b - a);
			fd = fc;
			fc = f(c);
		}
		else
		{
			a = c;
			c = d;
			d = a + k * (b - a);
			fc = fd;
			fd = f(d);
		}
	}
	return (a + b)/2; // ili nove vrijednosti a i b
}
*/

//Zlatni rez sa zadanom pocetnom tockom
/*
double Zlatni_rez(double h, double t,double e, std::function<double (double)> f)
{	
	double a,b;
	double k = 0.5*(sqrt(5)-1);
	//izracunaj prvo unimodalni interval
	unimodalni(h,t,a,b,f);
	double c = b - k * (b - a);
	double d = a + k * (b - a);
	double fc = f(c);
	double fd = f(d);
	while((b - a) > e)
	{
		if(fc < fd) {
			b = d;
			d = c;
			c = b - k * (b - a);
			fd = fc;
			fc = f(c);
		}
		else
		{
			a = c;
			c = d;
			d = a + k * (b - a);
			fc = fd;
			fd = f(d);
		}
	}
	return (a + b)/2; // ili nove vrijednosti a i b
}
*/


void myfunction (std::string i) {  // function:
  std::cout << ' ' << i <<std::endl;
}





int main(int argc, char* argv[]){
	/*
	std::ifstream myfile;
	myfile.open(argv[1]);
	std::string line;
	double preciznost;
	double tocka,a,b;
	//std::cout<<argc<<std::endl;
	if(argc<2)
	{
		std::cout<<"Nisi zadao ulazni txt file!!\nPrekidam program."<<std::endl;
		std::exit(0);
	}
	if(myfile.is_open())
	{
		while(getline (myfile,line))
		{
			std::vector<std::string> container = split(line,' ');
			if (container[0] == "Preciznost:")
			{
				preciznost=atof(container[1].c_str());			
			}
			if (container[0] == "Početna" && container[1] == "točka:")
			{
				tocka=atof(container[2].c_str());			
			}
			if (container[0] == "Početni" && container[1] == "interval:")
			{
				a=atof(container[2].c_str());
				b=atof(container[3].c_str());			
			}
			//std::cout<<container[2]<<std::endl;
			//std::for_each (container.begin(), container.end(), myfunction);
			
		}
	}
	std::cout<<preciznost<<" "<<tocka<<" "<<a<<" "<<b<<std::endl;
	//std::cout<<func3(tocka)<<std::endl;
	*/
	/*
	function1 func1;
	double ide = func1.function(1,2);
	std::cout<<ide<<" "<<func1.getNumbers()<<std::endl;
	ide = func1.function(1,3);
	std::cout<<ide<<" "<<func1.getNumbers()<<std::endl;
	func1.restartCount();
	std::cout<<func1.getNumbers()<<std::endl;

	function2 func2;
	ide = func2.function(1,2);
	std::cout<<ide<<" "<<func2.getNumbers()<<std::endl;
	ide = func2.function(1,3);
	std::cout<<ide<<" "<<func2.getNumbers()<<std::endl;
	func2.restartCount();
	std::cout<<func2.getNumbers()<<std::endl;

	std::vector<double>v;
	v.push_back(2);
	v.push_back(2);

	function3 func3;
	ide = func3.function(v);
	std::cout<<ide<<" "<<func3.getNumbers()<<std::endl;
	//v.push_back(3);
	ide = func3.function(v);
	std::cout<<ide<<" "<<func3.getNumbers()<<std::endl;
	func3.restartCount();
	std::cout<<func3.getNumbers()<<std::endl;

	function4 func4;
	ide = func4.function(v);
	std::cout<<ide<<" "<<func4.getNumbers()<<std::endl;
	ide = func4.function(v);
	std::cout<<ide<<" "<<func4.getNumbers()<<std::endl;
	std::cout<<func4.getNumbers()<<std::endl;
	ide = func4.function(2,2);
	std::cout<<ide<<" "<<func4.getNumbers()<<std::endl;
	ide = func4.function(2,2);
	std::cout<<ide<<" "<<func4.getNumbers()<<std::endl;
	func4.restartCount();
	std::cout<<func4.getNumbers()<<std::endl;

	function6 func6;
	ide = func6.function(v);
	std::cout<<ide<<" "<<func6.getNumbers()<<std::endl;
	v.push_back(3);
	ide = func6.function(v);
	std::cout<<ide<<" "<<func6.getNumbers()<<std::endl;
	func6.restartCount();
	std::cout<<func6.getNumbers()<<std::endl;
	*/

	std::vector<double>v;
	v.push_back(2);
	v.push_back(2);

	function3 func3(3);
	double ide = func3.function(3);
	std::cout<<ide<<" "<<func3.getNumbers()<<std::endl;
	
	double i,j,h=1;
	unimodalni(h,10,i,j,func3);
	std::cout<<i<<" "<<j<<std::endl;
	//std::cout<<Zlatni_rez(h,tocka,preciznost,func3)<<std::endl;
	return 0;
}
