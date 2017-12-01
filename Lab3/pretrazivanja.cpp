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
	public:
		function3():AbstractFunction(){}
		double function(std::vector<double> lista)
		{	
			if (lista.size()>2)
			{
				std::cout<<"Input vectors are not good!"<<std::endl;
				throw(std::invalid_argument( "Input vectors are not good!" ));
							
			}
			increase();
			//return (lista[0]-2)*(lista[0]+lista[1]))+sqrt(pow(lista[0],2)+pow(lista[1],2));
			return pow(lista[0]-2,2)+pow(lista[1]+3,2);
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


class function4: public AbstractFunction 
{
	public:
		function4():AbstractFunction(){}
		double function(std::vector<double> lista)
		{
			increase();
			int numDim = lista.size();
			double sum = 0;
			for(int i=0;i<numDim;i++)
			{
				sum +=pow((lista[i]),2);
			}
			return pow(lista[0]-3,2)+pow(lista[1],2);
		}
		void restartCount()
		{
			restartCounting();
		}
};

//function for subtracting vectors
void subtractSame(std::vector<double>& a,std::vector<double> b)
{
	for(int i=0;i<a.size();i++)
	{
		a[i]-=b[i];	
	}
}
void absoluteValue(std::vector<double>& a)
{
	//for(int i=0;i<a.size();i++) std::cout<<a[i]<<" ";
	for(int i=0;i<a.size();i++)
	{
		a[i]=abs(a[i]);	
	}
	//for(int i=0;i<a.size();i++) std::cout<<a[i]<<" ";
	return;
}

std::vector<double> subtract(std::vector<double> a,std::vector<double> b)
{
	std::vector<double> c;
	//std::cout<<"Vektori :";
	//for(int i=0;i<a.size();i++) std::cout<<a[i]<<" ";
	//for(int i=0;i<b.size();i++) std::cout<<b[i]<<" ";
	for(int i=0;i<a.size();i++)
	{
		c.push_back(a[i]-b[i]);	
	}
	//for(int i=0;i<c.size();i++) std::cout<<c[i]<<" ";
	//std::cout<<std::endl;
	return c;
}

//function for multiplying vectors with a constant
auto multiply =[](std::vector<double> container,double cons) -> std::vector<double> {for (int i=0 ; i<container.size();i++) {container[i]=cons*container[i];}return container;};

//function for dividing vectors with a constant
auto divide =[](std::vector<double> container,double cons) -> std::vector<double> {for (int i=0 ; i<container.size();i++) {container[i]=container[i]/cons;}return container;};

//function for summation of vectors
void addSame(std::vector<double>& a,std::vector<double> b)
{
	for(int i=0;i<a.size();i++)
	{
		a[i]+=b[i];	
	}
}

std::vector<double> add(std::vector<double> a,std::vector<double> b)
{
	std::vector<double> c;
	for(int i=0;i<a.size();i++)
	{
		c[i]=a[i]+b[i];	
	}
}


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

//UNIMODALNI ZA VISE DIMENZIJA
void unimodalniMulti(double h, std::vector<double> tocka, double& l, double &r, AbstractFunction& Class, std::vector<double> multiToOne)
{
	double lambda0 = tocka[0]/multiToOne[0];
	//std::cout<<"Ulazni multi-vektor: ";
	//for(int i=0;i<tocka.size();i++) std::cout<<tocka[i]<<" "<<std::endl;
	std::vector<double> vl(tocka.size(),0.0);
	std::vector<double> vr(tocka.size(),0.0);
	lambdaL = lambda0 - h;
	lambdaR = lambda0 + h;
	for(int i=0;i<tocka.size();i++) vl[i] = multiToOne[i]*lambdaL; 
	for(int i=0;i<tocka.size();i++) vr[i] = multiToOne[i]*lambdaR; 
	std::vector<double> m = tocka;
	//std::cout<<"multi-vektor vl: ";
	//for(int i=0;i<tocka.size();i++) std::cout<<vl[i]<<" "<<std::endl;
	//std::cout<<"multi-vektor vr: ";
	//for(int i=0;i<tocka.size();i++) std::cout<<vr[i]<<" "<<std::endl;
	//std::cout<<"multi-vektor m: ";
	//for(int i=0;i<tocka.size();i++) std::cout<<m[i]<<" "<<std::endl;
	double fl, fm, fr;
	int step = 1;

	fm = Class.function(tocka);
	fl = Class.function(vl);
	fr = Class.function(vr);
	
	//std::cout<<"fl"<<" "<<"fm"<<" "<<"fr"<<std::endl;
	//std::cout<<fl<<" "<<fm<<" "<<fr<<std::endl;
	if(fm < fr && fm < fl)
	{
		l = vl[0]/multiToOne[0];
		r = vr[0]/multiToOne[0];	
		return;
	}
	else if(fm > fr)
	{
		do
		{	vl = m;
			m = vr;
			fm = fr;
			//vr[dim] = tocka[dim] + h * (step *= 2);
			for(int i=0;i<tocka.size();i++) vr[i] = multiToOne[i]*(lambda0 + h * (step *= 2));
			fr = Class.function(vr);
		} while(fm > fr);
		l = vl[0]/multiToOne[0];
		r = vr[0]/multiToOne[0];
	}
	else
	{
		do
		{	vr = m;
			m = vl;
			fm = fl;
			//vl[dim] = tocka[dim] - h * (step *= 2);
			for(int i=0;i<tocka.size();i++) vl[i] = multiToOne[i]*(lambda0 - h * (step *= 2));
			fl = Class.function(vl);
			//std::cout<<"fl"<<" "<<"fm"<<" "<<"fr"<<std::endl;
			//std::cout<<fl<<" "<<fm<<" "<<fr<<std::endl;
		} while(fm > fl);
		l = vl[0]/multiToOne[0];
		r = vr[0]/multiToOne[0];
	}
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

//ZLATNI REZ ZA VISE DIMENZIJA
double Zlatni_rezMulti(double h, std::vector<double> t,double e, AbstractFunction& Class, std::vector<double> multiToOne)
{	
	double lambda0 = tocka[0]/multiToOne[0];
	double a,b;
	//izracunaj prvo unimodalni interval
	unimodalniMulti(h,t,a,b,Class,multiToOne);
	//std::cout<<a<<" "<<b<<std::endl;
	
	double k = 0.5*(sqrt(5)-1);
	std::vector<double> vc(t.size(),0.0);
	double c = b - k * (b - a);
	//vc[dim]=c;
	for (int i=0;i<tocka.size();i++) vc[i] = multiToOne[i]*c;
	std::vector<double> vd(t.size(),0.0);
	double d = a + k * (b - a);
	//vd[dim]=d;
	for (int i=0;i<tocka.size();i++) vd[i] = multiToOne[i]*d;
	double fc = Class.function(vc);
	double fd = Class.function(vd);
	while((b - a) > e)
	{
		if(fc < fd) {
			b = vd[0]/multiToOne[0];
			vd = vc;
			//vc[dim] = b - k * (b - a);
			for (int i=0;i<tocka.size();i++) vc[i] = multiToOne[i]*(b - k * (b - a));
			fd = fc;
			fc = Class.function(vc);
		}
		else
		{
			a = vc[0]/multiToOne[0];
			vc = vd;
			//vd[dim] = a + k * (b - a);
			for (int i=0;i<tocka.size();i++) vd[i] = multiToOne[i]*(a + k * (b - a));
			fc = fd;
			fd = Class.function(vd);
		}
	}
	return (a + b)/2; // ili nove vrijednosti a i b
}

//function to numerically derive function
double derive(AbstractFunction& Class, std::vector<double> x0, double delta = 1.0e-6)
{
	std::vector<double> x1; //= x0 - delta;
	std::vector<double> x2; //= x0 + delta;

	for(int i=0;i<x0.size();i++) x1[i]=x0[i] - delta;	
	for(int i=0;i<x0.size();i++) x2[i]=x0[i] + delta;
	
	double y1 = Class.function(x1);
	double y2 = Class.function(x2);
	return (y2 - y1) / (x2 - x1);
}

//function to numerically partially derive function
double derivePartially(AbstractFunction& Class, std::vector<double> x0, double delta = 1.0e-6, int coord)
{
	std::vector<double> x1; //= x0 - delta;
	std::vector<double> x2; //= x0 + delta;

	for(int i=0;i<x0.size();i++){
	if (i==coord)
		{
		x1[i]=x0[i] - delta;
		}
	else
		{
		x1[i]=0;
		}
	}	
	for(int i=0;i<x0.size();i++){
	if (i==coord)
		{
		x2[i]=x0[i] + delta;
		}
	else
		{
		x2[i]=0;
		}
	}
	
	double y1 = Class.function(x1);
	double y2 = Class.function(x2);
	return (y2 - y1) / (x2 - x1);
}

void gradientDescent(AbstractFunction& Class, std::vector<double> x0, std::vector<double>& result, double delta = 1.0e-6, int mode=1)
{
	std::vector<double> partialDerivations,x_old,x_new;
	double lambda;
	sumDerivations = 0;
	x_old = x0;
	int iter=0;
	do {
		for (int i=0;i<x0.size();i++)
		{
			partialDerivations[i]=derivePartially(Class, x_old, delta, i);
			sumDerivations+=partialDerivations[i];
		}

		//normalize vectors
		for (int i=0;i<x0.size();i++)
		{
			partialDerivations[i]/=sumDerivations;
		}
		if(mode==1)
		{
			lambda = Zlatni_rezMulti(1,x_old,delta,Class,partialDerivations);
		}
		for(int i=0;i<x0.size();i++) x_new[i]=x_old[i]+lambda*partialDerivations;
		if(Class.function(x_new)>Class.function(x_old))
		{
			iter++;
			if(iter >= 100)
			{
				std::cout<<"It is divergating"<<std::endl;
				return;
			}
		}
		if(abs(Class.function(x_old)-Class.function(x_new))<delta)
		{
			result=x_new;
			return;
		}
		x_old=x_new;
	}
	while(true);
}


void openFile(std::string name,std::vector<double>& tocka,std::vector<double>& minimumFunkcije,std::vector<double>& preciznost, std::vector<double>& pomaciFunkcije,double& leftPoint,double& rightPoint,double& distance)
{
	std::ifstream myfile;
	std::string line;
	myfile.open(name);
	if(myfile.is_open())
	{
		while(getline (myfile,line))
		{
			std::vector<std::string> container = split(line,' ');
			if (container[0] == "Preciznost:")
			{
				for(int i=1;i<container.size();i++) preciznost.push_back(atof(container[i].c_str()));
				//preciznost=atof(container[1].c_str());			
			}
			if (container[0] == "Početna" && container[1] == "točka:")
			{
				for(int i=2;i<container.size();i++) tocka.push_back(atof(container[i].c_str()));			
			}
			if (container[0] == "Minimum" && container[1] == "funkcije:")
			{
				for(int i=2;i<container.size();i++) minimumFunkcije.push_back(atof(container[i].c_str()));			
			}
			if (container[0] == "Pomaci" && container[1] == "funkcije:")
			{
				for(int i=2;i<container.size();i++) pomaciFunkcije.push_back(atof(container[i].c_str()));			
			}
			if (container[0] == "Početni" && container[1] == "interval:")
			{
				leftPoint=atof(container[2].c_str());
				rightPoint=atof(container[3].c_str());			
			}
			if (container[0] == "Razmak" && container[1] == "točaka:")
			{
				distance=atof(container[2].c_str());		
			}
			//std::cout<<container[2]<<std::endl;
			//std::for_each (container.begin(), container.end(), myfunction);
			
		}
	}
}
int main(int argc, char* argv[]){
	std::ifstream myfile;
	std::string line;
	std::vector<double> tocka,minimumFunkcije,preciznost,pomaciFunkcije;
	double leftPoint,rightPoint,distance;
	int zadatak;
	//std::cout<<argc<<std::endl;
	if(argc<3)
	{
		std::cout<<"Nisi zadao ulazni txt file!!\nPrekidam program."<<std::endl;
		std::exit(0);
	}
	zadatak=atof(argv[2]);
	openFile(argv[1],tocka,minimumFunkcije,preciznost,pomaciFunkcije,leftPoint,rightPoint,distance);
	
	//for(int i=0;i<tocka.size();i++) std::cout<<preciznost[i]<<" ";
	//for(int i=0;i<tocka.size();i++) std::cout<<tocka[i]<<" ";
	//for(int i=0;i<minimumFunkcije.size();i++) std::cout<<minimumFunkcije[i]<<" ";
	//std::cout<<leftPoint<<" "<<rightPoint<<std::endl;
	//std::cout<<func3(tocka)<<std::endl;
	
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
	
	double alfa=1;
	double beta=0.5;
	double gamma=2;	
	
	if(zadatak==1)
	{
		double i,j,h=1;
		function3 func3(minimumFunkcije);
		std::vector<double> temp_result;
		std::cout<<"Zadatak 1"<<std::endl;
		//Zlatni rez
		std::cout<<"Minimum funkcije je: "<<Zlatni_rezMulti(true,h,tocka,preciznost,func3,0)<<std::endl;
		std::cout<<"Broj poziva funkcije u Zlatnom rezu je: "<<func3.getNumbers()<<std::endl;
		//resetiraj brojac	
		func3.restartCount();
		//Koordinatne osi	
		temp_result=KoordiantneOsi(tocka,preciznost,func3);
		std::cout<<"Minimum funkcije je: ";
		for(int i=0;i<temp_result.size();i++) std::cout<<temp_result[i]<<" ";
		std::cout<<std::endl;
		std::cout<<"Broj poziva funkcije u Koordinanim osima je: "<<func3.getNumbers()<<std::endl;
		//resetiraj brojac	
		func3.restartCount();
	
		
		//Nelder Mead
		std::cout<<"Nelder-Mead: "<<std::endl;	
		temp_result=NelderMead(tocka,distance,alfa,beta,gamma,preciznost[0],func3);
		std::cout<<"Minimum funkcije je: ";
		for(int i=0;i<temp_result.size();i++) std::cout<<temp_result[i]<<" ";
		std::cout<<std::endl;
		std::cout<<"Broj poziva funkcije u Nelder-Meadu je: "<<func3.getNumbers()<<std::endl;
		//resetiraj brojac
		func3.restartCount();
	
	
		//Hook Jeeves
		std::cout<<"Hook-Jeeves: "<<std::endl;	
		temp_result=HookeJeeves(tocka,preciznost,pomaciFunkcije,func3);
		std::cout<<"Minimum funkcije je: ";
		for(int i=0;i<temp_result.size();i++) std::cout<<temp_result[i]<<" ";
		std::cout<<std::endl;
		std::cout<<"Broj poziva funkcije u Hook-Jeevesu je: "<<func3.getNumbers()<<std::endl;
		//resetiraj brojac
		func3.restartCount();
	}
	if(zadatak==2)
	{
		std::vector<std::vector<double>> rezultati;
		std::vector<int> brojPoziva;
		
		function1 func1;
		function2 func2;
		function3 func3(minimumFunkcije);
		function4 func4;
		std::vector<double> tocka1,tocka2,tocka4,temp,tempRes;
		tocka1.push_back(-1.9);
		tocka1.push_back(2);
		tocka2.push_back(0.1);
		tocka2.push_back(0.3);
		tocka4.push_back(5.1);
		tocka4.push_back(1.1);
	
		//funkcija 1
		//Koordinatne osi
		tempRes=KoordiantneOsi(tocka1,preciznost,func1);
		brojPoziva.push_back(func1.getNumbers());
		
		std::cout<<"Koordinatne osi - 1"<<std::endl;
		std::cout<<"Minimum: ";
		for(int k=0;k<tempRes.size();k++) std::cout<<std::setw(5)<<tempRes[k]<<" ";
		std::cout<<"Broj poziva: "<<func1.getNumbers()<<std::endl;
		func1.restartCount();
	
		//Nelder-Mead
		tempRes=NelderMead(tocka1,distance,alfa,beta,gamma,preciznost[0],func1);
		brojPoziva.push_back(func1.getNumbers());
		std::cout<<"Nelder-Mead - 1"<<std::endl;
		std::cout<<"Minimum: ";
		for(int k=0;k<tempRes.size();k++) std::cout<<std::setw(5)<<tempRes[k]<<" ";
		std::cout<<"Broj poziva: "<<func1.getNumbers()<<std::endl;
		func1.restartCount();

		//Hook-Jeeves
		tempRes=HookeJeeves(tocka1,preciznost,pomaciFunkcije,func1);
		brojPoziva.push_back(func1.getNumbers());
		std::cout<<"Hook-Jeeves - 1"<<std::endl;
		std::cout<<"Minimum: ";
		for(int k=0;k<tempRes.size();k++) std::cout<<std::setw(5)<<tempRes[k]<<" ";
		std::cout<<"Broj poziva: "<<func1.getNumbers()<<std::endl;
		func1.restartCount();

		//funkcija 2
		//Koordinatne osi
		tempRes=KoordiantneOsi(tocka2,preciznost,func2);
		brojPoziva.push_back(func2.getNumbers());
		std::cout<<"Koordinatne osi - 2"<<std::endl;
		std::cout<<"Minimum: ";
		for(int k=0;k<tempRes.size();k++) std::cout<<std::setw(5)<<tempRes[k]<<" ";
		std::cout<<"Broj poziva: "<<func2.getNumbers()<<std::endl;
		func2.restartCount();

		//Nelder-Mead
		tempRes=NelderMead(tocka2,distance,alfa,beta,gamma,preciznost[0],func2);
		brojPoziva.push_back(func2.getNumbers());
		std::cout<<"Nelder-Mead - 2"<<std::endl;
		std::cout<<"Minimum: ";
		for(int k=0;k<tempRes.size();k++) std::cout<<std::setw(5)<<tempRes[k]<<" ";
		std::cout<<"Broj poziva: "<<func2.getNumbers()<<std::endl;
		func2.restartCount();

		//Hook-Jeeves
		tempRes=HookeJeeves(tocka2,preciznost,pomaciFunkcije,func2);
		brojPoziva.push_back(func2.getNumbers());
		std::cout<<"Hook-Jeeves - 2"<<std::endl;
		std::cout<<"Minimum: ";
		for(int k=0;k<tempRes.size();k++) std::cout<<std::setw(5)<<tempRes[k]<<" ";
		std::cout<<"Broj poziva: "<<func2.getNumbers()<<std::endl;
		func2.restartCount();
	
		//funkcija 3
		//Koordinatne osi
		tempRes=KoordiantneOsi(tocka,preciznost,func3);
		brojPoziva.push_back(func3.getNumbers());
		std::cout<<"Koordinatne osi - 3"<<std::endl;
		std::cout<<"Minimum: ";
		for(int k=0;k<tempRes.size();k++) std::cout<<std::setw(5)<<tempRes[k]<<" ";
		std::cout<<"Broj poziva: "<<func3.getNumbers()<<std::endl;
		func3.restartCount();

		//Nelder-Mead
		tempRes=NelderMead(tocka,distance,alfa,beta,gamma,preciznost[0],func3);
		brojPoziva.push_back(func3.getNumbers());
		std::cout<<"Nelder-Mead - 3"<<std::endl;
		std::cout<<"Minimum: ";
		for(int k=0;k<tempRes.size();k++) std::cout<<std::setw(5)<<tempRes[k]<<" ";
		std::cout<<"Broj poziva: "<<func3.getNumbers()<<std::endl;
		func3.restartCount();

		//Hook-Jeeves
		tempRes=HookeJeeves(tocka,preciznost,pomaciFunkcije,func3);
		brojPoziva.push_back(func3.getNumbers());
		std::cout<<"Hook-Jeeves - 3"<<std::endl;
		std::cout<<"Minimum: ";
		for(int k=0;k<tempRes.size();k++) std::cout<<std::setw(5)<<tempRes[k]<<" ";
		std::cout<<"Broj poziva: "<<func3.getNumbers()<<std::endl;
		func3.restartCount();

		//funkcija 4
		//Koordinatne osi
		tempRes=KoordiantneOsi(tocka4,preciznost,func4);
		brojPoziva.push_back(func4.getNumbers());
		std::cout<<"Koordinatne osi - 4"<<std::endl;
		std::cout<<"Minimum: ";
		for(int k=0;k<tempRes.size();k++) std::cout<<std::setw(5)<<tempRes[k]<<" ";
		std::cout<<"Broj poziva: "<<func4.getNumbers()<<std::endl;
		func4.restartCount();
	
		//Nelder-Mead
		tempRes=NelderMead(tocka4,distance,alfa,beta,gamma,preciznost[0],func4);
		brojPoziva.push_back(func4.getNumbers());
		std::cout<<"Nelder-Mead - 4"<<std::endl;
		std::cout<<"Minimum: ";
		for(int k=0;k<tempRes.size();k++) std::cout<<std::setw(5)<<tempRes[k]<<" ";
		std::cout<<"Broj poziva: "<<func4.getNumbers()<<std::endl;
		func4.restartCount();

		//Hook-Jeeves
		tempRes=HookeJeeves(tocka4,preciznost,pomaciFunkcije,func4);
		brojPoziva.push_back(func4.getNumbers());
		std::cout<<"Hook-Jeeves - 4"<<std::endl;
		std::cout<<"Minimum: ";
		for(int k=0;k<tempRes.size();k++) std::cout<<std::setw(5)<<tempRes[k]<<" ";
		std::cout<<"Broj poziva: "<<func4.getNumbers()<<std::endl;
		func4.restartCount();
	}
	if(zadatak==3)
	{
		function4 func4;
		std::vector<double> tocka4,nelderMead,hookJeeves;
		std::vector<int> brojPoziva;
		tocka4.push_back(5);
		tocka4.push_back(5);
		
		//Nelder-Mead
		nelderMead=NelderMead(tocka4,distance,alfa,beta,gamma,preciznost[0],func4);
		brojPoziva.push_back(func4.getNumbers());
		func4.restartCount();

		//Hook-Jeeves
		hookJeeves=HookeJeeves(tocka4,preciznost,pomaciFunkcije,func4);
		brojPoziva.push_back(func4.getNumbers());
		func4.restartCount();
		
		std::cout<<"Tablica rezultata"<<std::endl;
		std::cout<<"Minimum: "<<std::setw(30)<<"Broj poziva"<<std::endl;
		std::cout<<"Nelder-Mead"<<std::endl;
		for(int i=0;i<nelderMead.size();i++)std::cout<< nelderMead[i]<<" ";
		std::cout<<" "<<std::setw(20);
		std::cout<<brojPoziva[0];
		std::cout<<std::endl;
		std::cout<<"Hook-Jeeves"<<std::endl;
		for(int i=0;i<hookJeeves.size();i++)std::cout<< hookJeeves[i]<<" ";
		std::cout<<" "<<std::setw(20);
		std::cout<<brojPoziva[1];
		std::cout<<std::endl;
		
	}
	if(zadatak==4)
	{
		function1 func1;
		std::vector<double> tocka1,temp;
		std::vector<std::vector<double>> rezultati;
		std::vector<int> brojPoziva;
		tocka1.push_back(0.5);
		tocka1.push_back(0.5);
		for(int t=0;t<20;t++)
		{
			rezultati.push_back(NelderMead(tocka1,t,alfa,beta,gamma,preciznost[0],func1));
			brojPoziva.push_back(func1.getNumbers());
			func1.restartCount();
		}
		std::cout<<"Tablica rezultata"<<std::endl;
		std::cout<<"Nelder-Mead"<<std::endl;
		std::cout<<"Udaljnost točaka: "<<std::setw(10)<<"Broj poziva"<<std::setw(15)<<"Minimum"<<std::endl;
		for(int t=0;t<20;t++)
		{
			std::cout<<t+1<<std::setw(20);
			std::cout<<brojPoziva[t]<<std::setw(20);
			temp=rezultati[t];
			for(int i=0;i<temp.size();i++) std::cout<<temp[i]<<" ";
			std::cout<<std::endl; 			
		}
	}
	if(zadatak==5)
	{
		srand (time(NULL));
		function6 func6;
		std::vector<double> tocka6,temp;
		std::cout<<"Hook-Jeeves"<<std::endl;
		std::cout<<"Početne točake: "<<std::setw(15)<<"Broj poziva"<<std::setw(15)<<"Minimum"<<std::endl;
		for(int j=0;j<5;j++)
		{
			tocka6.clear();
			tocka6.push_back(rand()%101-50); //raspon od -50 do 50
			tocka6.push_back(rand()%101-50); //raspon od -50 do 50
			//temp = HookeJeeves(tocka6,preciznost,pomaciFunkcije,func6);
			temp = NelderMead(tocka6,10,alfa,beta,gamma,preciznost[0],func6);
			for(int i=0;i<tocka6.size();i++) std::cout<<tocka6[i]<<" ";
			std::cout<<" "<<std::setw(20);
			std::cout<<func6.getNumbers()<<std::setw(20);		
			for(int i=0;i<temp.size();i++) std::cout<<temp[i]<<" ";
			std::cout<<std::endl;
			func6.restartCount();
			
		}
		
	}
	return 0;

	
}
