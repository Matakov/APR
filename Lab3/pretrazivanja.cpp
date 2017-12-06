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

#include "lab1.h"

//----------------------------- LABOS 1 MATRICE


// KRAJ KLASE



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

//Here are implitic functions
class function5: public AbstractFunction 
{
	public:
		function5():AbstractFunction(){}
		double function(std::vector<double> lista)
		{
			increase();
			return lista[1]-lista[0];
		}
		void restartCount()
		{
			restartCounting();
		}
};

class function6: public AbstractFunction 
{
	public:
		function6():AbstractFunction(){}
		double function(std::vector<double> lista)
		{
			increase();
			return 2-lista[0];
		}
		void restartCount()
		{
			restartCounting();
		}
};

class function7: public AbstractFunction 
{
	public:
		function7():AbstractFunction(){}
		double function(std::vector<double> lista)
		{
			increase();
			return 3-lista[0]-lista[1];
		}
		void restartCount()
		{
			restartCounting();
		}
};

class function8: public AbstractFunction 
{
	public:
		function8():AbstractFunction(){}
		double function(std::vector<double> lista)
		{
			increase();
			return 3+1.5*lista[0]-lista[1];
		}
		void restartCount()
		{
			restartCounting();
		}
};

class function9: public AbstractFunction 
{
	public:
		function9():AbstractFunction(){}
		double function(std::vector<double> lista)
		{
			increase();
			return lista[1]-1;
		}
		void restartCount()
		{
			restartCounting();
		}
};
//transformed function class
class transformed: public AbstractFunction
{
	protected:
		std::vector<AbstractFunction*> container;
		std::vector<AbstractFunction*> implicitFunc;
		std::vector<AbstractFunction*> explicitFunc;
		AbstractFunction* Class;
	public:
		transformed(AbstractFunction* Class, std::vector<AbstractFunction*> implicitFunc)
		{
			this->Class=Class;
			this->implicitFunc=implicitFunc;
		}
		transformed(AbstractFunction* Class, std::vector<AbstractFunction*> implicitFunc, std::vector<AbstractFunction*> explicitFunc)
		{
			this->Class=Class;
			this->implicitFunc=implicitFunc;
			this->explicitFunc=explicitFunc;
		}
		double function(std::vector<double> lista, double r)
		{
			increase();
			double output = this->Class->function(lista);
			for(int i=0;i<this->implicitFunc.size();i++)
			{
				if(this->implicitFunc[i]->function(lista)<=0)
				{
					output += 10e6;				
				}
				else
				{
					output += r*log(this->implicitFunc[i]->function(lista));
				}
			}
			if(this->explicitFunc.size()>0)
			{
				for(int i=0;i<this->explicitFunc.size();i++)
				{
					output += 1/r*pow(this->explicitFunc[i]->function(lista),2);
				}
			}
			//output += 
			return output;
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
void unimodalniMulti(double h, std::vector<double> tocka, double& l, double& r, AbstractFunction& Class, std::vector<double> multiToOne)
{
	double lambda0 = 0;
	//std::cout<<"Ulazni multi-vektor: ";
	//for(int i=0;i<tocka.size();i++) std::cout<<tocka[i]<<" "<<std::endl;
	std::vector<double> vl(tocka.size(),0.0);
	std::vector<double> vr(tocka.size(),0.0);
	double lambdaL = lambda0 - h;
	double lambdaR = lambda0 + h;
	for(int i=0;i<tocka.size();i++) vl[i] = tocka[i]+multiToOne[i]*lambdaL; 
	for(int i=0;i<tocka.size();i++) vr[i] = tocka[i]+multiToOne[i]*lambdaR; 
	std::vector<double> m = tocka;
	//std::cout<<"multi-vektor vl: ";
	//for(int i=0;i<tocka.size();i++) std::cout<<vl[i]<<" "<<std::endl;
	//std::cout<<"multi-vektor vr: ";
	//for(int i=0;i<tocka.size();i++) std::cout<<vr[i]<<" "<<std::endl;
	//std::cout<<"multi-vektor m: ";
	//for(int i=0;i<tocka.size();i++) std::cout<<m[i]<<" "<<std::endl;
	double fl, fm, fr;
	int step = 1;
	double move = 0;

	fm = Class.function(tocka);
	fl = Class.function(vl);
	fr = Class.function(vr);
	
	//std::cout<<"fl"<<" "<<"fm"<<" "<<"fr"<<std::endl;
	//std::cout<<fl<<" "<<fm<<" "<<fr<<std::endl;
	if(fm < fr && fm < fl)
	{
		l = lambdaL;
		r = lambdaR;	
		return;
	}
	else if(fm > fr)
	{
		do
		{	vl = m;
			m = vr;
			fm = fr;
			//vr[dim] = tocka[dim] + h * (step *= 2);
			lambdaR +=  h * (step *= 2);
			for(int i=0;i<tocka.size();i++) vr[i] = tocka[i]+multiToOne[i]*lambdaR;
			fr = Class.function(vr);
			//std::cout<<"multi-vektor vr: ";
			//for(int i=0;i<tocka.size();i++) std::cout<<vr[i]<<" "<<std::endl;
		} while(fm > fr);
		l = lambdaL;
		r = lambdaR;
	}
	else
	{
		do
		{	vr = m;
			m = vl;
			fm = fl;
			//vl[dim] = tocka[dim] - h * (step *= 2);
			lambdaL -=  h * (step *= 2);
			for(int i=0;i<tocka.size();i++) vl[i] = tocka[i]+multiToOne[i]*lambdaL;
			fl = Class.function(vl);
			//std::cout<<"multi-vektor vl: ";
			//for(int i=0;i<tocka.size();i++) std::cout<<vl[i]<<" "<<std::endl;
			//std::cout<<"fl"<<" "<<"fm"<<" "<<"fr"<<std::endl;
			//std::cout<<fl<<" "<<fm<<" "<<fr<<std::endl;
		} while(fm > fl);
		l = lambdaL;
		r = lambdaR;
	}
}

//TREBA RIJESITI UNIMODLANI I ZLATNI REZ!!!!

/*
Algoritam zlatnog reza

ulazne velicine:
- a, b: pocetne granice unimodalnog intervala
- e: preciznost

- t: pocetna tocka
- h: korak

- point: da li je dan unimodalni ili pocetna tocka
*/
void getC(const std::vector<double>& b,const std::vector<double>& a,double k,std::vector<double>& c)
{
	for(int i=0;i<b.size();i++) c[i]=b[i]-k*(b[i]-a[i]);
}

void getD(const std::vector<double>& b,const std::vector<double>& a,double k,std::vector<double>& d)
{
	for(int i=0;i<b.size();i++) d[i]=a[i]+k*(b[i]-a[i]);
}

double dist(const std::vector<double>& b,const std::vector<double>& a)
{
	double sum=0;
	for(int i=0;i<b.size();i++) sum+=(b[i]-a[i])*(b[i]-a[i]);
	return sqrt(sum);
}

//ZLATNI REZ ZA VISE DIMENZIJA
double Zlatni_rezMulti(double h, std::vector<double> tocka,double e, AbstractFunction& Class, std::vector<double> multiToOne)
{	
	double lambda0 = 0;
	//std::vector<double> a,b;
	double lambdaL,lambdaR;
	double lambdaC,lambdaD;
	double k = 0.5*(sqrt(5)-1);
	//izracunaj prvo unimodalni interval
	unimodalniMulti(h,tocka,lambdaL,lambdaR,Class,multiToOne);
	//std::cout<<"Unimodalni interval: ";
	//for(int i=0;i<tocka.size();i++) std::cout<<a<<" ";
	//std::cout<<":";
	//for(int i=0;i<tocka.size();i++) std::cout<<b<<" ";
	//std::cout<<lambdaL<<" : "<<lambdaR<<std::endl;
	
	lambdaC=lambdaR - k* (lambdaR-lambdaL);
	lambdaD=lambdaL + k* (lambdaR-lambdaL);
	
	std::vector<double> c(tocka);
	//getC(lambdaR,lambdaL,k,c);
	for(int i=0;i<tocka.size();i++) c[i]=tocka[i]+multiToOne[i]*(lambdaC);
	//vc[dim]=c;
	std::vector<double> d(tocka);
	//getD(lambdaR,lambdaL,k,c);
	for(int i=0;i<tocka.size();i++) d[i]=tocka[i]+multiToOne[i]*(lambdaD);
	//vd[dim]=d;
	double fc = Class.function(c);
	double fd = Class.function(d);
	while((lambdaR-lambdaL) > e)
	{
		if(fc < fd) {
			//b = d;
			//d = c;			
			lambdaR=lambdaD;
			lambdaD=lambdaC;
			//vc[dim] = b - k * (b - a);
			lambdaC=lambdaR - k* (lambdaR-lambdaL);
			for (int i=0;i<tocka.size();i++) c[i] = tocka[i]+lambdaC*multiToOne[i];
			//std::cout<<"A: "<<lambdaL<<", B: "<<lambdaR<<std::endl;
			fd = fc;
			fc = Class.function(c);
		}
		else
		{
			//a = c;
			//c = d;
			lambdaL=lambdaC;
			lambdaC=lambdaD;
			//vd[dim] = a + k * (b - a);
			lambdaD=lambdaL + k* (lambdaR-lambdaL);
			for (int i=0;i<tocka.size();i++) d[i] = tocka[i]+lambdaD*multiToOne[i];
			//std::cout<<"A: "<<lambdaL<<", B: "<<lambdaR<<std::endl;
			fc = fd;
			fd = Class.function(d);
		}
	}
	//std::vector<double> result(tocka);
	//for(int i=0;i<tocka.size();i++) result[i]=tocka[i]*(lambdaR+lambdaL)/2;
	//std::cout<<"Minimum:";
	//for(int i=0;i<tocka.size();i++) std::cout<<result[i]<<" ";
	//std::cout<<"A: "<<lambdaL<<", B: "<<lambdaR<<std::endl;
	return (lambdaR+lambdaL)/2; // ili nove vrijednosti a i b
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
	return (y2 - y1) / (2*delta);
}

//function to numerically partially derive function
double derivePartially(AbstractFunction& Class, std::vector<double> x0, int coord, double delta = 1.0e-6)
{
	std::vector<double> x1; //= x0 - delta;
	std::vector<double> x2; //= x0 + delta;

	//std::cout<<"x0: ";
	//for(int i=0;i<x0.size();i++) std::cout<<x0[i]<<" ";
	//std::cout<<std::endl;
	//std::cout<<"Coord: "<<coord<<std::endl;
	//std::cout<<"Delta: "<<delta<<std::endl;
	for(int i=0;i<x0.size();i++)
	{
		//std::cout<<i<<" "<<x0[i]<<std::endl;
		if (i==coord)
			{
				x1.push_back(x0[i] - delta);
			}
		else
			{
				x1.push_back(x0[i]);
			}
		//std::cout<<i<<" "<<x1[i]<<std::endl;
	}
	for(int i=0;i<x0.size();i++)
	{
		if (i==coord)
			{
				x2.push_back(x0[i] + delta);
			}
		else
			{
				x2.push_back(x0[i]);
			}
	}
	//for(int i=0;i<x0.size();i++) std::cout<<x1[i]<<" "<<x2[i]<<std::endl;
	double y1 = Class.function(x1);
	double y2 = Class.function(x2);
	//std::cout<<y1<<" "<<y2<<std::endl;
	//std::cout<<(y2 - y1) / (2*delta)<<std::endl;
	return (y2 - y1) / (2*delta);
}

//calculate second order partial derivations using derivePartially for first order
double secondOrderPartials(AbstractFunction& Class, std::vector<double> x0, int firstDer, int secondDer,double delta = 1.0e-6)
{
	std::vector<double> x1(x0); //= x0 - delta;
	std::vector<double> x2(x0); //= x0 + delta;
	//std::cout<<"Delta: "<<delta<<std::endl;
	//std::cout<<"firstDer: "<<firstDer<<std::endl;
	//std::cout<<"secondDer: "<<secondDer<<std::endl;

	for(int i=0;i<x0.size();i++)
	{
		if (i==secondDer)
			{
			x1[i]=x0[i] - delta;
			}
		else
			{
			x1[i]=x0[i];
			}
	}	
	for(int i=0;i<x0.size();i++)
	{
		if (i==secondDer)
		{
			x2[i]=x0[i] + delta;
		}
		else
		{
			x2[i]=x0[i];
		}
	}
	
	double y1 = derivePartially(Class,x1,firstDer,delta);
	double y2 = derivePartially(Class,x2,firstDer,delta);
	//std::cout<<y1<<" "<<y2<<std::endl;
	return (y2 - y1) / (2*delta);	

}


//GRADIENT DESCENT METHOD
void gradientDescent(AbstractFunction& Class, std::vector<double> x0, std::vector<double>& result,int& brSGD, double delta = 1.0e-6, int mode=1)
{
	brSGD=0;
	std::vector<double> partialDerivations(x0.size(),0.0),x_new(x0),x_old(x0);
	double lambda;
	double sumDerivations = 0;
	double sumGrad=0;
	int iter=0;
	double norm=0;
	do {
		norm=0;
		for (int i=0;i<x0.size();i++)
		{
			//std::cout<<i<<std::endl;
			partialDerivations[i]=derivePartially(Class, x_old, i,delta);
			sumDerivations+=partialDerivations[i]*partialDerivations[i];
		}
		brSGD++;
				
		if(mode==1)
		{
			sumDerivations=sqrt(sumDerivations);
			//std::cout<<sumDerivations<<std::endl;
			//normalize vectors
			for (int i=0;i<x0.size();i++)
			{
				partialDerivations[i]=partialDerivations[i]/sumDerivations;
			}
			//std::cout<<"Partial normalized: ";
			//for (int i=0;i<x0.size();i++) std::cout<<partialDerivations[i]<<" ";
			//std::cout<<std::endl;
			lambda = Zlatni_rezMulti(1,x_old,delta,Class,partialDerivations);
			//std::cout<<"Lambda: "<<lambda;
			//for (int i=0;i<x0.size();i++) std::cout<<lambda[i]<<" ";
			for(int i=0;i<x0.size();i++) x_new[i]=x_old[i]+lambda*partialDerivations[i];
		}
		else
		{
			for(int i=0;i<x0.size();i++) x_new[i]=x_old[i]-partialDerivations[i];
		}
		//std::cout<<"NEW X: ";
		//for(int i=0;i<x0.size();i++) std::cout<< x_new[i]<<" ";
		if(Class.function(x_new)>Class.function(x_old))
		{
			iter++;
			if(iter >= 100)
			{
				std::cout<<"It is divergating"<<std::endl;
				return;
			}
		}
		for(int i=0;i<x0.size();i++) norm+=partialDerivations[i]*partialDerivations[i];
		
		if(sqrt(norm)<1.0e-4)
		{
			result=x_old;
			return;
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

//NEWTON-RAPHSON METHOD
void NewtonRaphson(AbstractFunction& Class, std::vector<double> x0, std::vector<double>& result,int& brHessian,int& bgSGD, double delta = 1.0e-6,double eps = 1.0e-6, int mode=1)
{
	brHessian=0;
	bgSGD=0;
	std::vector<double> partialDerivations(x0.size(),0.0),x_old(x0),x(x0),dx(x0),temp_res(x0),dx_n(x0);
	std::vector<std::vector<double>> hessianMatrix,temp;
	for(int i=0;i<x0.size();i++) hessianMatrix.push_back(partialDerivations);
	double sumDerivations=0;
	double sumDx = 0;	
	bool cond=false;
	double norm,lambda;
	for(int i=0;i<x0.size();i++) x[i]=x0[i];
	//std::cout<<"X: "<<std::endl;
	//for(int i=0;i<x0.size();i++) std::cout<<x[i]<<" ";
	int iter=0;
	do
	{
		//Calculate G(x);
		for (int i=0;i<x0.size();i++)
		{
			partialDerivations[i]=derivePartially(Class, x, i, delta);
			sumDerivations+=partialDerivations[i];
		}
		//normalize vectors
		//for (int i=0;i<x0.size();i++)
		//{
		//	partialDerivations[i]/=sumDerivations;
		//}
		//std::cout<<"Normalized vectors: ";
		//for(int i=0;i<x0.size();i++) std::cout<<partialDerivations[i]<<" ";
		//std::cout<<"\n";
		//Calculate J(x); Hessian matrix
		for (int i=0;i<x0.size();i++)
		{
			for(int j=0;j<x0.size();j++)
			{
				hessianMatrix[i][j]=secondOrderPartials(Class,x,i,j,delta);
			}
		}
		brHessian++;
		bgSGD++;
		//resolve J*dx = -G;
		//transform vector of vectors into a matrix
		Matrica J(hessianMatrix);
		//std::cout<<"Hessian matrix:"<<std::endl;
		//J.printMatrix();
		
		for(int i=0;i<partialDerivations.size();i++) partialDerivations[i] = -partialDerivations[i];
		temp.clear();
		temp.push_back(partialDerivations);
		Matrica G(temp);
		Matrica GT=G.transpose();
		//std::cout<<"G matrix:"<<std::endl;
		//GT.printMatrix();
		//Decompose to L and U
		Matrica _4P=J.LUPdekompozicija();
		
		temp_res.clear();
		J.supstitucijaUnaprijed(_4P*GT,temp_res);
		//std::cout<<"temp_res matrix:"<<std::endl;
		//for(int i=0;i<dx.size();i++) std::cout<<temp_res[i]<<" ";
		//std::cout<<std::endl;
		J.supstitucijaUnazad(temp_res,dx);
		//std::cout<<"dx matrix:"<<std::endl;
		//for(int i=0;i<dx.size();i++) std::cout<<dx[i]<<" ";
		//std::cout<<std::endl;
		x_old=x;
		//x -= dx
		if(mode==1)
		{
			for (int i=0;i<x0.size();i++) sumDx+=dx[i]*dx[i];
			for (int i=0;i<x0.size();i++)
			{
				dx_n[i]=dx[i]/sqrt(sumDx);
			}
			//std::cout<<"dx_n matrix:"<<std::endl;
			//for(int i=0;i<dx_n.size();i++) std::cout<<dx_n[i]<<" ";
				
			for(int i=0;i<x.size();i++)
			{
				norm+=dx_n[i]*dx_n[i];	
			}
			norm = sqrt(norm);
			if(norm<eps) break;
		
			lambda = Zlatni_rezMulti(1,x,delta,Class,dx_n);
			//std::cout<<"Lambda: "<<lambda<<std::endl;
			//for (int i=0;i<x0.size();i++) std::cout<<lambda<<" ";
			//std::cout<<"X: "<<std::endl;
			//for(int i=0;i<x.size();i++) std::cout<<x[i]<<" ";
			//std::cout<<std::endl;
			for(int i=0;i<x.size();i++) x[i]=x[i]+lambda*dx_n[i];
			//std::cout<<"NEW X:"<<std::endl;
			//for(int i=0;i<x.size();i++) std::cout<<x[i]<<" ";
		}
		else
		{
			for(int i=0;i<x.size();i++)x[i]+=dx[i];
		}
		norm=0;
		for(int i=0;i<x.size();i++)
		{
			norm+=dx_n[i]*dx_n[i];	
		}
		norm = sqrt(norm);
		if(norm<eps) break;
		if(Class.function(x)>Class.function(x_old))
		{
			iter++;
			if(iter>=100)
			{
				std::cout<<"Divergira!"<<std::endl;
				result=x;
				return;			
			}		
		}
		if(iter>=1) break;
	}
	while(true);
	for(int i=0;i<x.size();i++) result.push_back(x[i]);
	return;
}

//calculate and return centroid of all points
void getCentroid(std::vector<std::vector<double>> array,std::vector<double>& c)
{
	//std::cout<<"Racunanje centroida"<<std::endl;
	std::vector<double> temp;
	//for(int j=0;j<c.size();j++) std::cout<<c[j]<<" ";
	//std::cout<<std::endl;
	for(int j=0;j<c.size();j++) c[j]=0;
	for(int i=0;i<array.size();i++)
	{
			temp=array[i];
			addSame(c,temp);
	}
	//std::cout<<"Array size: "<<array.size()<<std::endl;
	//for(int j=0;j<c.size();j++) std::cout<<c[j]<<" ";
	//std::cout<<std::endl;
	for(int j=0;j<c.size();j++) c[j]/=(array.size());
	
	//for(int j=0;j<c.size();j++) std::cout<<c[j]<<" ";
	//std::cout<<std::endl;
	return;
}

//calculate and return centroid 
void getCentroid(std::vector<std::vector<double>> array,std::vector<double>& c, double h)
{
	//std::cout<<"Racunanje centroida"<<std::endl;
	std::vector<double> temp;
	//for(int j=0;j<c.size();j++) std::cout<<c[j]<<" ";
	//std::cout<<std::endl;
	for(int j=0;j<c.size();j++) c[j]=0;
	for(int i=0;i<array.size();i++)
	{
		if(i!=h)
		{
			temp=array[i];
			addSame(c,temp);
		}
	}
	//std::cout<<"Array size: "<<array.size()<<std::endl;
	//for(int j=0;j<c.size();j++) std::cout<<c[j]<<" ";
	//std::cout<<std::endl;
	for(int j=0;j<c.size();j++) c[j]/=(array.size()-1);
	
	//for(int j=0;j<c.size();j++) std::cout<<c[j]<<" ";
	//std::cout<<std::endl;
	return;
}


//calculate reflex
void refleksija(std::vector<double> xc, std::vector<double> xh,double alfa,std::vector<double>& xr)
{
	//for(int i=0;i<xc.size();i++) std::cout<<xc[i]<<" ";
	//std::cout<<std::endl;
	//for(int i=0;i<xh.size();i++) std::cout<<xh[i]<<" ";
	//std::cout<<std::endl;

	for(int i=0;i<xc.size();i++)
	{
		xr[i]=(1+alfa)*xc[i]-alfa*xh[i];
	}
	return;
}

//return maximum
std::vector<double> getMaximum(std::vector<std::vector<double>> array, AbstractFunction& Class)
{
	std::vector<double> h=array[0];
	for(int i=0;i<array.size();i++)
	{
		if(Class.function(array[i])>Class.function(h)) h = array[i];
	}
	return h;
}

double getValue(std::vector<double> input, AbstractFunction& Class,std::map<std::vector<double>, double>& lookUpTable)
{
	if(lookUpTable.count(input))
	{
		return lookUpTable[input];
	}
	else
	{
		lookUpTable[input]=Class.function(input);
		return lookUpTable[input];
	}
}

//return max index and max2 index
void getMaximumIndex(std::vector<std::vector<double>> array, AbstractFunction& Class,std::map<std::vector<double>, double>& lookUpTable, int& h, int& h2)
{
	//lookUpTable.clear();
	h=0;
	double left,right;
	for(int i=0;i<array.size();i++)
	{
		left = getValue(array[i],Class,lookUpTable);
		right = getValue(array[h],Class,lookUpTable);
		if(left>right) h = i;
	}
	std::vector<double> values;
	for(std::map<std::vector<double>,double>::iterator it = lookUpTable.begin(); it != lookUpTable.end(); ++it) {
		values.push_back(it->second);
	}
	std::sort(values.begin(),values.end());
	double max = values[values.size()-1];
	double max2 = values[values.size()-2];
	for(int i=0;i<array.size();i++)
	{
		left = getValue(array[i],Class,lookUpTable);
		//if(left==max)  h = i;
		if(left==max2) h2 = i;
	}
	return;
}

//function to check the implicit values of a point
bool checkImplicit(std::vector<AbstractFunction*> implicitList,std::vector<double> x0)
{
	for(int i=0;i<implicitList.size();i++)
	{
		if(implicitList[i]->function(x0)<0)
		{
			return false;
		}
	}
	return true;
}

void printVector(std::vector<double> x)
{
	for(int i=0;i<x.size();i++)
	{
		std::cout<<x[i]<<" ";
	}
	std::cout<<std::endl;
}

//izracunaj ulazni skup tocki simpleksa
void startSimpleks(std::vector<double> x0,std::vector<double> explicitConditions,std::vector<AbstractFunction*> implicitList,std::vector<std::vector<double>>& array)
{
	//AbstractFunction* temp;
	double r;
	double temp;
	std::vector<double> tempX(x0.size(),explicitConditions[0]);
	std::vector<double> xc(x0);
	//std::cout<<"Function has been entered!"<<std::endl;
	//check explicit conditions
	//std::cout<<"Start vector size: "<<x0.size()<<std::endl;
	for(int i=0;i<x0.size();i++)
	{
		//std::cout<<i<<std::endl;
		if(x0[i]<explicitConditions[0])
		{
			x0[i]=explicitConditions[0];
		}
		if(x0[i]>explicitConditions[1]) 
		{
			x0[i]=explicitConditions[1];
		}
	}
	//std::cout<<"why does function break here?"<<std::endl;
	//for(int i=0;i<x0.size();i++) std::cout<<x0[i]<<" ";
	//check implicit conditions
	//std::cout<<"Number of implicit conditions: "<<implicitList.size()<<std::endl;
	for(int i=0;i<implicitList.size();i++)
	{
		temp=implicitList[i]->function(x0);
		//std::cout<<"Test: "<<temp<<std::endl;
		if(temp<0)
		{
			//NOT GOOD
			std::cout<<"Starting point does not satisfy implicit conditions."<<std::endl;
			return;
		}
		else
		{
			continue;
		}
	}
	for(int i=0;i<2*x0.size();i++)
	{
		tempX.clear();
		tempX.assign(x0.size(),explicitConditions[0]);
		/*
		for(int j=0;j<x0.size();j++)
		{
			do
			{
				r = ((double) rand() / (RAND_MAX));
				std::cout<<r<<std::endl;
				tempX[j] = explicitConditions[0] + r*(explicitConditions[1]-explicitConditions[0]);
				printVector(tempX);
				
			}
			while(!checkImplicit(implicitList,tempX));
		}
		*/
		do
		{
			for(int j=0;j<tempX.size();j++)
			{
				r = ((double) rand() / (RAND_MAX));
				//std::cout<<r<<std::endl;
				tempX[j] = explicitConditions[0] + r*(explicitConditions[1]-explicitConditions[0]);
				//printVector(tempX);
			}

		}
		while(!checkImplicit(implicitList,tempX));
		//printVector(tempX);
		for(int j=0;j<x0.size();j++)
		{
			tempX[j] = 0.5 * (tempX[j] + xc[j]);
		}
		//printVector(tempX);
		array.push_back(tempX);
		getCentroid(array,xc);
		//printVector(xc);
	}
	std::cout<<"Simplex size is: "<<array.size()<<std::endl;
	return;
}
bool checkStopingCriteria(std::vector<double> xc, std::vector<double> xh, double eps)
{
	double sum=0;
	for(int i=0;i<xc.size();i++)
	{
		sum+=(abs(xh[i]-xc[i]))*(abs(xh[i]-xc[i]));
	}
	if(sqrt(sum)<eps) return true;
	return false;
}

void printArray(std::vector<std::vector<double>> array)
{
	for(int i=0;i<array.size();i++)
	{
		printVector(array[i]);
	}
}

void Box(AbstractFunction& Class, std::vector<double> x0, std::vector<double>& result,std::vector<double> explicitConditions,std::vector<AbstractFunction*> implicitList,double alfa, double delta = 1.0e-6,double eps = 1.0e-6, int mode=1)
{
	std::vector<std::vector<double>> array;
	std::map<std::vector<double>, double> lookUpTable;
	std::vector<double> xc(x0.size(),0.0),xr(x0.size(),0.0),xc_old(x0.size(),0.0),xh(x0.size(),0.0);
	int iter=0;
	int h,h2;
	int z=0;
	//std::cout<<"Populating start array\n";
	std::cout<<"Alfa is: "<<alfa<<std::endl;
	//Calculate starting simplex
	startSimpleks(x0, explicitConditions, implicitList,array);
	for(int i=0;i<array.size();i++)
	{
		printVector(array[i]);
	}
	xc_old = x0;
	do
	{
		z++;
		std::cout<<"Iteracija: "<<z<<std::endl;
		getMaximumIndex(array,Class,lookUpTable,h,h2);	
		getCentroid(array,xc,h);
		//std::cout<<"Checkpoint 1\n";
		xh = array[h];
		refleksija(xc,xh,alfa,xr);
		//std::cout<<"Center: ";
		//printVector(xc);
		//std::cout<<"Reflex: ";
		//printVector(xr);
		//explicit conditions
		for(int i=0;i<xc.size();i++)
		{
			if(xr[i]<explicitConditions[0]) xr[i]=explicitConditions[0];
			if(xr[i]>explicitConditions[1]) xr[i]=explicitConditions[1];
		}
		//std::cout<<"Checkpoint 2\n";
		while(!checkImplicit(implicitList,xr))
		{	
			std::cout<<"Doing implicit"<<std::endl;
			for(int j=0;j<xr.size();j++)
			{
				xr[j] = 0.5 * (xr[j] + xc[j]);
			}
		}
		//std::cout<<"Checkpoint 3\n";
		if(Class.function(xr)>Class.function(array[h2]))
		{
			for(int j=0;j<xr.size();j++)
			{
				xr[j] = 0.5 * (xr[j] + xc[j]);
			}
		}
		array[h] = xr;
		//std::cout<<"Checkpoint 4\n";
		if(getValue(xc,Class,lookUpTable)>getValue(xc_old,Class,lookUpTable))
		{
			iter++;
			if(iter>=100)
			{
				std::cout<<"The solution is divering!"<<std::endl;
				return;
			}
		}
		//std::cout<<"Checkpoint 5\n";
		xc_old=xc;
		if(checkStopingCriteria(xc,xr,eps)) 
		{
			//std::cout<<"Why me?\n";
			result=xc;
			break;
		}
		printArray(array);
		//if(z>=4) break;
	}
	while(true);
	for(int i=0;i<xc.size();i++) result[i]=xc[i];
}

// ---------------------------------------------------------------------- LABOS 2 NELDER MEAD

//izracunaj ulazni skup tocki simpleksa
std::vector<std::vector<double>> tockeSimpleksa(std::vector<double> x0,double t)
{
	std::vector<double> ref(x0.size(),0.0);
	std::vector<std::vector<double>> array(x0.size()+1,ref);
	double a1 = t*(sqrt(x0.size()+1)+x0.size()-1)/(x0.size()*sqrt(2));
	double a2 = t*(sqrt(x0.size()+1)-1)/(x0.size()*sqrt(2));
	//std::cout<<"a1: "<<a1<<"\na2: "<<a2<<std::endl;
	std::vector<double> temp(x0.size(),a2);
	std::vector<double> temp2 = temp;
	//for(int j=0;j<temp2.size();j++) std::cout<<temp2[j]<<" ";
	//std::cout<<"<-That was prototype"<<std::endl;
	for(int i=0;i<x0.size();i++)
	{
		//std::cout<<i<<" ";
		temp2 = temp;
		temp2[i]=a1;
		//for(int j=0;j<temp2.size();j++) std::cout<<temp2[j]<<" ";
		//std::cout<<std::endl;
		addSame(array[i],temp2);
		addSame(array[i],x0);
	}
	array[array.size()-1]=x0;

	return array;
}

//vrati index najvece ulazne vrijednosti
int getMaximumIndex(std::vector<std::vector<double>> array, AbstractFunction& Class,std::map<std::vector<double>, double>& lookUpTable)
{
	int h=0;
	double left,right;
	for(int i=0;i<array.size();i++)
	{
		left = getValue(array[i],Class,lookUpTable);
		right = getValue(array[h],Class,lookUpTable);
		if(left>right) h = i;
	}
	return h;
}

//vrati ulaz koji daje najmanju vrijednost
std::vector<double> getMinimum(std::vector<std::vector<double>> array, AbstractFunction& Class)
{
	std::vector<double> l=array[0];
	for(int i=0;i<array.size();i++)
	{
		if(Class.function(array[i])<Class.function(l)) l = array[i];
	}
	return l;
}

//vrati index najmanje ulazne vrijednosti
int getMinimumIndex(std::vector<std::vector<double>> array, AbstractFunction& Class,std::map<std::vector<double>, double>& lookUpTable)
{
	int l=0;
	double left,right;
	for(int i=0;i<array.size();i++)
	{
		left = getValue(array[i],Class,lookUpTable);
		right = getValue(array[l],Class,lookUpTable);
		if(left<right) l = i;
	}
	return l;
}

//izracunaj ekspazivnu tocku
void ekspanzija(std::vector<double> xc, std::vector<double> xh,double gamma,double alfa,std::vector<double>& xe)
{
	std::vector<double> xr(xc.size(),0.0);
	refleksija(xc,xh,alfa,xr);
	for(int i=0;i<xc.size();i++)
	{
		xe[i]=(1-gamma)*xc[i]+gamma*xr[i];
	}
	return;
}

//izracunaj kontrakcijsku tocku
void kontrakcija(std::vector<double> xc, std::vector<double> xh,double beta,std::vector<double>& xk)
{
	for(int i=0;i<xc.size();i++)
	{
		xk[i]=(1-beta)*xc[i]+beta*xh[i];
	}
	return;
}

void pomakiPremaL(std::vector<std::vector<double>>& array,int l)
{
	for(int i=0;i<array.size();i++)
	{
		for(int j=0;j<array[i].size();j++)
		{
			array[i][j] = 0.5 * (array[i][j]+array[l][j]);
		}
	}
}


bool kriterijZaustavljanja(std::vector<std::vector<double>> array, std::vector<double> xc,double eps,AbstractFunction& Class,std::map<std::vector<double>, double>& lookUpTable)
{
	double fxc = getValue(xc,Class,lookUpTable);
	double sum=0;
	double temp;
	for(int i=0;i<array.size();i++)
	{
		temp = getValue(array[i],Class,lookUpTable);
		sum += pow((temp-fxc),2);
	}
	sum/=array.size();
	sum=sqrt(sum);
	if(sum<eps) return true;
	else return false;
}

/*
Nelder-Mead simpleks algoritam
Ulazne velicine: X0,razmak, alfa, beta, gama, epsilon,funkcija
*/

std::vector<double> NelderMead(std::vector<double> x0,double razmak, double alfa, double beta, double gamma, double epsilon, AbstractFunction& Class)
{
	std::vector<std::vector<double>> array;
	std::vector<double> xc(x0.size(),0.0),xr(x0.size(),0.0),xk(x0.size(),0.0),xe(x0.size(),0.0);
	int h,l;
	std::map<std::vector<double>, double> lookUpTable;
	bool checkVar = true;
	bool checkStop = false;
	
	//Izracunaj pocetne tocke simpleksa
	array=tockeSimpleksa(x0,razmak);
		
	do
	{
		h = getMaximumIndex(array,Class,lookUpTable);
		l = getMinimumIndex(array,Class,lookUpTable);
		//std::cout<<"Maks index: "<<h<<" Min index "<<l<<std::endl;	
		getCentroid(array,xc,h);
		//std::cout<<"Velicina centroida: "<<xc.size()<<std::endl;
		
		refleksija(xc,array[h],alfa,xr);
		//std::cout<<"Centroid: "<<std::endl;
		//for(int i=0;i<xc.size();i++) std::cout<<xc[i]<<" ";
		//std::cout<<std::endl;
		//std::cout<<"Refleksija: "<<std::endl;
		//for(int i=0;i<xc.size();i++) std::cout<<xr[i]<<" ";
		//std::cout<<std::endl;
		
		//std::cout<<"Vrijednost refleksije: "<<Class.function(xr)<<std::endl;
		//std::cout<<"Vrijednost najmanje: "<<Class.function(array[l])<<std::endl;
		
		//if(Class.function(xr) < Class.function(array[l]))
		if(getValue(xr,Class,lookUpTable) < getValue(array[l],Class,lookUpTable))
		{
			//std::cout<<"I am in if!"<<std::endl;
			ekspanzija(xc,array[h],gamma,alfa,xe);
			//if(Class.function(xe)<Class.function(array[l])) array[h] = xe;
			if(getValue(xe,Class,lookUpTable)<getValue(array[l],Class,lookUpTable)) array[h] = xe;
			else array[h] = xr;
		}
		else
		{
			checkVar=true;
			//std::cout<<"I am in else!"<<std::endl;
			for(int j=0;j<array.size();j++)
			{
				if(j!=h)
				{
					//if(Class.function(xr) <= Class.function(array[j])) checkVar=false;
					if(getValue(xr,Class,lookUpTable) <= getValue(array[j],Class,lookUpTable)) checkVar=false;
				}
			}
			if(checkVar)
			{
				//if(Class.function(xr) < Class.function(array[h])) array[h]=xr;
				if(getValue(xr,Class,lookUpTable) < getValue(array[h],Class,lookUpTable)) array[h]=xr;
				kontrakcija(xc,array[h],beta,xk);
				//if(Class.function(xk) < Class.function(array[h])) array[h]=xk;
				if(getValue(xk,Class,lookUpTable) < getValue(array[h],Class,lookUpTable)) array[h]=xk;
				else pomakiPremaL(array,l);
				
			}
			else array[h] = xr;
		}
		
		//kriterijZaustavljanja(array,xc,krit);
		checkStop = kriterijZaustavljanja(array,xc,epsilon,Class,lookUpTable);;//compareVectors(krit,epsilon);
	}
	while(!checkStop);
	return array[l];
	
	//return x0;
}

// ---------------------------------------------------------------------------------

//transformed without explicit
void transformationNM(AbstractFunction& Class, std::vector<double> x0, std::vector<double>& result,std::vector<AbstractFunction*> implicitList,double pomak,double alfa,double beta, double gamma,double eps = 1.0e-6)
{

	transformed theFunction(&Class,implicitList);
	std::vector<double> temp = NelderMead(x0,pomak,alfa,beta,gamma,eps,theFunction);
	result=temp;
	return;
}

//transformed with explicit
void transformationNM(AbstractFunction& Class, std::vector<double> x0, std::vector<double>& result,std::vector<AbstractFunction*> implicitList,std::vector<AbstractFunction*> explicitList,double pomak,double alfa,double beta, double gamma,double eps = 1.0e-6)
{

	transformed theFunction(&Class,implicitList,explicitList);
	std::vector<double> temp = NelderMead(x0,pomak,alfa,beta,gamma,eps,theFunction);
	for(int i=0;i<x0.size();i++) result[i]=temp[i];
	return;
}

void openFile(std::string name,std::vector<double>& tocka,std::vector<double>& minimumFunkcije,std::vector<double>& preciznost, std::vector<double>& pomaciFunkcije,double& leftPoint,double& rightPoint,double& distance, int& mode, double& alfa)
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
			if (container[0] == "Mode:")
			{
				mode=atof(container[1].c_str());		
			}
			if (container[0] == "Alfa:")
			{
				alfa=atof(container[1].c_str());		
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
	double leftPoint,rightPoint,distance,alfa;
	int zadatak;
	int mode;
	//std::cout<<argc<<std::endl;
	if(argc<3)
	{
		std::cout<<"Nisi zadao ulazni txt file ili zadatak koji se treba vrtiti\nPrekidam program."<<std::endl;
		std::exit(0);
	}
	zadatak=atof(argv[2]);
	openFile(argv[1],tocka,minimumFunkcije,preciznost,pomaciFunkcije,leftPoint,rightPoint,distance,mode,alfa);
	
	//std::cout<<"Alfa is: "<<alfa;
	//double alfa=1;
	double beta=0.5;
	double gamma=2;	
	
	if(zadatak==1)
	{
		function3 func3;
		std::vector<double> rezultat;
		int brCall;
		//std::cout<<func3.function(tocka)<<std::endl;
		std::cout<<"Gradient descent function 3:"<<std::endl;	
		gradientDescent(func3, tocka, rezultat,brCall,preciznost[0],mode);
		std::cout<<"Minimum: ";
		for(int k=0;k<rezultat.size();k++) std::cout<<std::setw(5)<<rezultat[k]<<" ";
		std::cout<<std::endl;
		std::cout<<"Broj poziva funkcije: "<<func3.getNumbers()<<std::endl;
		std::cout<<"Broj poziva SGD: "<<brCall<<std::endl;
		
		func3.restartCount();
		
	}
	if(zadatak==2)
	{
		function1 func1;
		function2 func2;
		int brCall1SGD,brCall2SGD;
		int brCall1H,brCall2H;
	
		std::vector<double> rezultat1;
		std::vector<double> rezultat2;

		
		//Gradient descent
		std::cout<<"Gradient descent function 1:"<<std::endl;
		gradientDescent(func1, tocka, rezultat1,brCall1SGD,preciznost[0],mode);
		std::cout<<"Minimum: ";
		for(int k=0;k<rezultat1.size();k++) std::cout<<std::setw(5)<<rezultat1[k]<<" ";
		std::cout<<std::endl;
		std::cout<<"Broj poziva funkcije: "<<func1.getNumbers()<<std::endl;
		std::cout<<"Broj poziva SGD: "<<brCall1SGD<<std::endl;
		func1.restartCount();
		
		std::cout<<"Gradient descent function 2:"<<std::endl;
		gradientDescent(func2, tocka, rezultat2,brCall2SGD,preciznost[0],mode);
		std::cout<<"Minimum: ";
		for(int k=0;k<rezultat2.size();k++) std::cout<<std::setw(5)<<rezultat2[k]<<" ";
		std::cout<<std::endl;
		std::cout<<"Broj poziva funkcije: "<<func2.getNumbers()<<std::endl;
		std::cout<<"Broj poziva SGD: "<<brCall2SGD<<std::endl;
		func2.restartCount();
		
		rezultat1.clear();
		rezultat2.clear();
		
		//NewtonRaphson
		std::cout<<"Newton-Raphson function 1:"<<std::endl;
		NewtonRaphson(func1, tocka, rezultat1,brCall1H,brCall1SGD,preciznost[0],preciznost[0],mode);
		std::cout<<"Minimum: ";
		for(int k=0;k<rezultat1.size();k++) std::cout<<std::setw(5)<<rezultat1[k]<<" ";
		std::cout<<std::endl;
		std::cout<<"Broj poziva funkcije: "<<func1.getNumbers()<<std::endl;
		std::cout<<"Broj poziva SGD: "<<brCall1SGD<<std::endl;
		std::cout<<"Broj poziva Hessian: "<<brCall1H<<std::endl;
		func1.restartCount();
		
		std::cout<<"Newton-Raphson function 2:"<<std::endl;
		NewtonRaphson(func2, tocka, rezultat2,brCall2H,brCall2SGD,preciznost[0],preciznost[0],mode);
		std::cout<<"Minimum: ";
		for(int k=0;k<rezultat2.size();k++) std::cout<<std::setw(5)<<rezultat2[k]<<" ";
		std::cout<<std::endl;
		std::cout<<"Broj poziva funkcije: "<<func2.getNumbers()<<std::endl;
		std::cout<<"Broj poziva SGD: "<<brCall2SGD<<std::endl;
		std::cout<<"Broj poziva Hessian: "<<brCall2H<<std::endl;
		func2.restartCount();
		
	}
	if(zadatak==3)
	{
		function1 func1;
		function2 func2;

		std::vector<double> rezultat1;
		std::vector<double> rezultat2;

		function5 impl1;
		function6 impl2;
		std::vector<AbstractFunction*> implicitCondition;
		implicitCondition.push_back(&impl1);	
		implicitCondition.push_back(&impl2);
		
		std::vector<double> explicitCondition;
		explicitCondition.push_back(-100);
		explicitCondition.push_back(100);
		std::cout<<"Doing box algorithm\n";
		
		std::cout<<"Box function 1:"<<std::endl;
		Box(func1,tocka,rezultat1,explicitCondition,implicitCondition,alfa,preciznost[0],preciznost[0],mode);
		std::cout<<"Minimum: ";
		for(int k=0;k<rezultat1.size();k++) std::cout<<std::setw(5)<<rezultat1[k]<<" ";
		std::cout<<"Broj poziva: "<<func1.getNumbers()<<std::endl;
		func1.restartCount();
		
		std::cout<<"Box function 2:"<<std::endl;
		Box(func2,tocka,rezultat2,explicitCondition,implicitCondition,alfa,preciznost[0],preciznost[0],mode);
		std::cout<<"Minimum: ";
		for(int k=0;k<rezultat2.size();k++) std::cout<<std::setw(5)<<rezultat2[k]<<" ";
		std::cout<<"Broj poziva: "<<func2.getNumbers()<<std::endl;
		func2.restartCount();		
	}
	if(zadatak==4)
	{
		function1 func1;
		function2 func2;

		std::vector<double> rezultat1;
		std::vector<double> rezultat2;

		function5 impl1;
		function6 impl2;
		std::vector<AbstractFunction*> implicitCondition;
		implicitCondition.push_back(&impl1);	
		implicitCondition.push_back(&impl2);
		
		//transformed Tfunc(&func2,implicitCondition);
		//std::cout<<Tfunc.function(tocka,1);
		std::cout<<"transformation function 1:"<<std::endl;
		transformationNM(func1,tocka,rezultat1,implicitCondition,1,alfa,beta,gamma,preciznost[0]);
		std::cout<<"Minimum: ";
		for(int k=0;k<rezultat1.size();k++) std::cout<<std::setw(5)<<rezultat1[k]<<" ";
		std::cout<<"Broj poziva: "<<func1.getNumbers()<<std::endl;
		func1.restartCount();


		std::cout<<"transformation function 2:"<<std::endl;
		transformationNM(func2,tocka,rezultat2,implicitCondition,1,alfa,beta,gamma,preciznost[0]);
		std::cout<<"Minimum: ";
		for(int k=0;k<rezultat2.size();k++) std::cout<<std::setw(5)<<rezultat2[k]<<" ";
		std::cout<<"Broj poziva: "<<func2.getNumbers()<<std::endl;
		func2.restartCount();	

	}
	if(zadatak==5)
	{
		function4 func4;

		function7 func7;
		function8 func8;
		function9 func9;

		std::vector<double> rezultat4;
		
		std::vector<AbstractFunction*> implicitCondition;
		implicitCondition.push_back(&func7);
		implicitCondition.push_back(&func8);

		std::vector<AbstractFunction*> explicitCondition;
		explicitCondition.push_back(&func9);

		std::cout<<"transformation function 4:"<<std::endl;
		transformationNM(func4,tocka,rezultat4,implicitCondition,explicitCondition,1,alfa,beta,gamma,preciznost[0]);
		std::cout<<"Minimum: ";
		for(int k=0;k<rezultat4.size();k++) std::cout<<std::setw(5)<<rezultat4[k]<<" ";
		std::cout<<"Broj poziva: "<<func4.getNumbers()<<std::endl;
		func4.restartCount();	
		
	}
	return 0;

	
}
