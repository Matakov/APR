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

//function for subtracting vectors
void subtractSame(std::vector<double>& a,std::vector<double> b)
{
	for(int i=0;i<a.size();i++)
	{
		a[i]-=b[i];	
	}
}
std::vector<double> absoluteValue(std::vector<double> a)
{
	std::vector<double> c;
	for(int i=0;i<a.size();i++)
	{
		c[i]=abs(a[i]);	
	}
	return c;
}

std::vector<double> subtract(std::vector<double> a,std::vector<double> b)
{
	std::vector<double> c;
	for(int i=0;i<a.size();i++)
	{
		c[i]=a[i]-b[i];	
	}
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

double sumSquareElements(std::vector<double> a)
{
	double sum=0;
	for(int i=0;i<a.size();i++) sum+=pow(a[i],2);
	return sum;
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
void unimodalniMulti(double h, std::vector<double> tocka, double& l, double &r, AbstractFunction& Class, int dim)
{
	for(int i=0;i<tocka.size();i++)
	{
		if(i!=dim) tocka[i]=0;
	}
	std::vector<double> vl(tocka.size(),0.0);
	std::vector<double> vr(tocka.size(),0.0);
	vl[dim] = tocka[dim] - h, vr[dim] = tocka[dim] + h; 
	std::vector<double> m = tocka;
	double fl, fm, fr;
	int step = 1;

	fm = Class.function(tocka);
	fl = Class.function(vl);
	fr = Class.function(vr);

	if(fm < fr && fm < fl)
	{
	l = vl[dim];
	r = vr[dim];	
	return;
	}
	else if(fm > fr)
		do
		{	vl = m;
			m = vr;
			fm = fr;
			vr[dim] = tocka[dim] + h * (step *= 2);
			fr = Class.function(vr);
		} while(fm > fr);
	else 
		do
		{	vr = m;
			m = vl;
			fm = fl;
			vl[dim] = tocka[dim] - h * (step *= 2);
			fl = Class.function(vl);
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

double Zlatni_rez(bool point,double h, double t,double e, AbstractFunction& Class)
{	
	double a,b;
	if(point)
	{
		//izracunaj prvo unimodalni interval
		unimodalni(h,t,a,b,Class);
	}
	else
	{
		a = h;
		b = t;
	}
	double k = 0.5*(sqrt(5)-1);
	double c = b - k * (b - a);
	double d = a + k * (b - a);
	double fc = Class.function(c);
	double fd = Class.function(d);
	while((b - a) > e)
	{
		if(fc < fd) {
			b = d;
			d = c;
			c = b - k * (b - a);
			fd = fc;
			fc = Class.function(c);
		}
		else
		{
			a = c;
			c = d;
			d = a + k * (b - a);
			fc = fd;
			fd = Class.function(d);
		}
	}
	return (a + b)/2; // ili nove vrijednosti a i b
}

//ZLATNI REZ ZA VISE DIMENZIJA
double Zlatni_rezMulti(bool point,double h, std::vector<double> t,double e, AbstractFunction& Class, int dim)
{	
	double a,b;
	if(point)
	{
		//izracunaj prvo unimodalni interval
		unimodalniMulti(h,t,a,b,Class,dim);
	}
	else
	{
		a = h;
		b = t[dim];
	}
	double k = 0.5*(sqrt(5)-1);
	std::vector<double> vc(t.size(),0.0);
	double c = b - k * (b - a);
	vc[dim]=c;
	std::vector<double> vd(t.size(),0.0);
	double d = a + k * (b - a);
	vd[dim]=d;
	double fc = Class.function(vc);
	double fd = Class.function(vd);
	while((b - a) > e)
	{
		if(fc < fd) {
			b = vd[dim];
			vd[dim] = vc[dim];
			vc[dim] = b - k * (b - a);
			fd = fc;
			fc = Class.function(vc);
		}
		else
		{
			a = vc[dim];
			vc[dim] = vd[dim];
			vd[dim] = a + k * (b - a);
			fc = fd;
			fd = Class.function(vd);
		}
	}
	return (a + b)/2; // ili nove vrijednosti a i b
}

//see if the left vector is smaller then the right
bool compareVectors(std::vector<double>a, std::vector<double> b)
{
	for(int i=0;i<a.size();i++)
	{
		if(a[i]>b[i]) return false;
	}
	return true;
}

/*
Pretraživanje po koordinatnim osima
*/
std::vector<double> KoordiantneOsi(std::vector<double> x0, std::vector<double> eps,AbstractFunction& Class)
{
	std::vector<double> x = x0;
	std::vector<double> ref,xS;
	double temp;
	double e;
	bool check;
	do{
		xS = x;
		for(int i=0;i<x.size();i++)
		{
			e=eps[i];
			temp=Zlatni_rezMulti(true,1,x,e,Class,i);
			x[i]=temp;
		}
		ref = subtract(x,xS);
		ref = absoluteValue(ref);
		check = compareVectors(ref,eps);
	}while(check);
	return x;
}


//izracunaj ulazni skup tocki simpleksa
std::vector<std::vector<double>> tockeSimpleksa(std::vector<double> x0,double t)
{
	std::vector<std::vector<double>> array;
	double a1 = t*(sqrt(x0.size()+1)+x0.size()-1)/(x0.size()*sqrt(2));
	double a2 = t*(sqrt(x0.size()+1)-1)/(x0.size()*sqrt(2));
	for(int i=0;i<=x0.size();i++)
	{
		std::vector<double> temp(x0.size(),a2);
		if(i<x0.size())	temp[i]=a1;
		array.push_back(temp);
	}
	return array;
}

//izracunaj i vrati centroid skupa vrijednosti
std::vector<double> getCentroid(std::vector<std::vector<double>> array)
{
	std::vector<double> c(array.size()-1,0.0);
	for(int i=0;i<(array.size()-1);i++)
	{
		for(int j=0;j<array.size();j++)
		{
			c[i]+=array[j][i];
		}
		c[i]/=(array.size()-1);
	}
}

//vrati ulaz koji daje najvecu vrijednost
std::vector<double> getMaximum(std::vector<std::vector<double>> array, AbstractFunction& Class)
{
	std::vector<double> h=array[0];
	for(int i=0;i<array.size();i++)
	{
		if(Class.function(array[i])>Class.function(h)) h = array[i];
	}
	return h;
}

//vrati index najvece ulazne vrijednosti
int getMaximumIndex(std::vector<std::vector<double>> array, AbstractFunction& Class)
{
	int h=0;
	for(int i=0;i<array.size();i++)
	{
		if(Class.function(array[i])>Class.function(array[h])) h = i;
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
int getMinimumIndex(std::vector<std::vector<double>> array, AbstractFunction& Class)
{
	int l=0;
	for(int i=0;i<array.size();i++)
	{
		if(Class.function(array[i])<Class.function(array[l])) l = i;
	}
	return l;
}

//izracunaj refleksivnu tocku
std::vector<double> refleksija(std::vector<double> xc, std::vector<double> xh,double alfa)
{
	std::vector<double> xr;
	for(int i=0;i<xc.size();i++)
	{
		xr[i]=(1+alfa)*xc[i]-alfa*xh[i];
	}
	return xr;
}

//izracunaj ekspazivnu tocku
std::vector<double> ekspanzija(std::vector<double> xc, std::vector<double> xh,double gamma,double alfa)
{
	std::vector<double> xr=refleksija(xc,xh,alfa);
	std::vector<double> xe;
	for(int i=0;i<xc.size();i++)
	{
		xe[i]=(1-gamma)*xc[i]+gamma*xr[i];
	}
	return xe;
}

//izracunaj kontrakcijsku tocku
std::vector<double> kontrakcija(std::vector<double> xc, std::vector<double> xh,double beta)
{
	std::vector<double> xk;
	for(int i=0;i<xc.size();i++)
	{
		xk[i]=(1-beta)*xc[i]+beta*xh[i];
	}
	return xk;
}

void pomakiPremaL(std::vector<std::vector<double>>& array,int l)
{
	for(int i=0;i<array.size();i++)
	{
		if(i==l) continue;
		for(int j=0;j<array[i].size();j++)
		{
			array[i][j] = 0.5 * (array[i][j]+array[l][j]);
		}
	}
}

std::vector<double> kriterijZaustavljanja(std::vector<std::vector<double>> array, std::vector<double> xc)
{
	std::vector<double> exitVector;
	for(int i=0;i<array.size()-1;i++)
	{
		for(int j=0;j<array.size();j++)
		{
			exitVector[i]+=pow((array[j][i]-xc[i]),2);
		}
		exitVector[i]=sqrt(exitVector[i]/array.size());
	}
	return exitVector;	
}

/*
Nelder-Mead simpleks algoritam
Ulazne velicine: X0, alfa, beta, gama, epsilon
*/
std::vector<double> NelderMead(std::vector<double> x0,double razmak, double alfa, double beta, double gamma, std::vector<double> epsilon, AbstractFunction& Class)
{
	std::vector<std::vector<double>> array = tockeSimpleksa(x0,razmak);
	std::vector<double> xc,xr,xk,xe,krit;
	int h,l;
	bool checkVar = true;
	bool checkStop = false;
	do
	{
		h = getMaximumIndex(array,Class);
		l = getMinimumIndex(array,Class);	
		xc = getCentroid(array);
		xr = refleksija(xc,array[h],alfa);
		if(Class.function(xr) < Class.function(array[l]))
		{
			xe = ekspanzija(xc,array[h],gamma,alfa);
			if(Class.function(xe)<Class.function(array[l])) array[h] = xe;
			else array[h] = xr;
		}
		else
		{
			for(int j=0;j<array.size();j++)
			{
				if(j==h) continue;
				if(Class.function(xr) <= Class.function(array[j])) checkVar=false;
			}
			if(checkVar)
			{
				if(Class.function(xr) < Class.function(array[h])) array[h]=xr;
				xk = kontrakcija(xc,array[h],beta);
				if(Class.function(xk) < Class.function(array[h])) array[h]=xk;
				else pomakiPremaL(array,l);
				
			}
			else array[h] = xr;
		}
		
		krit = kriterijZaustavljanja(array,xc);
		checkStop = compareVectors(krit,epsilon);
	}
	while(checkStop);
}


//Explore function for Hooks-Jeeves algorithm
std::vector<double> explore(std::vector<double> xP, std::vector<double> Dx,AbstractFunction& Class)
{
	std::vector<double> x=xP;
	double P,N;
	for(int i=0;i<x.size();i++)
	{
		P = Class.function(x);
		x[i] = x[i]+Dx[i];
		N = Class.function(x);
		if(N>P)
		{
			x[i]=x[i]-2*Dx[i];
			N = Class.function(x);
			if(N>P)
			{
				x[i]=x[i] + Dx[i];
			}
		}	
	}
	return x;
}

/*
Algoritam Hooke-Jeeves postupka
x0 - pocetna tocka
xB - bazna tocka 
xP - pocetna tocka pretrazivanja
xN - tocka dobivena pretrazivanjem
*/
/*
std::vector<double> HookeJeeves(std::vector<double> x0, std::vector<double> precision,std::vector<double> Dx,AbstractFunction& Class)
{
	std::vector<double> xB=x0;
	std::vector<double> xP=x0;
	std::vector<double> xN;
	do
	{
		xN = explore(xP,Dx,Class);
		if(Class.function(xN)<Class.function(xP))
		{
			xN = multiply(xN,2);
			xP = subtract(xN,xB);
			xB=xN;
		}
		else
		{
			Dx = divide(Dx,2);
			xP = xB;		
		}
	}while();
	return xB;

}
*/

void myfunction (std::string i) {  // function:
  std::cout << ' ' << i <<std::endl;
}





int main(int argc, char* argv[]){
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
	std::cout<<Zlatni_rez(true,h,tocka,preciznost,func3)<<std::endl;
	return 0;
}
