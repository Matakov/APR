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
	private:
		double b;
		std::vector<double> bList;
	public:
		function3():AbstractFunction(){}
		function3(double input):b(input),AbstractFunction(){}
		function3(std::vector<double> input):bList(input),AbstractFunction(){}
		/*
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
		*/
		double function(std::vector<double> lista)
		{
			increase();
			int numDim = lista.size();
			double sum = 0;
			for(int i=0;i<numDim;i++)
			{
				sum +=pow((lista[i]-bList[i]),2);
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
void unimodalniMulti(double h, std::vector<double> tocka, double& l, double &r, AbstractFunction& Class, int dim)
{
	for(int i=0;i<tocka.size();i++)
	{
		if(i!=dim) tocka[i]=0;
	}
	//std::cout<<"Ulazni multi-vektor: ";
	//for(int i=0;i<tocka.size();i++) std::cout<<tocka[i]<<" "<<std::endl;
	std::vector<double> vl(tocka.size(),0.0);
	std::vector<double> vr(tocka.size(),0.0);
	vl[dim] = tocka[dim] - h, vr[dim] = tocka[dim] + h; 
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
	l = vl[dim];
	r = vr[dim];	
	return;
	}
	else if(fm > fr)
	{
		do
		{	vl = m;
			m = vr;
			fm = fr;
			vr[dim] = tocka[dim] + h * (step *= 2);
			fr = Class.function(vr);
		} while(fm > fr);
		l = vl[dim];
		r = vr[dim];
	}
	else
	{
		do
		{	vr = m;
			m = vl;
			fm = fl;
			vl[dim] = tocka[dim] - h * (step *= 2);
			fl = Class.function(vl);
			//std::cout<<"fl"<<" "<<"fm"<<" "<<"fr"<<std::endl;
			//std::cout<<fl<<" "<<fm<<" "<<fr<<std::endl;
		} while(fm > fl);
		l = vl[dim];
		r = vr[dim];
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
double Zlatni_rezMulti(bool point,double h, std::vector<double> t,std::vector<double> e, AbstractFunction& Class, int dim)
{	
	double a,b;
	if(point)
	{
		//izracunaj prvo unimodalni interval
		unimodalniMulti(h,t,a,b,Class,dim);
		//std::cout<<a<<" "<<b<<std::endl;
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
	while((b - a) > e[dim])
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
	std::vector<double> x(x0);
	std::vector<double> ref,xS;
	double temp;
	bool check=false;
	do{
		xS = x;
		for(int i=0;i<x.size();i++)
		{
			//std::cout<<i<<std::endl;
			temp=Zlatni_rezMulti(true,1,x,eps,Class,i);
			x[i]=temp;
			//for(int i=0;i<x.size();i++) std::cout<<x[i]<<std::endl;
			//for(int i=0;i<xS.size();i++) std::cout<<xS[i]<<std::endl;
		}
		ref = subtract(x,xS);
		//for(int i=0;i<ref.size();i++) std::cout<<ref[i]<<std::endl;
		absoluteValue(ref);
		check = compareVectors(ref,eps);
	}while(!check);
	return x;
}


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
	
	/*
	std::cout<<"Veličina polja: "<<array.size()<<std::endl;
	for(int i=0;i<array.size();i++)
	{
		std::cout<<i<<" ";
		temp2 = array[i];
		for(int j=0;j<temp2.size();j++)
		{
			//std::cout<<j<<" ";
			std::cout<<temp2[j]<<" ";
		}
		std::cout<<std::endl;
	}
	*/
	return array;
}

//izracunaj i vrati centroid skupa vrijednosti
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

//izracunaj refleksivnu tocku
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
	/*
	for(int i=0;i<array.size()-1;i++)
	{
		for(int j=0;j<array.size();j++)
		{
			exitVector[i]+=pow((array[j][i]-xc[i]),2);
		}
		exitVector[i]=sqrt(exitVector[i]/array.size());
	}
	return;
	*/
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
	
	/*
	std::cout<<"Velicina polja: "<<array.size()<<std::endl;
	std::vector<double> temp2;
	for(int i=0;i<array.size();i++)
	{
		std::cout<<i<<" ";
		temp2 = array[i];
		for(int j=0;j<temp2.size();j++)
		{
			//std::cout<<j<<" ";
			std::cout<<temp2[j]<<" ";
		}
		std::cout<<std::endl;
	}
	*/	
	
	
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


//Explore function for Hooks-Jeeves algorithm
std::vector<double> explore(std::vector<double> xP, std::vector<double> Dx,AbstractFunction& Class,std::map<std::vector<double>, double>& lookUpTable)
{
	std::vector<double> x=xP;
	double P,N;
	for(int i=0;i<x.size();i++)
	{
		P = getValue(x,Class,lookUpTable);//Class.function(x);
		x[i] = x[i]+Dx[i];
		N = getValue(x,Class,lookUpTable);//Class.function(x);
		if(N>P)
		{
			x[i]=x[i]-2*Dx[i];
			N = getValue(x,Class,lookUpTable);//Class.function(x);
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

std::vector<double> HookeJeeves(std::vector<double> x0, std::vector<double> precision,std::vector<double> Dx,AbstractFunction& Class)
{
	std::vector<double> xB=x0;
	std::vector<double> xP=x0;
	std::vector<double> xN,temp,tempX;
	std::map<std::vector<double>, double> lookUpTable;
	tempX=divide(precision,2);
	bool check = false;
	do
	{
		xN = explore(xP,Dx,Class,lookUpTable);
		//if(Class.function(xN)<Class.function(xB))
		if(getValue(xN,Class,lookUpTable)<getValue(xB,Class,lookUpTable))
		{
			temp = multiply(xN,2);
			xP = subtract(temp,xB);
			xB=xN;
		}
		else
		{
			Dx = divide(Dx,2);
			xP = xB;		
		}
		check = compareVectors(Dx,tempX);
	}while(!check);
	return xB;

}


void myfunction (std::string i) {  // function:
  std::cout << ' ' << i <<std::endl;
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
		std::vector<double> tocka1,tocka2,tocka4,temp;
		tocka1.push_back(-1.9);
		tocka1.push_back(2);
		tocka2.push_back(0.1);
		tocka2.push_back(0.3);
		tocka4.push_back(5.1);
		tocka4.push_back(1.1);
	
		//funkcija 1
		//Koordinatne osi
		rezultati.push_back(KoordiantneOsi(tocka1,preciznost,func1));
		brojPoziva.push_back(func1.getNumbers());
		func1.restartCount();
		std::cout<<"Koordinatne osi - 1"<<std::endl;
	
		//Nelder-Mead
		rezultati.push_back(NelderMead(tocka1,distance,alfa,beta,gamma,preciznost[0],func1));
		brojPoziva.push_back(func1.getNumbers());
		func1.restartCount();

		//Hook-Jeeves
		rezultati.push_back(HookeJeeves(tocka1,preciznost,pomaciFunkcije,func1));
		brojPoziva.push_back(func1.getNumbers());
		func1.restartCount();

		//funkcija 2
		//Koordinatne osi
		rezultati.push_back(KoordiantneOsi(tocka2,preciznost,func2));
		brojPoziva.push_back(func2.getNumbers());
		func2.restartCount();
		std::cout<<"Koordinatne osi - 2"<<std::endl;

		//Nelder-Mead
		rezultati.push_back(NelderMead(tocka2,distance,alfa,beta,gamma,preciznost[0],func2));
		brojPoziva.push_back(func2.getNumbers());
		func2.restartCount();

		//Hook-Jeeves
		rezultati.push_back(HookeJeeves(tocka2,preciznost,pomaciFunkcije,func2));
		brojPoziva.push_back(func2.getNumbers());
		func2.restartCount();
	
		//funkcija 3
		//Koordinatne osi
		rezultati.push_back(KoordiantneOsi(tocka,preciznost,func3));
		brojPoziva.push_back(func3.getNumbers());
		func3.restartCount();
		std::cout<<"Koordinatne osi - 3"<<std::endl;

		//Nelder-Mead
		rezultati.push_back(NelderMead(tocka,distance,alfa,beta,gamma,preciznost[0],func3));
		brojPoziva.push_back(func3.getNumbers());
		func3.restartCount();

		//Hook-Jeeves
		rezultati.push_back(HookeJeeves(tocka,preciznost,pomaciFunkcije,func3));
		brojPoziva.push_back(func3.getNumbers());
		func3.restartCount();

		//funkcija 4
		//Koordinatne osi
		rezultati.push_back(KoordiantneOsi(tocka4,preciznost,func4));
		brojPoziva.push_back(func4.getNumbers());
		func4.restartCount();
		std::cout<<"Koordinatne osi - 4"<<std::endl;
	
		//Nelder-Mead
		rezultati.push_back(NelderMead(tocka4,distance,alfa,beta,gamma,preciznost[0],func4));
		brojPoziva.push_back(func4.getNumbers());
		func4.restartCount();

		//Hook-Jeeves
		rezultati.push_back(HookeJeeves(tocka4,preciznost,pomaciFunkcije,func4));
		brojPoziva.push_back(func4.getNumbers());
		func4.restartCount();
	
		std::cout<<"Tablica rezultata"<<std::endl;
		std::cout<<"Nadjeni minimum/Broj poziva"<<std::endl;
		std::cout<<"--------------------------------------------------------------------------"<<std::setw(5)<<std::endl;
		for(int i=0;i<4;i++)
		{
			for(int j=0;j<3;j++)
			{
				//print Minimum
				temp = rezultati[i*3+j];
				for(int k=0;k<temp.size();k++) std::cout<<std::setw(12)<<temp[k]<<" ";
				//print number of calls
				std::cout<<brojPoziva[i*3+j];
			}
			std::cout<<std::endl;
		}
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
		std::cout<<"Nelder-Mead"<<std::endl;
		for(int i=0;i<nelderMead.size();i++)std::cout<< nelderMead[i]<<" ";
		std::cout<<brojPoziva[0];
		std::cout<<std::endl;
		std::cout<<"Hook-Jeeves"<<std::endl;
		for(int i=0;i<hookJeeves.size();i++)std::cout<< hookJeeves[i]<<" ";
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
