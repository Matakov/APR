#import "lab2.h"

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
	double lambda0 = 0;
	//std::cout<<"Ulazni multi-vektor: ";
	//for(int i=0;i<tocka.size();i++) std::cout<<tocka[i]<<" "<<std::endl;
	std::vector<double> vl(tocka.size(),0.0);
	std::vector<double> vr(tocka.size(),0.0);
	double lambdaL = lambda0 - h;
	double lambdaR = lambda0 + h;
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
double Zlatni_rezMulti(double h, std::vector<double> tocka,double e, AbstractFunction& Class, std::vector<double> multiToOne)
{	
	double lambda0 = 0;
	//std::vector<double> a,b;
	double a,b;
	//izracunaj prvo unimodalni interval
	unimodalniMulti(h,tocka,a,b,Class,multiToOne);
	//std::cout<<a<<" "<<b<<std::endl;
	
	double k = 0.5*(sqrt(5)-1);
	std::vector<double> vc(tocka.size(),0.0);
	double c = b - k * (b - a);
	//vc[dim]=c;
	for (int i=0;i<tocka.size();i++) vc[i] = multiToOne[i]*c;
	std::vector<double> vd(tocka.size(),0.0);
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


