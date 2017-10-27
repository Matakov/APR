#include<iostream>
#include<iomanip>
#include<sstream>
#include<functional>
#include<fstream>
#include <algorithm>
#include <vector>
#include <stdlib.h>
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

/*
Postupak trazenja unimodalnog intervala

Ulazne velicine:
- tocka: pocetna tocka pretrazivanja
- h: pomak pretrazivanja
- f: ciljna funkcija

Izlazne vrijednosti:
- unimodalni interval [l, r]
*/
void unimodalni(double h, double tocka, double& l, double &r, std::function<double (double)> f)
{
	l = tocka - h, r = tocka + h; 
	double m = tocka;
	double fl, fm, fr;
	int step = 1;

	fm = f(tocka);
	fl = f(l);
	fr = f(r);

	if(fm < fr && fm < fl)	return;
	else if(fm > fr)
		do
		{	l = m;
			m = r;
			fm = fr;
			r = tocka + h * (step *= 2);
			fr = f(r);
		} while(fm > fr);
	else 
		do
		{	r = m;
			m = l;
			fm = fl;
			l = tocka - h * (step *= 2);
			fl = f(l);
		} while(fm > fl);
}

/*
Algoritam zlatnog reza

ulazne velicine:
- a, b: pocetne granice unimodalnog intervala
- e: preciznost

double Zlatni_rez(double a,double b,double e, std::function<double (double)> f)
{	
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
//*/

//Rosenbrockova 'banana' funkcija
auto func1 =[](double a, double b) -> double {return 100*pow(b-pow(a,2),2)+pow((1-a),2);};

auto func2 =[](double a, double b) -> double {return pow(a-4,2)+ 4*pow(b-2);};

auto func3 =[](double a, double b) -> double {return pow((a-b),2);};

//Jakobovićeva funkcija
auto func4 =[](double a, double b) -> double {return abs((a-b)*(a+b))+sqrt(pow(a,2)+pow(b,2));};

//Schaffer's function
auto func5 =[](std::vector<double> container) -> double {return 0.5+pow(sin(sqrt()),2);};

void myfunction (std::string i) {  // function:
  std::cout << ' ' << i <<std::endl;
}





int main(int argc, char* argv[]){
	std::ifstream myfile;
	myfile.open(argv[1]);
	std::string line;
	double preciznost;
	double tocka,a,b;
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
	std::cout<<func3(tocka)<<std::endl;
	double i,j,h=1;
	unimodalni(h,tocka,i,j,func3);
	std::cout<<i<<" "<<j<<std::endl;
	std::cout<<Zlatni_rez(h,tocka,preciznost,func3)<<std::endl;
	return 0;
}
