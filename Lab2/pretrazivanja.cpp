#include<iostream>
#include<iomanip>
#include<sstream>

/*
Postupak trazenja unimodalnog intervala

Ulazne velicine:
- tocka: pocetna tocka pretrazivanja
- h: pomak pretrazivanja
- f: ciljna funkcija

Izlazne vrijednosti:
- unimodalni interval [l, r]
*/
double unimodalni(double h, double tocka, double& l, double &r)
{
	l = tocka - h, r = tocka + h; 
	double m = tocka;
	double fl, fm, fr;
	int step = 1;

	double fm = f(tocka);
	double fl = f(l);
	double fr = f(r);

	if(fm < fr && fm < fl)
		return;
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
*/
double Zlatni_rez(double a,double b,double e)
{	
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


//Zlatni rez sa zadanom pocetnom tockom

double Zlatni_rez(double t,double e)
{	
	double a,b;
	//izracunaj prvo unimodalni interval
	unimodalni(h,t,a,b);
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






int main(){



	return 0;
}
