#include "genetski.h"
//https://stackoverflow.com/questions/22746429/c-decimal-to-binary-converting
//https://www.programiz.com/cpp-programming/examples/binary-decimal-convert
void TurnToBinary(std::vector<double>& x,double borderLeft, double borderRight,double numBytes)
{
	std::vector<double> temp;
	for(int i=0;i<x.size();i++)
	{
		temp.push_back((pow(2,numBytes)-1)*(x[i]-borderLeft)/(borderRight-borderLeft));	
	}
	x=temp;
	return;
}

void TurnFromBinary(std::vector<double>& b,double borderLeft, double borderRight,double numBytes)
{
	std::vector<double> temp,x;
	for(int i=0;i<x.size();i++)
	{
		temp.push_back(borderLeft+b[i]/(pow(2,numBytes)-1)*(borderRight-borderLeft));	
	}
	b=temp;
	return;
}
/*
DRŽATI U DECIMALNOM FORMATU DO SAMOG KRIŽANJA -> TADA PREBACITI U BINARNI STRING I TO KRIŽATI I MUTIRATI!!!!!!
*/
void windowing(std::map<std::vector<double>, double>& valueMap,double a, double b)
{
	std::vector<double> values;
	for(std::map<std::vector<double>,double>::iterator it = valueMap.begin(); it != valueMap.end(); ++it) {
		values.push_back(it->second);
	}
	std::sort(values.begin(),values.end());
	double max = values[values.size()-1];
	//double max2 = values[values.size()-2];
	double min = values[0];
	/*	
	for(int i=0;i<array.size();i++)
	{
		left = getValue(array[i],Class,valueMap);
		if(left==max)  h = i;
		if(left==max2) h2 = i;
	}
	*/
	for(std::map<std::vector<double>,double>::iterator it = valueMap.begin(); it != valueMap.end(); ++it)
	{
		valueMap[it->first]=a+(b-a)*(valueMap[it->first]-min)/(max-min);
		//it->second;
	}
	return;
}


void translacijaDobrote(std::map<std::vector<double>, double>& valueMap)
{
	std::vector<double> values;
	for(std::map<std::vector<double>,double>::iterator it = valueMap.begin(); it != valueMap.end(); ++it) {
		values.push_back(it->second);
	}
	std::sort(values.begin(),values.end());
	//double max = values[values.size()-1];
	//double max2 = values[values.size()-2];
	double min = values[0];
	/*	
	for(int i=0;i<array.size();i++)
	{
		left = getValue(array[i],Class,valueMap);
		if(left==max)  h = i;
		if(left==max2) h2 = i;
	}
	*/
	for(std::map<std::vector<double>,double>::iterator it = valueMap.begin(); it != valueMap.end(); ++it)
	{
		valueMap[it->first]=valueMap[it->first]-min;
		//it->second;
	}
	return;
}

//Dohvati vrijednost
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

//Funkcija za evaluaciju populacije
void evaluatePopulace(AbstractFunction& Class, std::vector<std::vector<double>>& array, std::map<std::vector<double>, double>& valueMap, double numBytes,double borderLeft, double borderRight)
{
	std::vector<double> temp;
	double tempValue;
	for(int i=0;i<array.size();i++)
	{
		temp=array[i];
		tempValue=getValue(array[i],Class,valueMap);
	}
	return;
}


//Funkcija za kreiranje populacije
void createPopulace(std::vector<std::vector<double>>& array,double brPop,double borderLeft, double borderRight, double vectorSize, double numBytes)
{
	double r=0;
	std::vector<double> temp(vectorSize,0.0);
	for(int i=0;i<brPop;i++)
	{
		for(int j=0;j<vectorSize;j++)
		{
			r = ((double) rand() / (RAND_MAX));
			temp[j]=borderLeft+r*(borderRight-borderLeft);
		}
		array.push_back(temp);
	}
	return;
}

void geneticAlgorithm(AbstractFunction& Class,std::vector<double>& result,double brPop,double prob, double brEval,double borderLeft, double borderRight, double vectorSize, double numBytes=1)
{
	std::vector<std::vector<double>> array;
	std::map<std::vector<double>, double> valueMap;
	int iter=0;
	createPopulace(array,brPop,borderLeft,borderRight,vectorSize, numBytes);
	evaluatePopulace(Class,array,valueMap, numBytes,borderLeft,borderRight);
	do
	{
		iter++;
		



		if(iter>=brEval) break;
	}
	while(true);
	//petlja moze stati kada funkcija cilja padne ispod 1e-6
}



int main(int argc, char* argv[])
{
	double brPopulacije,preciznost,vjerojatnost,brEval,borderLeft, borderRight, velicinaVektora;
	int zadatak=1;
	if(argc<3)
	{
		std::cout<<"Nisi zadao ulazni txt file ili zadatak koji se treba vrtiti\nPrekidam program."<<std::endl;
		std::exit(0);
	}
	zadatak=atof(argv[2]);
	openFile(argv[1], brPopulacije, preciznost, vjerojatnost, brEval, borderLeft, borderRight, velicinaVektora);	

	return 0;
}
