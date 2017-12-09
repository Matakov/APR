#include "genetski.h"
//https://stackoverflow.com/questions/22746429/c-decimal-to-binary-converting
//https://www.programiz.com/cpp-programming/examples/binary-decimal-convert

/*
Ulazni parametri:
		-referenca na vektor
		-lijeva granica
		-desna granica
		-broj bitova
Izlazni parametri:
		-referencirani vektor
*/
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
/*
Ulazni parametri:
		-referenca na vektor
		-lijeva granica
		-desna granica
		-broj bitova
Izlazni parametri:
		-referencirani vektor
*/
void TurnFromBinary(std::vector<double>& b,double borderLeft, double borderRight,double numBytes)
{
	std::vector<double> temp;
	for(int i=0;i<b.size();i++)
	{
		temp.push_back(borderLeft+b[i]/(double)(pow(2,numBytes)-1)*(borderRight-borderLeft));	
	}
	b=temp;
	return;
}
/*
Ulazni parametri:
		-referenca na vektor
		-ulazni vektor
		-broj bitova
Izlazni parametri:
		-referencirani vektor
*/
void getToBinaryString(std::vector<double>& result,std::vector<double> x,double numBytes)
{
	result.clear();
	std::vector<double> temp,temp2;
	int remainder;//i=1,step=1;
	//std::cout<<"Numbers are: ";
	//for(int i=0;i<x.size();i++) std::cout<<x[i]<<" ";
	//std::cout<<std::endl;
	for(int i=0;i<x.size();i++)
	{
		temp2.clear();
		while(x[i]!=0)
		{
			remainder=(int)x[i]%2;
			//std::cout << "Number: " << x[i] << "/2, Remainder = " << remainder << ", Quotient = " << (int)x[i]/2 << std::endl;
			x[i]=(int)x[i]/2;
			temp2.push_back(remainder);
		}
		while(temp2.size()<numBytes)
		{
			temp2.push_back(0);
		}
		std::reverse(temp2.begin(), temp2.end());
		for(int k=0;k<temp2.size();k++)
		{
			result.push_back(temp2[k]);
		}
	}

	return;
}
/*
Ulazni parametri:
		-referenca na vektor
		-ulazni vektor
		-broj bitova
Izlazni parametri:
		-referencirani vektor
*/
void getFromBinaryString(std::vector<double>& result,std::vector<double> x,double numBytes)
{
	int t=x.size()/numBytes;
	//std::cout<<t<<std::endl;
	double number;
	//std::vector<double> temp;
	for(int i=0;i<t;i++)
	{
		number=0;
		//temp.clear();
		for(int j=0;j<numBytes;j++) 	//temp[j]=x[t*numBytes+j];
		{
			//std::cout<<"N: "<<i<<" "<<j<<" "<<x[i*numBytes+j]<<" ";
			number+=x[i*numBytes+j]*pow(2,numBytes-(j+1));
		}
		//std::cout<<std::endl;
		result.push_back(number);
	}
}
/*
Ulazni parametri:
		-referenca na polje
		-lijeva granica
		-desna granica
		-broj bitova
Izlazni parametri:
		-referencirano polje
*/
void binarizeArray(std::vector<std::vector<double>>&array,double borderLeft, double borderRight,double numBytes)
{
	std::vector<std::vector<double>> temp(array);
	for(int i=0;i<temp.size();i++)
	{
		TurnToBinary(temp[i],borderLeft, borderRight,numBytes);
		getToBinaryString(array[i],temp[i],numBytes);
	}
	//array=temp;
}
/*
Ulazni parametri:
		-referenca na polje
		-lijeva granica
		-desna granica
		-broj bitova
Izlazni parametri:
		-referencirano polje
*/
void debinarizeArray(std::vector<std::vector<double>>&array,double borderLeft, double borderRight,double numBytes)
{
	std::vector<std::vector<double>> temp;
	std::vector<double> temp_vector;
	for(int i=0;i<array.size();i++)
	{
		temp_vector.clear();
		//temp[i].clear();
		getFromBinaryString(temp_vector,array[i],numBytes);
		//for(int j=0;j<temp_vector.size();j++) std::cout<<temp_vector[j]<<" ";
		//std::cout<<std::endl;
		TurnFromBinary(temp_vector,borderLeft, borderRight,numBytes);
		//TurnToBinary(temp[i],borderLeft, borderRight,numBytes);
		//getToBinaryString(temp[i],array[i],numBytes);
		//for(int j=0;j<temp_vector.size();j++) std::cout<<temp_vector[j]<<" ";
		temp.push_back(temp_vector);
	}
	array=temp;
}
/*
DRŽATI U DECIMALNOM FORMATU DO SAMOG KRIŽANJA -> TADA PREBACITI U BINARNI STRING I TO KRIŽATI I MUTIRATI!!!!!!
*/
/*
Ulazni parametri:
		-referenca na mapu vrijednosti
		-parametar a
		-parametar b
Izlazni parametri:
		-referencirani vektor
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

/*
Ulazni parametri:
		-referenca na mapu
Izlazni parametri:
		-referencirana mapa
*/
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
void evaluatePopulace(AbstractFunction& Class, std::vector<std::vector<double>>& array, std::map<std::vector<double>, double>& valueMap)
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
/*
ulazni parametri:
		- vecičina populacije
		- lijeva granica
		- desna granica
		- dimenzionalnost podataka
		- broj bitova


Izlazni parametri:
		- populacija
*/
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

bool checkIfInArray(std::vector<double> vec, double x)
{
	for(int i=0;i<vec.size();i++)
	{
		if(x==vec[i]) return false;
	}
	return true;	
}

void nTurnirSelecion(std::vector<std::vector<double>>&selectedarray,std::vector<std::vector<double>> const &array, std::map<std::vector<double>, double>& valueMap, double n=3, double mode=1)
{
	std::vector<double> temp;
	bool truth = true;
	double r;
	if(mode==1)
	{
		double leftBorder=0;
		double sumGodness=0;
		double iter=0;
		std::map<std::vector<double>, double>::iterator it;
		for(it=valueMap.begin();it!=valueMap.end();++it)
		{
			sumGodness+=it->second;
		}
		do
		{
			truth = true;
			r =(double)(int)(array.size()*((double) rand() / (RAND_MAX)));
			iter=0;
			for(it=valueMap.begin();it!=valueMap.end();++it)
			{
				iter++;				
				if((leftBorder<r)&&(r<leftBorder+it->second/sumGodness))
				{
					truth = checkIfInArray(temp,iter);
					if(truth) temp.push_back(iter);	
				}
				leftBorder+=it->second/sumGodness;
			}
		}
		while(temp.size()<n);
	}
	else
	{
		//for(int k=0;k<n;k++)
		do
		{
			r =(double)(int)(array.size()*((double) rand() / (RAND_MAX)));
			truth = checkIfInArray(temp,r);			
			if(truth) temp.push_back(r);
		}
		while(temp.size()<n);
	}
	for(int i=0;i<temp.size();i++)
	{
		selectedarray.push_back(array[i]);
	}
}


void geneticAlgorithm(AbstractFunction& Class,std::vector<double>& result,double brPop,double prob, double brEval,double borderLeft, double borderRight, double vectorSize, double numBytes=1, double n=3, double mode=1)
{
	std::vector<std::vector<double>> array;
	std::vector<std::vector<double>> selectedarray;
	std::map<std::vector<double>, double> valueMap;
	int iter=0;
	createPopulace(array,brPop,borderLeft,borderRight,vectorSize, numBytes);
	evaluatePopulace(Class,array,valueMap);
	do
	{
		//treba izbrisati i ponovno izracunati valueMap
		evaluatePopulace(Class,array,valueMap);
		iter++;
		selectedarray.clear();
		nTurnirSelecion(selectedarray,array,valueMap,n,mode);



		if(iter>=brEval) break;
	}
	while(true);
	//petlja moze stati kada funkcija cilja padne ispod 1e-6
}



int main(int argc, char* argv[])
{
	double brPopulacije,brojBitova,vjerojatnost,brEval,borderLeft, borderRight, velicinaVektora;
	int zadatak=1;
	if(argc<3)
	{
		std::cout<<"Nisi zadao ulazni txt file ili zadatak koji se treba vrtiti\nPrekidam program."<<std::endl;
		std::exit(0);
	}
	zadatak=atof(argv[2]);
	openFile(argv[1], brPopulacije, brojBitova, vjerojatnost, brEval, borderLeft, borderRight, velicinaVektora);	

	std::vector<double> number,x;
	number.push_back(15);
	number.push_back(9);
	//std::vector<double> result;
	//getToBinaryString(result ,number , 10);
	//std::cout<<"Binary numbers: ";
	//for(int i=0;i<result.size();i++) std::cout<<result[i]<<" ";
	//std::cout<<std::endl;
	//getFromBinaryString(x,result,10);
	//std::cout<<"Decimal numbers: ";
	//for(int i=0;i<x.size();i++) std::cout<<x[i]<<" ";

	std::vector<std::vector<double>> array;
	
	/*
	createPopulace(array,brPopulacije,borderLeft,borderRight,velicinaVektora,brojBitova);
	std::cout<<"Decimal numbers: ";	
	for(int i=0;i<array.size();i++)
	{
		for(int j=0;j<array[i].size();j++) std::cout<<array[i][j]<<" ";
		std::cout<<std::endl;
	}
	std::vector<std::vector<double>> array2(array);
	binarizeArray(array2,borderLeft,borderRight,brojBitova);
	for(int i=0;i<array2.size();i++)
	{
		for(int j=0;j<array2[i].size();j++) std::cout<<array2[i][j]<<" ";
		std::cout<<std::endl;
	}
	std::vector<std::vector<double>> array3(array2);
	debinarizeArray(array3,borderLeft,borderRight,brojBitova);
	std::cout<<array3.size()<<std::endl;
	for(int i=0;i<array3.size();i++)
	{
		for(int j=0;j<array3[i].size();j++) std::cout<<array3[i][j]<<" ";
		std::cout<<std::endl;
	}
	*/
	std::cout<<(double)(int)(16*((double) rand() / (RAND_MAX)))<<std::endl;
	std::cout<<(double)(int)(16*((double) rand() / (RAND_MAX)))<<std::endl;
	return 0;
}
