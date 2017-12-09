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

void binarizeVector(std::vector<double>& vector,double borderLeft, double borderRight,double numBytes)
{
	std::vector<double> temp(vector);
	TurnToBinary(temp,borderLeft, borderRight,numBytes);
	getToBinaryString(vector,temp,numBytes);
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

void debinarizeVector(std::vector<double>& vector,double borderLeft, double borderRight,double numBytes)
{
	std::vector<double> temp(vector);
	getFromBinaryString(temp,vector,numBytes);
	TurnFromBinary(temp,borderLeft, borderRight,numBytes);
	//TurnToBinary(temp,borderLeft, borderRight,numBytes);
	//getToBinaryString(vector,temp,numBytes);
	vector=temp;
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
		if(x==vec[i]) return true;
	}
	return false;	
}

/*
N turnirska selekcija
Ulazni parametri:
		- referenca na selektirano polje
		- referenca na populaciju
		- referenca na mapu vrijednosti
		- broj natjecatelja
		- nacin rada
Izlazni parametri:
		- selektirano polje
		- mapa vrijednosti
*/
void nTurnirSelecion(std::vector<std::vector<double>>& selectedarray,std::vector<std::vector<double>> const &array, std::map<std::vector<double>, double>& valueMap, double n=3, double mode=1)
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
					if(!truth) temp.push_back(iter);	
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
			if(!truth) temp.push_back(r);
		}
		while(temp.size()<n);
	}
	for(int i=0;i<temp.size();i++)
	{
		selectedarray.push_back(array[temp[i]]);
	}
}
/*
Function to eliminate worst element in a selected array
Ulazno/izlazni parametri:
		-referenca na selektirano polje
		-referenca na populaciju
		-referenca na mapu vrijednosti polja
*/
void eliminateWorst(std::vector<std::vector<double>>& selectedarray,std::vector<std::vector<double>>& array,  std::map<std::vector<double>, double>& valueMap)
{
	double maxIndex=0;
	double max=0;
	for(int i=0;i<selectedarray.size();i++)
	{
		if(max<valueMap[selectedarray[i]])
		{
			maxIndex=i;
			max=valueMap[selectedarray[i]];
		}		
	}
	valueMap.erase(selectedarray[maxIndex]);
	array.erase(std::remove(array.begin(), array.end(), selectedarray[maxIndex]), array.end());	
	return;
}

/*
Pronadji najbolja 2 u selektiranom polju i vrati ih
*/
//return max index and max2 index
void getMaximumIndex(std::vector<std::vector<double>> selectedArray,std::map<std::vector<double>, double>& lookUpTable, double& h, double& h2)
{
	//lookUpTable.clear();
	h=0;
	double left;
	std::vector<double> values;
	for(std::map<std::vector<double>,double>::iterator it = lookUpTable.begin(); it != lookUpTable.end(); ++it) {
		values.push_back(it->second);
	}
	std::sort(values.begin(),values.end());
	double max = values[0];
	double max2 = values[1];
	for(int i=0;i<selectedArray.size();i++)
	{
		//left = getValue(array[i],Class,lookUpTable);
		left = lookUpTable[selectedArray[i]];
		if(left==max)  h = i;
		if(left==max2) h2 = i;
	}
	return;
}

//Funkcija koja uz ulazni parametar od bitnovne vrijenosti generira nasumični vektor
std::vector<double> generateR(double numBytes)
{
	std::vector<double> temp;	
	double r=(double) rand() / (RAND_MAX);
	for(int i=0;i<numBytes;i++)
	{
		r=(double) rand() / (RAND_MAX);
		if(r<0.5) temp.push_back(0);
		else temp.push_back(1);
	}
	return temp;
}

/*
Funkcija koja radi križanje i mutaciju jedinki
Ulazni parametri:
		-referenca dijeteta
		-referenca setektiranih jedinki
		-mapa vrijednosti
		-vrsta križanja
		-da li je aritmetički ili binarno
		-granice vrijednosti
*/
void crossover(std::vector<double>& child,std::vector<std::vector<double>>& selectedarray,std::map<std::vector<double>, double>& valueMap, double typeOfCrossover, double numBytes, double borderLeft, double borderRight)
{
	child.clear();
	double r;
	double best1,best2;
	std::vector<double> R;
	if(numBytes>1) //Binarni prikaz
	{
		//Binarize array
		getMaximumIndex(selectedarray,valueMap,best1,best2);
		binarizeArray(selectedarray,borderLeft,borderRight,numBytes);
		//one point crossover
		if(typeOfCrossover==1)
		{
			r =(double)(int)(numBytes*((double) rand() / (RAND_MAX)));
			for(int i=0;i<(selectedarray[0].size()/numBytes);i++)
			{
				for(int j=0;j<numBytes;j++)
				{
					if(j<r)
					{
						child.push_back(selectedarray[best1][i*numBytes+j]);
					}
					else
					{
						child.push_back(selectedarray[best1][i*numBytes+j]);					
					}
				}
			}
			//array.push_back(child);
		}
		//uniform crossover
		if(typeOfCrossover==2)
		{
			r =(double)(int)(numBytes*((double) rand() / (RAND_MAX)));
			for(int i=0;i<(selectedarray[0].size()/numBytes);i++)
			{
				R = generateR(numBytes);
				for(int j=0;j<numBytes;j++)
				{

					child.push_back(selectedarray[best1][i*numBytes+j]*selectedarray[best2][i*numBytes+j]+R[j]*((1-selectedarray[best1][i*numBytes+j])*selectedarray[best2][i*numBytes+j]+selectedarray[best1][i*numBytes+j]*(1-selectedarray[best2][i*numBytes+j])));

				}
			}
			//array.push_back(child);
		}
		//debinarize array
		debinarizeArray(selectedarray,borderLeft,borderRight,numBytes);
	}
	else	//Aritmeticki prikaz
	{
		getMaximumIndex(selectedarray,valueMap,best1,best2);
		//one point crossover
		if(typeOfCrossover==1)
		{
			r =(double)(int)(selectedarray[0].size()*((double) rand() / (RAND_MAX)));
			for(int j=0;j<selectedarray[0].size();j++)
			{
				if(j<r)
				{
					child.push_back(selectedarray[best1][j]);
				}
				else
				{
					child.push_back(selectedarray[best2][j]);
				}
			}
			//array.push_back(child);
		}
		//Arithmetic Recombination
		if(typeOfCrossover==1)
		{
			r =(double) rand() / (RAND_MAX);
			for(int j=0;j<selectedarray[0].size();j++)
			{
				child.push_back(r*selectedarray[best1][j]+(1-r)*selectedarray[best2][j]);

			}
			//array.push_back(child);
		}
	}	
	return;
}

/*
Funkcija koja radi križanje i mutaciju jedinki
Ulazni parametri:
		-referenca setektiranih jedinki
		-referenca na populaciju
		-vrsta križanja
		-da li je aritmetički ili binarno
*/
void mutation(std::vector<double>& child, double probability, double typeOfMutation, double numBytes, double borderLeft, double borderRight)
{
	double temp,l,r=(double)(int)(numBytes*((double) rand() / (RAND_MAX)));
	if(numBytes>1) //Binarna mutacija
	{
		binarizeVector(child, borderLeft, borderRight, numBytes);
		if(typeOfMutation==1) //Bit Flip Mutation
		{
			for(int i=0;i<child.size()/numBytes;i++)
			{
				r=(double)(int)(numBytes*((double) rand() / (RAND_MAX)));
				child[i*numBytes+r]=1-child[i*numBytes+r];
			}
		}
		if(typeOfMutation==2) //Swap Mutation
		{
			for(int i=0;i<child.size()/numBytes;i++)
			{
				r=(double)(int)(numBytes*((double) rand() / (RAND_MAX)));
				l=(double)(int)(numBytes*((double) rand() / (RAND_MAX)));
				if(l>r)
				{
					temp=l;
					l=r;
					r=temp;
				}
				temp = child[i*numBytes+l]; 
				child[i*numBytes+l]=child[i*numBytes+r];
				child[i*numBytes+r]=child[i*numBytes+l];
			}
		}
		if(typeOfMutation==3) //Inversion Mutation
		{
			for(int i=0;i<child.size()/numBytes;i++)
			{
				r=(double)(int)(numBytes*((double) rand() / (RAND_MAX)));
				l=(double)(int)(numBytes*((double) rand() / (RAND_MAX)));
				if(l>r)
				{
					temp=l;
					l=r;
					r=temp;
				}
				std::reverse(child.begin()+i*numBytes+l,child.begin()+i*numBytes+r);
			}
		}
		if(typeOfMutation==4) //Random Shiffle
		{
			for(int i=0;i<child.size()/numBytes;i++)
			{
				r=(double)(int)(numBytes*((double) rand() / (RAND_MAX)));
				l=(double)(int)(numBytes*((double) rand() / (RAND_MAX)));
				if(l>r)
				{
					temp=l;
					l=r;
					r=temp;
				}
				std::random_shuffle(child.begin()+i*numBytes+l,child.begin()+i*numBytes+r);
			}
		}
		if(typeOfMutation==5) //Boundary mutation
		{
			for(int i=0;i<child.size()/numBytes;i++)
			{
				r=(double) rand() / (RAND_MAX);
				if(r<0.5)
				{
					for(int j=0;j<numBytes;j++)
					{
						child[i*numBytes+j]=0;
					}
				}
				else
				{
					for(int j=0;j<numBytes;j++)
					{
						child[i*numBytes+j]=1;
					}
				}
			}
			
		}
		debinarizeVector(child, borderLeft, borderRight, numBytes);
	}
	else //Mutacija sa pomičnim korakom
	{
		if(typeOfMutation==1) //Uniform mutation
		{
			for(int i=0;i<child.size();i++)
			{
				child[i]=((borderRight-borderLeft)*((double) rand() / (RAND_MAX)));
			}			
		}
		if(typeOfMutation==3) //Gaussian mutation
		{
			for(int i=0;i<child.size();i++)
			{
				std::default_random_engine generator;
  				std::normal_distribution<double> distribution(borderRight-borderLeft,3);
				r = distribution(generator);
				if(r>borderRight) r=borderRight;
				if(r<borderLeft) r=borderLeft;	
				child[i]=r;
			}			
		}
		if(typeOfMutation==3) //Boundary mutation
		{
			for(int i=0;i<child.size();i++)
			{
				r=(double) rand() / (RAND_MAX);
				if(r<0.5)
				{
					child[i]=borderLeft;
				}
				else
				{
					child[i]=borderRight;
				}
			}			
		}
	}
	return;
}

void geneticAlgorithm(AbstractFunction& Class,std::vector<double>& result,double brPop,double prob, double brEval,double borderLeft, double borderRight, double vectorSize, double numBytes=1, double n=3, double mode=1,double typeOfCrossover=1, double typeOfMutation=1, double probability=0.2)
{
	std::vector<std::vector<double>> array;
	std::vector<std::vector<double>> selectedarray;
	std::vector<double> child;
	std::map<std::vector<double>, double> valueMap;
	int iter=0;
	double best,best2;
	createPopulace(array,brPop,borderLeft,borderRight,vectorSize, numBytes);
	evaluatePopulace(Class,array,valueMap);
	do
	{
		//valueMap.clear();
		//treba izbrisati i ponovno izracunati valueMap
		evaluatePopulace(Class,array,valueMap);
		iter++;
		selectedarray.clear();
		child.clear();
		nTurnirSelecion(selectedarray,array,valueMap,n,mode);
		eliminateWorst(selectedarray,array,valueMap);
		//crossoverNMutation(selectedarray,array,valueMap,typeOfCrossover, numBytes, probability, borderLeft, borderRight);
		crossover(child,selectedarray,valueMap, typeOfCrossover, numBytes, borderLeft, borderRight);
		mutation(child, probability, typeOfMutation, numBytes, borderLeft, borderRight);
		array.push_back(child);
		//Treba implementirati elitizam
		if(Class.function(child)<1e-6) break;
		if(iter>=brEval) break;
	}
	while(true);
	//petlja moze stati kada funkcija cilja padne ispod 1e-6
	getMaximumIndex(array, valueMap,best,best2);
	result=array[best];
	return;
}


int main(int argc, char* argv[])
{
	double brPopulacije,brojBitova,vjerojatnost,brEval,borderLeft, borderRight, velicinaVektora, vrstaMutacija, vrstaKrizanja;
	int zadatak=1;
	if(argc<3)
	{
		std::cout<<"Nisi zadao ulazni txt file ili zadatak koji se treba vrtiti\nPrekidam program."<<std::endl;
		std::exit(0);
	}
	zadatak=atof(argv[2]);
	openFile(argv[1], brPopulacije, brojBitova, vjerojatnost, brEval, borderLeft, borderRight, velicinaVektora, vrstaKrizanja, vrstaMutacija);	

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
