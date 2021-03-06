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
void TurnToBinary(std::vector<double>& x,double borderLeft, double borderRight,double numBites)
{
	std::vector<double> temp;
	for(int i=0;i<x.size();i++)
	{
		temp.push_back((pow(2,numBites)-1)*(x[i]-borderLeft)/(borderRight-borderLeft));	
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
void TurnFromBinary(std::vector<double>& b,double borderLeft, double borderRight,double numBites)
{
	std::vector<double> temp;
	for(int i=0;i<b.size();i++)
	{
		temp.push_back(borderLeft+b[i]/(double)(pow(2,numBites)-1)*(borderRight-borderLeft));	
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
void getToBinaryString(std::vector<double>& result,std::vector<double> x,double numBites)
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
		while(temp2.size()<numBites)
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
void getFromBinaryString(std::vector<double>& result,std::vector<double> x,double numBites)
{
	int t=x.size()/numBites;
	//std::cout<<t<<std::endl;
	double number;
	//std::vector<double> temp;
	for(int i=0;i<t;i++)
	{
		number=0;
		//temp.clear();
		for(int j=0;j<numBites;j++) 	//temp[j]=x[t*numBites+j];
		{
			//std::cout<<"N: "<<i<<" "<<j<<" "<<x[i*numBites+j]<<" ";
			number+=x[i*numBites+j]*pow(2,numBites-(j+1));
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
void binarizeArray(std::vector<std::vector<double>>&array,double borderLeft, double borderRight,double numBites)
{
	std::vector<std::vector<double>> temp(array);
	for(int i=0;i<temp.size();i++)
	{
		TurnToBinary(temp[i],borderLeft, borderRight,numBites);
		getToBinaryString(array[i],temp[i],numBites);
	}
	//array=temp;
}

void binarizeVector(std::vector<double>& vector,double borderLeft, double borderRight,double numBites)
{
	std::vector<double> temp(vector);
	TurnToBinary(temp,borderLeft, borderRight,numBites);
	getToBinaryString(vector,temp,numBites);
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
void debinarizeArray(std::vector<std::vector<double>>&array,double borderLeft, double borderRight,double numBites)
{
	std::vector<std::vector<double>> temp;
	std::vector<double> temp_vector;
	for(int i=0;i<array.size();i++)
	{
		temp_vector.clear();
		//temp[i].clear();
		getFromBinaryString(temp_vector,array[i],numBites);
		//for(int j=0;j<temp_vector.size();j++) std::cout<<temp_vector[j]<<" ";
		//std::cout<<std::endl;
		TurnFromBinary(temp_vector,borderLeft, borderRight,numBites);
		//TurnToBinary(temp[i],borderLeft, borderRight,numBites);
		//getToBinaryString(temp[i],array[i],numBites);
		//for(int j=0;j<temp_vector.size();j++) std::cout<<temp_vector[j]<<" ";
		temp.push_back(temp_vector);
	}
	array=temp;
}

void debinarizeVector(std::vector<double>& vector,double borderLeft, double borderRight,double numBites)
{
	std::vector<double> temp;
	getFromBinaryString(temp,vector,numBites);
	TurnFromBinary(temp,borderLeft, borderRight,numBites);
	//TurnToBinary(temp,borderLeft, borderRight,numBites);
	//getToBinaryString(vector,temp,numBites);
	
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
double getValue(std::vector<double> input, AbstractFunction& Class,std::map<std::vector<double>, double>& lookUpTable, double problem)
{
	if(lookUpTable.count(input))
	{
		return lookUpTable[input];
	}
	else
	{
		if(problem==1) lookUpTable[input]=-Class.function(input);
		else lookUpTable[input]=Class.function(input);
		return lookUpTable[input];
	}
}

//Funkcija za evaluaciju populacije
void evaluatePopulace(AbstractFunction& Class, std::vector<std::vector<double>>& array, std::map<std::vector<double>, double>& valueMap, double problem)
{
	std::vector<double> temp;
	double tempValue;
	for(int i=0;i<array.size();i++)
	{
		temp=array[i];
		tempValue=getValue(array[i],Class,valueMap,problem);
	}
	return;
}


//Funkcija za kreiranje populacije
/*
ulazni parametri:
		- veličina populacije
		- lijeva granica
		- desna granica
		- dimenzionalnost podataka
		- broj bitova


Izlazni parametri:
		- populacija
*/
void createPopulace(std::vector<std::vector<double>>& array,double brPop,double borderLeft, double borderRight, double vectorSize, double numBites)
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
	//std::cout<<n<<std::endl;
	selectedarray.clear();
	if(mode==1) //N tournament with roulette wheel
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
			leftBorder=0;
			r =(double) rand() / (RAND_MAX);
			iter=0;
			for(it=valueMap.begin();it!=valueMap.end();++it)
			{
								
				if((leftBorder<r)&&(r<leftBorder+it->second/sumGodness))
				{
					truth = checkIfInArray(temp,iter);
					if(!truth) temp.push_back(iter);	
				}
				leftBorder+=it->second/sumGodness;
				iter++;
			}
		}
		while(temp.size()<n);
	}
	else	// select random N for tournament
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
	//std::cout<<temp.size()<<std::endl;
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
void eliminateWorst(std::vector<double>& worst,std::vector<std::vector<double>>& array,  std::map<std::vector<double>, double>& valueMap)
{
	valueMap.erase(worst);
	array.erase(std::remove(array.begin(), array.end(), worst), array.end());
	//selectedarray.erase(selectedarray.begin() + maxIndex);	
	return;
}

//get best chromosome
void findBest(std::vector<double>& best,std::map<std::vector<double>, double>& valueMap)
{
	std::map<std::vector<double>, double>::iterator it;
	double bestValue=std::numeric_limits<int>::min();
	for(it=valueMap.begin();it!=valueMap.end();++it)
	{
		if(it->second>bestValue) bestValue=it->second;
	}
	for(it=valueMap.begin();it!=valueMap.end();++it)
	{
		if(it->second==bestValue) best=it->first;
	}	
	return;
}


/*
Select parents from a selected set
*/
void selectParents(std::vector<std::vector<double>>& selectedArray,std::vector<double>& parent1,std::vector<double>& parent2,std::vector<double>& worst,std::map<std::vector<double>, double>& valueMap)
{
	std::vector<double> values;
	double best,best2,outcast;
	parent1.clear();
	parent2.clear();
	worst.clear();

	for(std::map<std::vector<double>,double>::iterator it = valueMap.begin(); it != valueMap.end(); ++it)
	{
		for(int i=0;i<selectedArray.size();i++)
		{
			if(it->first==selectedArray[i]) values.push_back(it->second);
		}
	}
	std::sort(values.begin(),values.end());
	//best=values[0];
	//best2=values[1];
	//outcast=values[values.size()-1];
	best=values[values.size()-1];
	best2=values[values.size()-2];
	outcast=values[0];
	
	for(std::map<std::vector<double>,double>::iterator it = valueMap.begin(); it != valueMap.end(); ++it)
	{
		if(it->second==best) parent1=it->first;
		if(it->second==best2) parent2=it->first;
		if(it->second==outcast) worst=it->first;
	}

	return;
}


//Funkcija koja uz ulazni parametar od bitnovne vrijenosti generira nasumični vektor
std::vector<double> generateR(double numBites)
{
	std::vector<double> temp;	
	double r=(double) rand() / (RAND_MAX);
	for(int i=0;i<numBites;i++)
	{
		r=(double) rand() / (RAND_MAX);
		if(r<0.5) temp.push_back(0);
		else temp.push_back(1);
	}
	return temp;
}

void printVector(std::vector<double> vector)
{
	for(int i=0;i<vector.size();i++) std::cout<<vector[i]<<" ";
	std::cout<<std::endl;
}


/*
Crossover of binary functions

*/
void crossoverBinary(std::vector<double>& child,std::vector<double>& parent1,std::vector<double>& parent2, double typeOfCrossover, double numBites, double borderLeft, double borderRight)
{
	binarizeVector(parent1, borderLeft, borderRight, numBites);	
	binarizeVector(parent2, borderLeft, borderRight, numBites);
	std::vector<double> R;
	double r;
	int num;
	if(typeOfCrossover==0) typeOfCrossover=(double)(int)(2*((double) rand() / (RAND_MAX)))+1;
	if(typeOfCrossover==1)
	{
		num = parent1.size()/numBites;
		//std::cout<<num<<std::endl;
		for(int i=0;i<num;i++)
		{
			r =(double)(int)(numBites*((double) rand() / (RAND_MAX)));
			for(int j=0;j<numBites;j++)
			{
				if(j<r)
				{
					child.push_back(parent1[i*numBites+j]);
				}
				else
				{
					child.push_back(parent2[i*numBites+j]);					
				}
			}
		}
		//array.push_back(child);
	}
	//uniform crossover
	if(typeOfCrossover==2)
	{
		r =(double)(int)(numBites*((double) rand() / (RAND_MAX)));
		for(int i=0;i<(parent1.size()/numBites);i++)
		{
			R = generateR(numBites);
			for(int j=0;j<numBites;j++)
			{
				child.push_back(parent1[i*numBites+j]*parent2[i*numBites+j]+R[j]*((1-parent1[i*numBites+j])*parent2[i*numBites+j]+parent1[i*numBites+j]*(1-parent2[i*numBites+j])));
			}
		}
		//array.push_back(child);
	}
	//std::cout<<"Child: ";
	//printVector(child);
	//debinarize array
	debinarizeVector(child,borderLeft,borderRight,numBites);
	//debinarizeArray(selectedarray,borderLeft,borderRight,numBites);
	//std::cout<<"Child: ";
	//printVector(child);
	return;
}

void crossoverArithmetic(std::vector<double>& child,std::vector<double>& parent1,std::vector<double>& parent2, double typeOfCrossover, double numBites, double borderLeft, double borderRight)
{
	double r;
	if(typeOfCrossover==0) typeOfCrossover=(double)(int)(2*((double) rand() / (RAND_MAX)))+1;
	if(typeOfCrossover==1) //one point crossover
	{
		r =(double)(int)(parent1.size()*((double) rand() / (RAND_MAX)));
		for(int j=0;j<parent1.size();j++)
		{
			if(j<r)
			{
				child.push_back(parent1[j]);
			}
			else
			{
				child.push_back(parent2[j]);
			}
		}
		//array.push_back(child);
	}
	if(typeOfCrossover==2) //Arithmetic Recombination
	{
		r =(double) rand() / (RAND_MAX);
		for(int j=0;j<parent1.size();j++)
		{
			child.push_back(r*parent1[j]+(1-r)*parent2[j]);
		}
		//array.push_back(child);
	}
	return;
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
void crossover(std::vector<double>& child,std::vector<double>& parent1,std::vector<double>& parent2, double typeOfCrossover, double numBites, double borderLeft, double borderRight)
{
	child.clear();
	if(numBites>1) //Binarni prikaz
	{
		//Binarize array
		//binarizeArray(selectedarray,borderLeft,borderRight,numBites);
		crossoverBinary(child,parent1,parent2,typeOfCrossover,numBites,borderLeft,borderRight);
		//std::cout<<best1<<" "<<best2<<std::endl;
		//printPopulace(selectedarray);
		//one point crossover
		
		
	}
	else	//Aritmeticki prikaz
	{
		crossoverArithmetic(child,parent1,parent2,typeOfCrossover,numBites,borderLeft,borderRight);
	}	
	return;
}


/*
Mutation of a binary child
*/
void mutationBinary(std::vector<double>& child, double probability, double typeOfMutation, double numBites, double borderLeft, double borderRight)
{
	double temp,l,r=(double)(int)(numBites*((double) rand() / (RAND_MAX)));
	if(typeOfMutation==0) typeOfMutation=(double)(int)(4*((double) rand() / (RAND_MAX)))+1;
	binarizeVector(child, borderLeft, borderRight, numBites);
	if(typeOfMutation==1) //Bit Flip Mutation
	{
		for(int i=0;i<child.size()/numBites;i++)
		{
			r=(double)(int)(numBites*((double) rand() / (RAND_MAX)));
			child[i*numBites+r]=1-child[i*numBites+r];
		}
	}
	if(typeOfMutation==2) //Swap Mutation
	{
		for(int i=0;i<child.size()/numBites;i++)
		{
			r=(double)(int)(numBites*((double) rand() / (RAND_MAX)));
			l=(double)(int)(numBites*((double) rand() / (RAND_MAX)));
			if(l>r)
			{
				temp=l;
				l=r;
				r=temp;
			}
			temp = child[i*numBites+l]; 
			child[i*numBites+l]=child[i*numBites+r];
			child[i*numBites+r]=child[i*numBites+l];
		}
	}
	if(typeOfMutation==3) //Inversion Mutation
	{
		for(int i=0;i<child.size()/numBites;i++)
		{
			r=(double)(int)(numBites*((double) rand() / (RAND_MAX)));
			l=(double)(int)(numBites*((double) rand() / (RAND_MAX)));
			if(l>r)
			{
				temp=l;
				l=r;
				r=temp;
			}
			std::reverse(child.begin()+i*numBites+l,child.begin()+i*numBites+r);
		}
	}
	if(typeOfMutation==4) //Random Shiffle
	{
		for(int i=0;i<child.size()/numBites;i++)
		{
			r=(double)(int)(numBites*((double) rand() / (RAND_MAX)));
			l=(double)(int)(numBites*((double) rand() / (RAND_MAX)));
			if(l>r)
			{
				temp=l;
				l=r;
				r=temp;
			}
			std::random_shuffle(child.begin()+i*numBites+l,child.begin()+i*numBites+r);
		}
	}
	if(typeOfMutation==5) //Boundary mutation
	{
		for(int i=0;i<child.size()/numBites;i++)
		{
			r=(double) rand() / (RAND_MAX);
			if(r<0.5)
			{
				for(int j=0;j<numBites;j++)
				{
					child[i*numBites+j]=0;
				}
			}
			else
			{
				for(int j=0;j<numBites;j++)
				{
					child[i*numBites+j]=1;
				}
			}
		}
		
	}
	debinarizeVector(child, borderLeft, borderRight, numBites);	
	return;
}


/*
Mutation of arithmetic unit
*/
void mutationArithmetic(std::vector<double>& child, double probability, double typeOfMutation, double numBites, double borderLeft, double borderRight)
{
	double temp,l,r=(double)(int)(numBites*((double) rand() / (RAND_MAX)));
	if(typeOfMutation==0) typeOfMutation=(double)(int)(3*((double) rand() / (RAND_MAX)))+1;
	if(typeOfMutation==1) //Uniform mutation
	{
		for(int i=0;i<child.size();i++)
		{
			child[i]=((borderRight-borderLeft)*((double) rand() / (RAND_MAX))-borderLeft);
		}			
	}
	if(typeOfMutation==2) //Gaussian mutation
	{
		for(int i=0;i<child.size();i++)
		{
			std::default_random_engine generator;
  			std::normal_distribution<double> distribution((borderRight-borderLeft)/2,3);
			r = distribution(generator);
			if(r>borderRight) r=borderRight;
			if(r<borderLeft) r=borderLeft;	
			child[i]=r;
		}			
	}
	/*
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
	*/
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
void mutation(std::vector<double>& child, double probability, double typeOfMutation, double numBites, double borderLeft, double borderRight)
{
	double mute=(double) rand() / (RAND_MAX);
	if(numBites>1) //Binarna mutacija
	{
		//if(mute>=probability) return;
		if(mute>=(1-pow(1-probability,numBites))) return;
		mutationBinary(child,probability,typeOfMutation,numBites,borderLeft,borderRight);
	}
	else //Mutacija sa pomičnim korakom
	{
		mutationArithmetic(child,probability,typeOfMutation,numBites,borderLeft,borderRight);
	}
	return;
}
void printPopulace(std::vector<std::vector<double>> array)
{
	for(int i=0;i<array.size();i++)
	{
		for(int j=0;j<array[i].size();j++) std::cout<<array[i][j]<<" ";
		std::cout<<std::endl;
	}
	
}

void printPopulaceFitness(std::map<std::vector<double>, double>& valueMap)
{
	std::map<std::vector<double>, double>::iterator it;
	for(it=valueMap.begin();it!=valueMap.end();++it)
	{
		for(int j=0;j<it->first.size();j++) std::cout<<it->first[j]<<" ";
		std::cout<<" : "<<it->second;
		std::cout<<std::endl;
	}
	
}

void removeDuplicates(std::vector<std::vector<double>>& array)
{
	std::vector<std::vector<double>>::iterator it;
	std::sort (array.begin(), array.end());
	it = std::unique (array.begin(), array.end());
	array.resize( std::distance(array.begin(),it));
	return;
}

void printVectorValue(std::vector<double>& worst, std::map<std::vector<double>, double>& valueMap)
{
	for(int j=0;j<worst.size();j++) std::cout<<worst[j]<<" ";
	std::cout<<": "<<valueMap[worst]<<std::endl;
	return;
}

void printSelectedPopulaceFitness(std::vector<std::vector<double>>& array, std::map<std::vector<double>, double>& valueMap)
{
	for(int i=0;i<array.size();i++) printVectorValue(array[i],valueMap);
	return;
}

///////////////////////////////////////////
/*
NIJE DOBRO MORAM OBRNUTI REZULTAT POPULACIJE TE MORAM POMAKNUTI DA SE NADJE NAJBOLJI!!!!!
*/

void geneticAlgorithm(AbstractFunction& Class,std::vector<double>& result,double brPop,double prob, double brEval,double borderLeft, double borderRight, double vectorSize, double numBites=1, double n=3, double mode=1,double typeOfCrossover=1, double typeOfMutation=1, double probability=0.4, double problem=1)
{
	//Problem=1: Minimizacija
	//Problem=2: Maksimizacija
	//std::cout<<brEval<<std::endl;
	std::vector<std::vector<double>> array;
	std::vector<std::vector<double>> selectedarray;
	std::vector<double> child,unit,temp,parent1,parent2,worst;
	std::map<std::vector<double>, double> valueMap;
	int iter=0;
	double best,best2,ran,r,tempValue;
	createPopulace(array,brPop,borderLeft,borderRight,vectorSize, numBites);
	std::cout<<"Population:"<<std::endl;
	printPopulace(array);
	//evaluatePopulace(Class,array,valueMap,problem);
	//translacijaDobrote(valueMap);
	do
	{
		iter++;
		evaluatePopulace(Class,array,valueMap,problem);
		translacijaDobrote(valueMap);
		//valueMap.clear();
		//treba izbrisati i ponovno izracunati valueMap
		//evaluatePopulace(Class,array,valueMap,problem);
		
		//std::cout<<"Iteration: "<<iter<<" Population size: "<<array.size()<<std::endl;		
		//std::cout<<"Population:"<<std::endl;
		//printPopulace(array);
		//printPopulaceFitness(valueMap);
		
		
		selectedarray.clear();
		nTurnirSelecion(selectedarray,array,valueMap,n,mode);

		//////////
		//std::cout<<"Selected chromosomes from population:"<<std::endl;
		//printPopulace(selectedarray);
		//printSelectedPopulaceFitness(selectedarray,valueMap);		
		//////////

		if(selectedarray.size()==2)
		{
			parent1=selectedarray[0];
			parent2=selectedarray[1];
			ran=(double) rand() / (RAND_MAX);
			if(ran<0.5) eliminateWorst(parent1,array,valueMap);
			else eliminateWorst(parent2,array,valueMap);
		}
		else
		{
			selectParents(selectedarray,parent1,parent2,worst,valueMap);
			eliminateWorst(worst,array,valueMap);
		}

		//////////
		//std::cout<<"Population without worst selected:"<<std::endl;
		//printPopulace(array);
		//std::cout<<"Worst: "<<std::endl;
		//printVector(worst);
		//printVectorValue(worst,valueMap); //Value ce uvijek biti 0 jer ga se brise!!!	
		//std::cout<<"Parents: "<<std::endl;
		//printVector(parent1);
		//printVector(parent2);
		//printVectorValue(parent1,valueMap);
		//printVectorValue(parent2,valueMap);
		//////////
		
		crossover(child, parent1, parent2, typeOfCrossover, numBites, borderLeft, borderRight);

		//std::cout<<"Child: "<<std::endl;
		//printVector(child);

		mutation(child, probability, typeOfMutation, numBites, borderLeft, borderRight);		
		//Remove Duplicates
		removeDuplicates(array);
		if(array.size()<3) break;
		/*
		if(valueMap.find(child)!=valueMap.end())
		{
			if (problem==1) valueMap[child]=-Class.function(child);
			else valueMap[child]=Class.function(child)
		}
		*/
		array.push_back(child);	
				
		/*
		TREBA DRUGACIJE NAPRAVITI OVO!!!!!		
		*/
		while(array.size()<brPop)
		{
			removeDuplicates(array);
			//child.clear();
			ran=(double) rand() / (RAND_MAX);
			if(ran<0.5)
			{
				nTurnirSelecion(selectedarray,array,valueMap,n,mode);
				selectParents(selectedarray,parent1,parent2,worst,valueMap);
				crossover(child, parent1, parent2, typeOfCrossover, numBites, borderLeft, borderRight);
				mutation(child, probability, typeOfMutation, numBites, borderLeft, borderRight);
				array.push_back(child);
			}
			else
			{
				temp.clear();
				for(int j=0;j<vectorSize;j++)
				{
					r = ((double) rand() / (RAND_MAX));
					temp.push_back(borderLeft+r*(borderRight-borderLeft));
				}
				array.push_back(temp);
			}
		}
		
		//valueMap.clear();
		//evaluatePopulace(Class,array,valueMap,problem);
		//translacijaDobrote(valueMap);
		//std::cout<<"Population with child:"<<std::endl;
		//printPopulace(array);
		//Treba implementirati elitizam
		if(iter%100==0)
		{
			unit.clear();
			std::cout<<"Iteration: "<<iter<<" ";//<<std::endl;
			findBest(unit,valueMap);
			std::cout<<"Best unit: ";
			printVector(unit);
			std::cout<<"Best value: "<<valueMap[unit]<<std::endl;
			//printPopulaceFitness(valueMap);
		}
		//else
		//{
		//	std::cout<<"Iteration: "<<iter<<" Population size: "<<array.size()<<std::endl;
		//}
		valueMap.clear();
		if(Class.function(child)<1e-6) break;
		if(iter>=brEval) break;
	}
	while(true);
	//petlja moze stati kada funkcija cilja padne ispod 1e-6
	findBest(unit,valueMap);
	result=unit;
	//std::cout<<valueMap[unit]<<std::endl;
	return;
}


int main(int argc, char* argv[])
{
	double brPopulacije,brojBitova,vjerojatnost,brEval,borderLeft, borderRight, velicinaVektora, vrstaMutacija, vrstaKrizanja,velicinaTurnira,problem;
	int zadatak=1;
	if(argc<3)
	{
		std::cout<<"Nisi zadao ulazni txt file ili zadatak koji se treba vrtiti\nPrekidam program."<<std::endl;
		std::exit(0);
	}
	zadatak=atof(argv[2]);
	openFile(argv[1], brPopulacije, brojBitova, vjerojatnost, brEval, borderLeft, borderRight, velicinaVektora, vrstaKrizanja, vrstaMutacija,velicinaTurnira,problem);	

	//int *array;
	//array=(int *) malloc(1*sizeof(int));
	 /* initialize random seed: */
  	srand (time(NULL));
	//free(array);
	
	//std::vector<double> number,x;
	//number.push_back(15);
	//number.push_back(9);
	//x.push_back(1);
	//x.push_back(9);
	
	//std::vector<double> result;
	//getToBinaryString(result ,number , 10);
	//std::cout<<"Binary numbers: ";
	//for(int i=0;i<result.size();i++) std::cout<<result[i]<<" ";
	//std::cout<<std::endl;
	//getFromBinaryString(x,result,10);
	//std::cout<<"Decimal numbers: ";
	//for(int i=0;i<x.size();i++) std::cout<<x[i]<<" ";

	/*
	std::vector<std::vector<double>> array;
	array.push_back(x);
	array.push_back(x);
	array.push_back(number);
	array.push_back(number);
	array.push_back(x);
	array.push_back(number);
	array.push_back(x);
	printPopulace(array);
	removeDuplicates(array);
	std::cout<<"Removed duplicates: "<<std::endl;
	printPopulace(array);
	*/	
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
	//std::cout<<(double)(int)(16*((double) rand() / (RAND_MAX)))<<std::endl;
	//std::cout<<(double)(int)(16*((double) rand() / (RAND_MAX)))<<std::endl;

	//std::vector<double> result;
	//function3 func3;

	//geneticAlgorithm(func3,result,brPopulacije,vjerojatnost,brEval,borderLeft,borderRight,velicinaVektora,brojBitova,3,1,vrstaKrizanja,vrstaMutacija,vjerojatnost);
	if(zadatak==0)
	{
		function3 func6;
		std::vector<double> result6;
		geneticAlgorithm(func6,result6,brPopulacije,vjerojatnost,brEval,borderLeft,borderRight,velicinaVektora,brojBitova,velicinaTurnira,1,vrstaKrizanja,vrstaMutacija,vjerojatnost,problem);		


	}
	if(zadatak==1)
	{
		function3 func3;
		function6 func6;
		function7 func7;

		std::vector<double> result3;
		std::vector<double> result6;
		std::vector<double> result7;

		geneticAlgorithm(func3,result3,brPopulacije,vjerojatnost,brEval,borderLeft,borderRight,5,brojBitova,velicinaTurnira,1,vrstaKrizanja,vrstaMutacija,vjerojatnost,problem);

		geneticAlgorithm(func6,result6,brPopulacije,vjerojatnost,brEval,borderLeft,borderRight,2,brojBitova,velicinaTurnira,1,vrstaKrizanja,vrstaMutacija,vjerojatnost,problem);	
		
		geneticAlgorithm(func7,result7,brPopulacije,vjerojatnost,brEval,borderLeft,borderRight,2,brojBitova,velicinaTurnira,1,vrstaKrizanja,vrstaMutacija,vjerojatnost,problem);
	}

	if(zadatak==2)
	{
		function6 func6;
		function7 func7;
	
		std::vector<double> dim;
	
		dim.push_back(1);
		dim.push_back(3);
		dim.push_back(6);
		dim.push_back(10);

		std::vector<double> result6;
		std::vector<double> result7;

		for(int i=0;i<dim.size();i++)
		{
			geneticAlgorithm(func6,result6,brPopulacije,vjerojatnost,brEval,borderLeft,borderRight,dim[i],brojBitova,velicinaTurnira,1,vrstaKrizanja,vrstaMutacija,vjerojatnost,problem);
			geneticAlgorithm(func7,result7,brPopulacije,vjerojatnost,brEval,borderLeft,borderRight,dim[i],brojBitova,velicinaTurnira,1,vrstaKrizanja,vrstaMutacija,vjerojatnost,problem);
			std::cout<<"Funkcija 6 dimenzija: "<<dim[i]<<" rezultat: ";
			printVector(result6);
			std::cout<<"Funkcija 7 dimenzija: "<<dim[i]<<" rezultat: ";
			printVector(result7);
			result6.clear();
			result7.clear();
		}
	}
	if(zadatak==3)
	{
		function6 func6;
		function7 func7;
		double n = ceil((log(1+(borderRight-borderLeft)*pow(10,4)))/(log(2)));
		std::vector<double> dim;
		double median63,median66,median73,median76;
		double brpog63=0,brpog66=0,brpog73=0,brpog76=0;
	
		dim.push_back(3);
		dim.push_back(6);

		std::vector<double> result6;
		std::vector<double> result7;
		std::vector<double> func6result3;
		std::vector<double> func7result3;		
		std::vector<double> func6result6;
		std::vector<double> func7result6;
		

		for(int i=0;i<dim.size();i++)
		{
			geneticAlgorithm(func6,result6,brPopulacije,vjerojatnost,brEval,borderLeft,borderRight,dim[i],n,velicinaTurnira,1,vrstaKrizanja,vrstaMutacija,vjerojatnost,problem);	
			geneticAlgorithm(func7,result7,brPopulacije,vjerojatnost,brEval,borderLeft,borderRight,dim[i],n,velicinaTurnira,1,vrstaKrizanja,vrstaMutacija,vjerojatnost,problem);
			if(i==0)
			{
				func6result3.push_back(func6.function(result6));
				func7result3.push_back(func7.function(result7));
			}
			else
			{
				func6result6.push_back(func6.function(result6));
				func7result6.push_back(func7.function(result7));
			}
			result6.clear();
			result7.clear();
		}


		median63 = CalcMHWScore(func6result3);
		median66 = CalcMHWScore(func6result6);
		median73 = CalcMHWScore(func7result3);
		median76 = CalcMHWScore(func7result6);
		for(int j=0;j<10;j++)
		{
			if(func6result3[j]<1e-6) brpog63++;
			if(func6result6[j]<1e-6) brpog66++;
			if(func6result3[j]<1e-6) brpog73++;
			if(func7result6[j]<1e-6) brpog76++;
		}
		std::cout<<"Broj pogodaka za funkciju 6, dimenzija 3: "<<brpog63<<" Median: "<<median63<<std::endl;
		std::cout<<"Broj pogodaka za funkciju 6, dimenzija 6: "<<brpog66<<" Median: "<<median66<<std::endl;
		std::cout<<"Broj pogodaka za funkciju 7, dimenzija 3: "<<brpog73<<" Median: "<<median73<<std::endl;
		std::cout<<"Broj pogodaka za funkciju 7, dimenzija 6: "<<brpog76<<" Median: "<<median76<<std::endl;
		printVector(func6result3);
		printVector(func7result3);
	}	
	if(zadatak==4)
	{
		function7 func6;
		std::vector<double> result6;
		std::vector<double> temp,medianPop,medianMul;;
		std::vector<std::vector<double>> paramPop;
		std::vector<std::vector<double>> paramMut;
		//double medianPop,medianMul;
		for(int i=20;i<200;i+=20)
		{
			temp.clear();
			for(int j=0;j<10;j++)
			{				
				geneticAlgorithm(func6,result6,i,vjerojatnost,brEval,borderLeft,borderRight,velicinaVektora,brojBitova,velicinaTurnira,1,vrstaKrizanja,vrstaMutacija,vjerojatnost,problem);
				temp.push_back(func6.function(result6));
				result6.clear();	
			}
			paramPop.push_back(temp);
		}
		std::ofstream myfile;
  		myfile.open ("example.dat");
		for(int j=0;j<9;j++)
		{
			myfile << j << " ";
	  		for(int i=0;i<paramPop.size();i++)
			{
				myfile << paramPop[i][j] <<" ";
			}
			myfile << "\n";
		}
  		myfile.close();
		/*
		for(double i=0.1;i<=1;i+=0.1)
		{	
			temp.clear();
			for(int j=0;j<10;j++)
			{
				geneticAlgorithm(func6,result6,brPopulacije,i,brEval,borderLeft,borderRight,velicinaVektora,brojBitova,velicinaTurnira,1,vrstaKrizanja,vrstaMutacija,vjerojatnost,problem);
				temp.push_back(func6.function(result6));
				result6.clear();	
			}
			paramMut.push_back(temp);
		}
		for(int i=0;i<paramPop.size();i++)
		{
			medianPop.push_back(CalcMHWScore(paramPop[i]));
			medianMul.push_back(CalcMHWScore(paramMut[i]));
		}
		*/
	}
	if(zadatak==5)
	{
		function7 func7;
		std::vector<double> result7;
		std::vector<double> turnirrulet;
		std::vector<double> turnir;
		std::vector<double> rulet;
		double medianturnirruletp,medianrulet,medianturnir;

		for(int i=3;i<10;i+=2)	// parametar 1 određuje da li je turnir ili nije, ovo je turnir i rulet
		{
			std::cout<<i<<std::endl;
			geneticAlgorithm(func7,result7,brPopulacije,vjerojatnost,brEval,borderLeft,borderRight,velicinaVektora,brojBitova,i,1,vrstaKrizanja,vrstaMutacija,vjerojatnost,problem);
			turnirrulet.push_back(func7.function(result7));
			result7.clear();	
		}
		medianturnirruletp=CalcMHWScore(turnirrulet);
		for(int i=3;i<10;i+=2)	// ovo je rulet
		{
			std::cout<<i<<std::endl;
			geneticAlgorithm(func7,result7,brPopulacije,vjerojatnost,brEval,borderLeft,borderRight,velicinaVektora,brojBitova,2,1,vrstaKrizanja,vrstaMutacija,vjerojatnost,problem);
			rulet.push_back(func7.function(result7));
			result7.clear();	
		}
		medianrulet=CalcMHWScore(rulet);
		for(int i=3;i<10;i+=2)	// ovo je turnir
		{
			std::cout<<i<<std::endl;
			geneticAlgorithm(func7,result7,brPopulacije,vjerojatnost,brEval,borderLeft,borderRight,velicinaVektora,brojBitova,i,2,vrstaKrizanja,vrstaMutacija,vjerojatnost,problem);
			turnir.push_back(func7.function(result7));
			result7.clear();	
		}
		medianturnir=CalcMHWScore(turnir);
	}
	return 0;
}
