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

#include "utility.h"

void geneticAlgorithm(AbstractFunction&, std::vector<double>& ,double ,double , double, double, double, double, double );

//Funkcija dobrote
void windowing(std::map<std::vector<double>, double>& ,double , double );
void translacijaDobrote(std::map<std::vector<double>, double>& );

//Create populace
void createPopulace(std::vector<std::vector<double>>& ,double ,double , double , double ,double );

//Evaluate populace
void evaluatePopulace(AbstractFunction& , std::vector<std::vector<double>>& , std::map<std::vector<double>, double>&, double );
//Get chromosome value
double getValue(std::vector<double> , AbstractFunction& ,std::map<std::vector<double>, double>& );

//Binary transformation
void TurnToBinary(std::vector<double>& ,double , double ,double );
void TurnFromBinary(std::vector<double>& ,double , double ,double );
