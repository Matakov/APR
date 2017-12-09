#include<iostream>
#include<iomanip>
#include<sstream>
#include<functional>
#include<fstream>
#include<algorithm>
#include<vector>
#include<stdlib.h>
#include<math.h>
#include<map>
#include<time.h>
#include<algorithm>
#include<random>

#include "utility.h"

void geneticAlgorithm(AbstractFunction& ,std::vector<double>& ,double ,double , double ,double , double , double , double , double , double ,double , double , double );

//Funkcija dobrote
void windowing(std::map<std::vector<double>, double>& ,double , double );
void translacijaDobrote(std::map<std::vector<double>, double>& );

//Create populace
void createPopulace(std::vector<std::vector<double>>& ,double ,double , double , double ,double );

//Evaluate populace
void evaluatePopulace(AbstractFunction& , std::vector<std::vector<double>>& , std::map<std::vector<double>, double>& );
//Get chromosome value
double getValue(std::vector<double> , AbstractFunction& ,std::map<std::vector<double>, double>& );

//Binary transformation
void TurnToBinary(std::vector<double>& ,double , double ,double );
void TurnFromBinary(std::vector<double>& ,double , double ,double );

void getToBinaryString(std::vector<double>& ,std::vector<double> ,double );
void getFromBinaryString(std::vector<double>& ,std::vector<double> ,double );

void binarizeArray(std::vector<std::vector<double>>&,double , double ,double );
void debinarizeArray(std::vector<std::vector<double>>&,double , double ,double );

void binarizeVector(std::vector<double>& ,double , double ,double );
void debinarizeVector(std::vector<double>& ,double , double ,double );

//check if already in array
bool checkIfInArray(std::vector<double>, double );

//parent selection
void nTurnirSelecion(std::vector<std::vector<double>>&, std::vector<std::vector<double>> const &, std::map<std::vector<double>, double>&, double, double );

//removing the worst unit in a populace
void eliminateWorst(std::vector<std::vector<double>>& ,std::vector<std::vector<double>>& ,  std::map<std::vector<double>, double>& );

//generating stohastic unit
std::vector<double> generateR(double );

//crossover operation
void crossover(std::vector<double>& ,std::vector<std::vector<double>>& ,std::map<std::vector<double>, double>& , double , double , double , double );

//mutation operation
void mutation(std::vector<double>& , double , double , double );
