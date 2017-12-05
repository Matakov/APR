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

class Matrica{
    	private:
		int row,column;
		double** arrayPointer;
		std::vector<std::vector<double>> array;

	public:
		//class constructors
		Matrica();
		Matrica(int l, int k);
		Matrica(const Matrica& B);
		Matrica(const std::vector<std::vector<double>> B);
		Matrica(std::string filename);
		//class destructor
		~Matrica(void);
		//function to display array
		void printMatrix();
		//get number of rows
		int getRow() const;
		//set number od rows
		bool setRow(int r);
		//get number of columns
		int getColumn() const ;
		//set number of columns
		bool setColumn(int c);
		//get element
		double getElement(int i,int j);
		//set element
		bool setElement(int i,int j, double elem);
		//print to file
		bool printToFile(std::string filename);
		//adding operators
		Matrica& operator+=(const Matrica& B);
		//subtraction operators
		Matrica& operator-=(const Matrica& B);
		//copy array for assigning operator
		double** copyArray();
		//assigning operator
		Matrica& operator=(Matrica B);
		//adding operators
		Matrica operator+(const Matrica& B);
		//subtraction operators
		Matrica operator-(const Matrica& B);
		//multiply scalar
		Matrica& operator*=(const double b);
		double multiplyVectors(const Matrica&B, int k, int l);
		//multiplying operators
		Matrica operator*(const Matrica& B);
		//transpose operators
		Matrica transpose();
		//transpose self
		void transposeSelf();
		//comparing (==) operator
		bool operator==(const Matrica& B);
		//get column vector
		Matrica getColumnVector(int column);
		//get row vector
		Matrica getRowVector(int row);
		//LU decomposition
		void LUdekompozicija();
		//supstitucija unaprijed
		void supstitucijaUnaprijed(const Matrica& b, std::vector<double>& result);
		//supstitucija unazad
		void supstitucijaUnazad(const std::vector<double>& y, std::vector<double>& result);
		//napravi jedinicnu matricu
		Matrica JedinicnaMatrica();
		//permutiraj red
		void permutirajRed(Matrica& jedinicna,int pivot);
		//LUP dekompozicija matrice
		Matrica LUPdekompozicija();
		//Multiply matrix with scalar
		Matrica multiplyScalar(int b);
		//Divide matrix with scalar
		Matrica divideScalar(int b);
};
