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

//----------------------------- LABOS 1 MATRICE
extern bool x = false;

class Matrica{
    private:
	int row,column;
	double** arrayPointer;
	std::vector<std::vector<double>> array;

    public:
	//class constructors
	Matrica()
	{
	row=0;
	column=0;
	arrayPointer=nullptr;
	} 	
	Matrica(int l, int k): row(l),column(k)
	{
		/*
		arrayPointer = new double*[row];
		if(row)
		{
			arrayPointer[0] = new double[row*column];
			for (int i=1;1<row;++i)
				arrayPointer[i] = arrayPointer[0] + i*column;
		}
		*/
		arrayPointer = new double*[row];
		for(int i=0;i<row;i++)
		{
			arrayPointer[i]= new double[column];
		}
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<column;j++)
			{
				arrayPointer[i][j]= 0;
			}
		}
	}
	Matrica(const Matrica& B)
	{
		row=B.getRow();
		column=B.getColumn();
		arrayPointer = new double*[row];
		for(int i=0;i<row;i++)
		{
			arrayPointer[i]= new double[column];
		}
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<column;j++)
			{
				arrayPointer[i][j]= B.arrayPointer[i][j];
			}
		}
	}

	Matrica(const std::vector<std::vector<double>> B)
	{
		row=B.size();
		column=B[0].size();
		arrayPointer = new double*[row];
		for(int i=0;i<row;i++)
		{
			arrayPointer[i]= new double[column];
		}
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<column;j++)
			{
				arrayPointer[i][j]=B[i][j];
			}
		}
	}
	
	Matrica(std::string filename)
	{
		std::fstream myfile(filename);
		//myfile.open(filename);
		//int k=0;
		//int l=0;
		std::vector<double> vec;
		std::string line;
		double variable;
		if(myfile.is_open())
			{
			while(!myfile.eof())
				{
				vec.clear();
				std::getline(myfile, line);
				std::istringstream iss(line);
				while (iss >> variable)
				{
					//std::cout << variable << std::endl;
					vec.push_back(variable);
				}
				if(!vec.empty()) array.push_back(vec);
			}			
		}
		//number of rows is determined by the number of rows in a file
		row = array.size();
		//number of columns is determind by the nuber of columns in a row
		column = array[0].size();

		//add to arrayPointer
		arrayPointer = new double*[row];
		for(int i=0;i<row;++i)
		{
			arrayPointer[i]= new double[column];
		}
		for(int i=0;i<row;++i)
		{
			for(int j=0;j<column;++j)
			{
				arrayPointer[i][j]=array[i][j];
			}
		}		
		
	}
	//class destructor
	~Matrica(void)
	{
		if(arrayPointer){
		for(int i=0;i<row;i++)
		{
			delete []arrayPointer[i];
		}
		delete []arrayPointer;
		}
		//std::cout<<"Object has been deleted"<<std::endl;
	}
	//function to display array
	void printMatrix()
	{
		/*
		//std::cout<<array.size()<<std::endl;
		for(auto const& line: array)
		{
			//std::cout<<line.size()<<std::endl;
			for(auto const& elem: line)
			{
				std::cout.width(6);
				std::cout<< elem << " ";
			}
			std::cout<<std::endl;
		}
		*/
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<column;j++)
			{
				std::cout.width(6);
				std::cout<< arrayPointer[i][j] << " ";
			}
			std::cout<<std::endl;
		}
	}
	//get number of rows
	int getRow() const {return row;}
	//set number od rows
	bool setRow(int r)
	{
		//if 
		if(r != row)
		try
		{	
			//need to allocate new space
			double** newArrayPointer = new double*[r];
			for(int i=0;i<r;++i)
			{
				newArrayPointer[i]= new double[column];
			}
			//need to copy elements from old array to new
			for(int i=0;i<r;i++)
			{
				for(int j=0;j<column;j++)
				{
					if(i<row)
					{
						newArrayPointer[i][j]=this->arrayPointer[i][j];
					}
					else
					{
						newArrayPointer[i][j]=0;
					}
					
				}
			}
			//print array
			/*
			for(int i=0;i<std::max(row,r);i++)
			{
				for(int j=0;j<column;j++)
				{
					std::cout.width(5);
					std::cout<<newArrayPointer[i][j];
				}
				std::cout<<std::endl;
			}
			*/
			this->row = r;
			delete []this->arrayPointer;
			this->arrayPointer= newArrayPointer;
			return true;
		}
		catch(int e)
		{
			return false;
		}
	}
	//get number of columns
	int getColumn() const {return column;}
	//set number of columns
	bool setColumn(int c)
	{
		try
		{	
			//need to allocate new space
			double** newArrayPointer = new double*[this->row];
			for(int i=0;i<row;++i)
			{
				newArrayPointer[i]= new double[c];
			}
			//need to copy elements from old array to new
			for(int i=0;i<this->row;i++)
			{
				for(int j=0;j<c;j++)
				{
					if(j<column)
					{
						newArrayPointer[i][j]=this->arrayPointer[i][j];
					}
					else
					{
						newArrayPointer[i][j]=0;
					}
					
				}
			}
			//print array
			/*
			for(int i=0;i<std::max(row,r);i++)
			{
				for(int j=0;j<column;j++)
				{
					std::cout.width(5);
					std::cout<<newArrayPointer[i][j];
				}
				std::cout<<std::endl;
			}
			*/
			this->column = c;
			delete []this->arrayPointer;
			this->arrayPointer= newArrayPointer;
			return true;
		}
		catch(std::exception& e)
		{
			std::cout<<e.what()<<std::endl;
			return false;
		}
	}
	//get element
	double getElement(int i,int j){return arrayPointer[i][j];}
	//set element
	bool setElement(int i,int j, double elem)
	{
		if(i >= row || j >= column)
		{
			std::cout<<"Matrix indexes are not corrent"<<std::endl;
			return false;
		}
		try
		{
			arrayPointer[i][j]=elem;
			return true;
		}
		catch(std::exception& e)
		{
			std::cout<<e.what()<<std::endl;
			return false;		
		}
	}
	//print to file
	bool printToFile(std::string filename)
	{
		std::ofstream myfile(filename);
		if(myfile.is_open())
		{
			for(int i=0;i<row;i++)
			{
				for(int j=0;j<column;j++)
				{
					myfile<<arrayPointer[i][j]<<" ";
				}
				myfile<<std::endl;
			}
			myfile.close();
			return true;
		}
		else
		{
			return false;
		}
		
	}
	//adding operators
	Matrica& operator+=(const Matrica& B)
	{	
		if((this->getRow() != B.getRow()) || (this->getColumn() != B.getColumn()))
		{
			std::cout<<"Matrix dimensions are not the same!"<<std::endl;
			exit(1);
		}
		for(int i=0;i<this->row;i++)
		{
			for(int j=0;j<this->column;j++)
			{
				//std::cout<<this->arrayPointer[i][j]<<" "<<B.arrayPointer[i][j]<<std::endl;
				this->arrayPointer[i][j]= this->arrayPointer[i][j] + B.arrayPointer[i][j];
				
			}		
		}
		return *this;
	}
	
	//subtraction operators
	Matrica& operator-=(const Matrica& B)
	{	
		if((this->getRow() != B.getRow()) || (this->getColumn() != B.getColumn()))
		{
			std::cout<<"Matrix dimensions are not the same!"<<std::endl;
			exit(1);
		}
		for(int i=0;i<this->row;i++)
		{
			for(int j=0;j<this->column;j++)
			{
				//std::cout<<this->arrayPointer[i][j]<<" "<<B.arrayPointer[i][j]<<std::endl;
				this->arrayPointer[i][j]= this->arrayPointer[i][j] - B.arrayPointer[i][j];
				
			}		
		}
		return *this;
	}
	//copy array for assigning operator
	double** copyArray()
	{
		double** newArrayPointer = new double*[this->row];
		for(int i=0;i<this->row;++i)
		{
			newArrayPointer[i]= new double[this->column];
		}
		for(int i=0;i<this->row;i++)
		{
			for(int j=0;j<this->column;j++)
			{
				newArrayPointer[i][j]=this->arrayPointer[i][j];
			}
		}
		return newArrayPointer;
	}
	//assigning operator
	Matrica& operator=(Matrica B)
	{	
		if(arrayPointer) delete[] this->arrayPointer;
		this->row=B.getRow();
		this->column=B.getColumn();
		//this->arrayPointer=B.copyArray();
		//std::swap(*this,B);
		double** arrayPointer = new double*[row];
		for(int i=0;i<row;++i)
		{
			this->arrayPointer[i]= new double[column];
		}
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<column;j++)
			{
				this->arrayPointer[i][j]=B.arrayPointer[i][j];
			}
		}
		return *this;
	}
	//adding operators
	Matrica operator+(const Matrica& B)
	{	
		if((this->getRow() != B.getRow()) || (this->getColumn() != B.getColumn()))
		{
			std::cout<<"Matrix dimensions are not the same!"<<std::endl;
			exit(1);
		}
		Matrica A(this->getRow(),this->getColumn());
		for(int i=0;i<this->row;i++)
		{
			for(int j=0;j<this->column;j++)
			{
				//std::cout<<this->arrayPointer[i][j]<<" "<<B.arrayPointer[i][j]<<std::endl;
				A.arrayPointer[i][j]=this->arrayPointer[i][j]+B.arrayPointer[i][j];
				//std::cout<<A.arrayPointer[i][j]<<std::endl;
			}		
		}
		return A;
	}
	//subtraction operators
	Matrica operator-(const Matrica& B)
	{	
		if((this->getRow() != B.getRow()) || (this->getColumn() != B.getColumn()))
		{
			std::cout<<"Matrix dimensions are not the same!"<<std::endl;
			exit(1);
		}
		Matrica A(this->getRow(),this->getColumn());
		for(int i=0;i<this->row;i++)
		{
			for(int j=0;j<this->column;j++)
			{
				//std::cout<<this->arrayPointer[i][j]<<" "<<B.arrayPointer[i][j]<<std::endl;
				A.arrayPointer[i][j]=this->arrayPointer[i][j]-B.arrayPointer[i][j];
				//std::cout<<A.arrayPointer[i][j]<<std::endl;
			}		
		}
		return A;
	}
	//multiply scalar
	Matrica& operator*=(const double b)
	{	
		for(int i=0;i<this->row;i++)
		{
			for(int j=0;j<this->column;j++)
			{
				//std::cout<<this->arrayPointer[i][j]<<" "<<B.arrayPointer[i][j]<<std::endl;
				this->arrayPointer[i][j]= this->arrayPointer[i][j] * b;
				
			}		
		}
		return *this;
	}
	double multiplyVectors(const Matrica&B, int k, int l)
	{
		double score=0;
		for(int i=0;i<this->getColumn();i++)
		{
			score+=this->arrayPointer[k][i]*B.arrayPointer[i][l];	
		}
		return score;
	}
	//multiplying operators
	Matrica operator*(const Matrica& B)
	{	
		if(this->getColumn() != B.getRow())
		{
			std::cout<<"Matrix dimensions are not the same!"<<std::endl;
			exit(1);
		}
		//std::cout<<"Pass"<<std::endl;
		//std::cout<<this->getRow()<<" "<<B.getColumn()<<std::endl;
		Matrica A(this->getRow(),B.getColumn());
			
		for(int i=0;i<this->getRow();i++)
		{
			for(int j=0;j<B.getColumn();j++)
			{
				double score=0;
				for(int k=0;k<this->getColumn();k++)
				{
					//std::cout<<this->arrayPointer[i][k]<<" "<<B.arrayPointer[k][j]<<std::endl;
					score+=this->arrayPointer[i][k]*B.arrayPointer[k][j];	
				}
				A.arrayPointer[i][j]=score;
			}
		}			
		return A;
		
	}
	//transpose operators
	Matrica transpose()
	{	
		Matrica A(this->column,this->row);
		//std::cout<<this->row<<" "<<this->column<<std::endl;
		//std::cout<<A.getRow()<<" "<<A.getColumn()<<std::endl;
		
		for(int i=0;i<this->row;i++)
		{
			for(int j=0;j<this->column;j++)
			{
				A.arrayPointer[j][i]=this->arrayPointer[i][j];
			}		
		}
		
		return A;
	
	}
	//transpose self
	void transposeSelf()
	{	
		Matrica A(this->column,this->row);
		//std::cout<<this->row<<" "<<this->column<<std::endl;
		//std::cout<<A.getRow()<<" "<<A.getColumn()<<std::endl;
		
		for(int i=0;i<this->row;i++)
		{
			for(int j=0;j<this->column;j++)
			{
				A.arrayPointer[j][i]=this->arrayPointer[i][j];
			}		
		}
		
		for(int i=0;i<this->row;i++)
		{
			for(int j=0;j<this->column;j++)
			{
				this->arrayPointer[i][j]=A.arrayPointer[i][j];
			}		
		}
		
	
	}
	//comparing (==) operator
	bool operator==(const Matrica& B)
	{
		for(int i=0;i<this->row;i++)
		{
			for(int j=0;j<this->column;j++)
			{
				if(this->arrayPointer[i][j]!=B.arrayPointer[i][j])
				{
					return false;
				}
			}
		}
		return true;
	}
	
	//get column vector
	Matrica getColumnVector(int column)
	{
		Matrica A(1,this->row);
		for(int i=0;i<this->row;i++)
		{
			A.arrayPointer[0][i]=this->arrayPointer[i][column];
		}
		return A;
	}
	//get row vector
	Matrica getRowVector(int row)
	{
		Matrica A(1,this->column);
		for(int i=0;i<this->column;i++)
		{
			A.arrayPointer[0][i]=this->arrayPointer[row][i];
		}
		return A;
	}
	//LU decomposition
	void LUdekompozicija()
	{
		for(int i=0;i<this->getRow()-1;i++)
		{
			for(int j=i+1;j<this->getRow();j++)
			{
				if(this->arrayPointer[i][i]==0 || fabs(this->arrayPointer[i][i])<1e-7)
				{
					std::cout<<"Dijelis sa nulom, nema rijesenja."<<std::endl;
					x=true;
					//std::exit(0);
					break;
				}
				this->arrayPointer[j][i]=this->arrayPointer[j][i]/this->arrayPointer[i][i];
				for(int k=i+1;k<this->getRow();k++)
				{
					this->arrayPointer[j][k]-=this->arrayPointer[j][i]*this->arrayPointer[i][k];
				}
			}
		}
	}
	
	//supstitucija unaprijed
	void supstitucijaUnaprijed(const Matrica& b, std::vector<double>& result)
	{
		Matrica sup(this->getRow(),this->getRow());
		for(int i=0;i<this->getRow();i++)
		{
			for(int j=0;j<this->getRow();j++)
			{
				if(i>j)
				{
					sup.arrayPointer[i][j]=this->arrayPointer[i][j];
				}
				else if(i==j)
				{
					sup.arrayPointer[i][j]=1;
				}
				else
				{
					sup.arrayPointer[i][j]=0;
				}
			}
		}
		Matrica y(this->getRow(),1);
		for(int i=0;i<this->getRow();i++)
		{
			y.arrayPointer[i][0]=b.arrayPointer[i][0];
		}
		for(int i=1;i<this->getRow();i++)
		{
			for(int j=0;j<i;j++)
			{
				y.arrayPointer[i][0]-=(sup.arrayPointer[i][j]*y.arrayPointer[j][0]);
//				std::cout<<i<<" "<<j<<" "<<y.arrayPointer[i][0]<<" "<<sup.arrayPointer[i][j]<<" "<<y.arrayPointer[i][0]<<std::endl;
			}
		}
		for(int i=0;i<this->getRow();i++) result[i]=y.arrayPointer[i][0];
		return;
	}
	//supstitucija unazad
	void supstitucijaUnazad(const std::vector<double>& y, std::vector<double>& result)
	{
		Matrica sup(this->getRow(),this->getRow());
		for(int i=0;i<this->getRow();i++)
		{
			for(int j=0;j<this->getRow();j++)
			{
				if(i<=j)
				{
					sup.arrayPointer[i][j]=this->arrayPointer[i][j];
				}
				else
				{
					sup.arrayPointer[i][j]=0;
				}
			}
		}
		Matrica x(this->getRow(),1);
		for(int i=0;i<this->getRow();i++)
		{
			x.arrayPointer[i][0]=y[i];
		}
		for(int i=this->getRow()-1;i>=0;i--)
		{
			for(int j=this->getRow()-1;j>i;j--)
			{
				x.arrayPointer[i][0]-=(sup.arrayPointer[i][j]*x.arrayPointer[j][0]);
			}
			x.arrayPointer[i][0]/=this->arrayPointer[i][i];
			sup.arrayPointer[i][i]/=sup.arrayPointer[i][i];
		}
		for(int i=0;i<y.size();i++) result[i]=x.arrayPointer[i][0];
		return;
	}
	
	//napravi jedinicnu matricu
	Matrica JedinicnaMatrica()
	{
		Matrica jed(this->getRow(),this->getColumn());
		for(int i=0;i<jed.getRow();i++)
		{
			for(int j=0;j<jed.getColumn();j++)
			{
				//std::cout<<i<<" "<<j<<std::endl;
				if(i==j)
				{
					jed.arrayPointer[i][j]=1;
				}
			}
		}
		return jed;
	}
	//permutiraj red
	void permutirajRed(Matrica& jedinicna,int pivot)
	{
		double maxVal=std::numeric_limits<double>::min();
		int pos=0;
		int col=pivot;
		//Matrica jedinicna = this->JedinicnaMatrica();
		//pronadji pivot element
		for(int i=pivot;i<this->getRow();i++)
		{
			if(std::abs(this->arrayPointer[i][col])>maxVal)
			{
				maxVal=std::abs(this->arrayPointer[i][col]);
				pos=i;
			}
		}
		//std::cout<<pos<<" "<<maxVal<<std::endl;
		//permutiraj red
		double* help = this->arrayPointer[pos];
		double* helpRow = new double[this->getColumn()];
		for (int i=0;i<this->getColumn();i++)
		{
			helpRow[i]=*(help+i);
		}
		for (int i=0;i<this->getColumn();i++)
		{
			this->arrayPointer[pos][i]=this->arrayPointer[pivot][i];
		}
		for (int i=0;i<this->getColumn();i++)
		{
			this->arrayPointer[pivot][i]=helpRow[i];
		}
		//std::cout<<*(help+2)<<std::endl;
		//std::cout<<help<<std::endl;
		/*
		this->arrayPointer[pos]=this->arrayPointer[pivot];
		this->arrayPointer[pivot]=help;
		help = jedinicna.arrayPointer[pos];
		jedinicna.arrayPointer[pos]=jedinicna.arrayPointer[pivot];
		jedinicna.arrayPointer[pivot]=help;
		#*/
		help = jedinicna.arrayPointer[pos];
		for (int i=0;i<this->getColumn();i++)
		{
			helpRow[i]=*(help+i);
		}
		for (int i=0;i<this->getColumn();i++)
		{
			jedinicna.arrayPointer[pos][i]=jedinicna.arrayPointer[pivot][i];
		}
		for (int i=0;i<this->getColumn();i++)
		{
			jedinicna.arrayPointer[pivot][i]=helpRow[i];
		}
		
		//this->printMatrix();
		//return jedinicna;
	}
	//LUP dekompozicija matrice
	Matrica LUPdekompozicija()
	{
		Matrica jedinicna = this->JedinicnaMatrica();
		for(int i=0;i<this->getRow();i++)
		{
			this->permutirajRed(jedinicna,i);
			//this->printMatrix();
			//jedinicna.printMatrix();
			for(int j=i+1;j<this->getColumn();j++)
			{
				this->arrayPointer[j][i]=this->arrayPointer[j][i]/this->arrayPointer[i][i];
				for(int k=i+1;k<this->getRow();k++)
				{
					this->arrayPointer[j][k]-=this->arrayPointer[j][i]*this->arrayPointer[i][k];
				}
			}
		}
		return jedinicna;
	}
	//Multiply matrix with scalar
	Matrica multiplyScalar(int b)
	{
		Matrica A(this->getRow(),this->getColumn());
		for(int i=0;i<this->getRow();i++)
		{
			for(int j=0;j<this->getColumn();j++)
			{
				A.arrayPointer[i][j]=this->arrayPointer[i][j]*b;
			}
		}
		return A;
	}
	//Divide matrix with scalar
	Matrica divideScalar(int b)
	{
		Matrica A(this->getRow(),this->getColumn());
		for(int i=0;i<this->getRow();i++)
		{
			for(int j=0;j<this->getColumn();j++)
			{
				A.arrayPointer[i][j]=this->arrayPointer[i][j]/b;
			}
		}
		return A;
	}
};

// KRAJ KLASE



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

//Rosenbrockova 'banana' funkcija
class function1: public AbstractFunction 
{
	public:
		function1():AbstractFunction(){}
		double function(double a, double b)
		{
			increase();
			return 100*pow(b-pow(a,2),2)+pow((1-a),2);
		}
		double function(std::vector<double> lista)
		{	
			if (lista.size()>2)
			{
				std::cout<<"Input vectors are not good!"<<std::endl;
				throw(std::invalid_argument( "Input vectors are not good!" ));
							
			}
			increase();
			return 100*pow(lista[1]-pow(lista[0],2),2)+pow((1-lista[0]),2);
		}
		void restartCount()
		{
			restartCounting();
		}
};

class function2: public AbstractFunction 
{
	public:
		function2():AbstractFunction(){}
		double function(double a, double b)
		{
			increase();
			return pow(a-4,2)+ 4*pow(b-2,2);
		}
		double function(std::vector<double> lista)
		{	
			if (lista.size()>2)
			{
				std::cout<<"Input vectors are not good!"<<std::endl;
				throw(std::invalid_argument( "Input vectors are not good!" ));
							
			}
			increase();
			return pow(lista[0]-4,2)+ 4*pow(lista[1]-2,2);
		}
		void restartCount()
		{
			restartCounting();
		}
};

class function3: public AbstractFunction 
{
	public:
		function3():AbstractFunction(){}
		double function(std::vector<double> lista)
		{	
			if (lista.size()>2)
			{
				std::cout<<"Input vectors are not good!"<<std::endl;
				throw(std::invalid_argument( "Input vectors are not good!" ));
							
			}
			increase();
			//return (lista[0]-2)*(lista[0]+lista[1]))+sqrt(pow(lista[0],2)+pow(lista[1],2));
			return pow(lista[0]-2,2)+pow(lista[1]+3,2);
		}
		double function(double a, double b)
		{	
			increase();
			return abs((a-b)*(a+b))+sqrt(pow(a,2)+pow(b,2));
		}
		void restartCount()
		{
			restartCounting();
		}
};


class function4: public AbstractFunction 
{
	public:
		function4():AbstractFunction(){}
		double function(std::vector<double> lista)
		{
			increase();
			int numDim = lista.size();
			double sum = 0;
			for(int i=0;i<numDim;i++)
			{
				sum +=pow((lista[i]),2);
			}
			return pow(lista[0]-3,2)+pow(lista[1],2);
		}
		void restartCount()
		{
			restartCounting();
		}
};

//Here are implitic functions
class function5: public AbstractFunction 
{
	public:
		function5():AbstractFunction(){}
		double function(std::vector<double> lista)
		{
			increase();
			return lista[1]-lista[0];
		}
		void restartCount()
		{
			restartCounting();
		}
};

class function6: public AbstractFunction 
{
	public:
		function6():AbstractFunction(){}
		double function(std::vector<double> lista)
		{
			increase();
			return 2-lista[0];
		}
		void restartCount()
		{
			restartCounting();
		}
};

void openFile(std::string name,std::vector<double>& tocka,std::vector<double>& minimumFunkcije,std::vector<double>& preciznost, std::vector<double>& pomaciFunkcije,double& leftPoint,double& rightPoint,double& distance)
{
	std::ifstream myfile;
	std::string line;
	myfile.open(name);
	if(myfile.is_open())
	{
		while(getline (myfile,line))
		{
			std::vector<std::string> container = split(line,' ');
			if (container[0] == "Preciznost:")
			{
				for(int i=1;i<container.size();i++) preciznost.push_back(atof(container[i].c_str()));
				//preciznost=atof(container[1].c_str());			
			}
			if (container[0] == "Početna" && container[1] == "točka:")
			{
				for(int i=2;i<container.size();i++) tocka.push_back(atof(container[i].c_str()));			
			}
			if (container[0] == "Minimum" && container[1] == "funkcije:")
			{
				for(int i=2;i<container.size();i++) minimumFunkcije.push_back(atof(container[i].c_str()));			
			}
			if (container[0] == "Pomaci" && container[1] == "funkcije:")
			{
				for(int i=2;i<container.size();i++) pomaciFunkcije.push_back(atof(container[i].c_str()));			
			}
			if (container[0] == "Početni" && container[1] == "interval:")
			{
				leftPoint=atof(container[2].c_str());
				rightPoint=atof(container[3].c_str());			
			}
			if (container[0] == "Razmak" && container[1] == "točaka:")
			{
				distance=atof(container[2].c_str());		
			}
			//std::cout<<container[2]<<std::endl;
			//std::for_each (container.begin(), container.end(), myfunction);
			
		}
	}
}
int main(int argc, char* argv[]){
	std::ifstream myfile;
	std::string line;
	std::vector<double> tocka,minimumFunkcije,preciznost,pomaciFunkcije;
	double leftPoint,rightPoint,distance;
	int zadatak;
	//std::cout<<argc<<std::endl;
	/*	
	if(argc<3)
	{
		std::cout<<"Nisi zadao ulazni txt file!!\nPrekidam program."<<std::endl;
		std::exit(0);
	}
	zadatak=atof(argv[2]);
	openFile(argv[1],tocka,minimumFunkcije,preciznost,pomaciFunkcije,leftPoint,rightPoint,distance);
	*/
	//for(int i=0;i<tocka.size();i++) std::cout<<preciznost[i]<<" ";
	//for(int i=0;i<tocka.size();i++) std::cout<<tocka[i]<<" ";
	//for(int i=0;i<minimumFunkcije.size();i++) std::cout<<minimumFunkcije[i]<<" ";
	//std::cout<<leftPoint<<" "<<rightPoint<<std::endl;
	//std::cout<<func3(tocka)<<std::endl;
	
	/*
	function1 func1;
	double ide = func1.function(1,2);
	std::cout<<ide<<" "<<func1.getNumbers()<<std::endl;
	ide = func1.function(1,3);
	std::cout<<ide<<" "<<func1.getNumbers()<<std::endl;
	func1.restartCount();
	std::cout<<func1.getNumbers()<<std::endl;

	function2 func2;
	ide = func2.function(1,2);
	std::cout<<ide<<" "<<func2.getNumbers()<<std::endl;
	ide = func2.function(1,3);
	std::cout<<ide<<" "<<func2.getNumbers()<<std::endl;
	func2.restartCount();
	std::cout<<func2.getNumbers()<<std::endl;

	std::vector<double>v;
	v.push_back(2);
	v.push_back(2);

	function3 func3;
	ide = func3.function(v);
	std::cout<<ide<<" "<<func3.getNumbers()<<std::endl;
	//v.push_back(3);
	ide = func3.function(v);
	std::cout<<ide<<" "<<func3.getNumbers()<<std::endl;
	func3.restartCount();
	std::cout<<func3.getNumbers()<<std::endl;

	function4 func4;
	ide = func4.function(v);
	std::cout<<ide<<" "<<func4.getNumbers()<<std::endl;
	ide = func4.function(v);
	std::cout<<ide<<" "<<func4.getNumbers()<<std::endl;
	std::cout<<func4.getNumbers()<<std::endl;
	ide = func4.function(2,2);
	std::cout<<ide<<" "<<func4.getNumbers()<<std::endl;
	ide = func4.function(2,2);
	std::cout<<ide<<" "<<func4.getNumbers()<<std::endl;
	func4.restartCount();
	std::cout<<func4.getNumbers()<<std::endl;

	function6 func6;
	ide = func6.function(v);
	std::cout<<ide<<" "<<func6.getNumbers()<<std::endl;
	v.push_back(3);
	ide = func6.function(v);
	std::cout<<ide<<" "<<func6.getNumbers()<<std::endl;
	func6.restartCount();
	std::cout<<func6.getNumbers()<<std::endl;
	*/
	
	double alfa=1;
	double beta=0.5;
	double gamma=2;	
	
	function1 func1;
	function2 func2;
	std::vector<AbstractFunction*> vector;
	vector.push_back(&func1);
	vector.push_back(&func2);
	
	for(int i=0;i<vector.size();i++)
	{
		std::cout<<"Result"<<vector[i]->function(1,2)<<std::endl;
	}
	return 0;

	
}
