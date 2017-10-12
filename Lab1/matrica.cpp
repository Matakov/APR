#include<iostream>
#include<iomanip>
#include<string>
#include<fstream>
#include<vector>
#include<sstream>
#include<exception>
#include<cmath>
#include <limits>

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
				this->arrayPointer[j][i]=this->arrayPointer[j][i]/this->arrayPointer[i][i];
				for(int k=i+1;k<this->getRow();k++)
				{
					this->arrayPointer[j][k]-=this->arrayPointer[j][i]*this->arrayPointer[i][k];
				}
			}
		}
	}
	
	//supstitucija unaprijed
	Matrica supstitucijaUnaprijed(const Matrica& b)
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
		return y;
	}
	//supstitucija unazad
	Matrica supstitucijaUnazad(const Matrica& y)
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
			x.arrayPointer[i][0]=y.arrayPointer[i][0];
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
		return x;
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





int main()
{
	/*
	Matrica A("B.txt");
	Matrica B("A.txt");
	Matrica C("C.txt");
	//Matrica B("A.txt");
	//Matrica C=A;
	//A.printMatrix();
	//A*=2;
	//std::cout<<A.getRow()<<" "<<A.getColumn()<<std::endl;
	//A.printMatrix();
	//std::cout<<B.getRow()<<" "<<B.getColumn()<<std::endl;
	//B.printMatrix();
	//C.printMatrix();
	//Matrica C=A*B;
	//C.printMatrix();
	//B.LUdekompozicija();
	//std::cout<<"Matrica:"<<std::endl;
	//B.printMatrix();
	//std::cout<<"Vektor b:"<<std::endl;
	//C.printMatrix();
	//Matrica y=B.supstitucijaUnaprijed(C);
	//y.printMatrix();
	//Matrica x=B.supstitucijaUnazad(y);
	//x.printMatrix();
	//A.printMatrix();
	Matrica D=B.JedinicnaMatrica();
	//A.permutirajRed(D,0);
	//D.printMatrix();
	Matrica P=B.LUPdekompozicija();
	B.printMatrix();
	P.printMatrix();
	Matrica y=B.supstitucijaUnaprijed(P*C);
	y.printMatrix();
	Matrica x=B.supstitucijaUnazad(y);
	x.printMatrix();
	*/
	//ZADATAK 1
	std::cout<<"Zadatak 1: "<<std::endl;	
	Matrica A("1.txt");
	Matrica B(A);
	std::cout<<"Matrica A: "<<std::endl;
	A.printMatrix();
	std::cout<<"Matrica B: "<<std::endl;
	B.printMatrix();
	std::cout<<"Jednakost matrica A i B: "<<std::endl;
	std::cout<<(A==B)<<std::endl;
	B.divideScalar(15015);
	B.multiplyScalar(15015);
	std::cout<<"Jednakost matrica A i B: "<<std::endl;
	std::cout<<(A==B)<<std::endl;
	
	
	//ZADATAK 2
	std::cout<<std::endl;
	std::cout<<"Zadatak 2: "<<std::endl;
	Matrica _2a("2a.txt");
	Matrica _2b("2b.txt");
	//LU	
	_2a.LUdekompozicija();
	std::cout<<"Matrica A nakon LU dekompozicije: "<<std::endl;
	_2a.printMatrix();
	Matrica _2ya=_2a.supstitucijaUnaprijed(_2b);
	std::cout<<"Y vektor nakon supstitucije unaprijed: "<<std::endl;
	_2ya.printMatrix();
	Matrica _2xa=_2a.supstitucijaUnazad(_2ya);
	std::cout<<"X vektor nakon supstitucije unazad: "<<std::endl;
	_2xa.printMatrix();
	//LUP
	Matrica _2a2("2a.txt");
	Matrica _2P=_2a2.LUPdekompozicija();
	std::cout<<"Matrica A nakon LUP dekompozicije: "<<std::endl;
	_2a2.printMatrix();
	std::cout<<"Jedinična matrica nakon LUP dekompozicije: "<<std::endl;	
	_2P.printMatrix();
	Matrica _2yb=_2a2.supstitucijaUnaprijed(_2P*_2b);
	std::cout<<"Y vektor nakon supstitucije unaprijed: "<<std::endl;
	_2yb.printMatrix();
	Matrica _2xb=_2a2.supstitucijaUnazad(_2yb);
	std::cout<<"X vektor nakon supstitucije unazad: "<<std::endl;
	_2xb.printMatrix();
	
	//ZADATAK 3
	std::cout<<std::endl;
	std::cout<<"Zadatak 3: "<<std::endl;
	Matrica _3a("3a.txt");
	Matrica _3b("slob.txt");
	//LU	
	_3a.LUdekompozicija();
	std::cout<<"Matrica A nakon LU dekompozicije: "<<std::endl;
	_3a.printMatrix();
	Matrica _3ya=_3a.supstitucijaUnaprijed(_3b);
	std::cout<<"Y vektor nakon supstitucije unaprijed: "<<std::endl;
	_3ya.printMatrix();
	Matrica _3xa=_3a.supstitucijaUnazad(_3ya);
	std::cout<<"X vektor nakon supstitucije unazad: "<<std::endl;
	_3xa.printMatrix();
	//LUP
	Matrica _3a2("3a.txt");
	Matrica _3P=_3a2.LUPdekompozicija();
	std::cout<<"Matrica A nakon LUP dekompozicije: "<<std::endl;
	_3a2.printMatrix();
	std::cout<<"Jedinična matrica nakon LUP dekompozicije: "<<std::endl;	
	_3P.printMatrix();
	Matrica _3yb=_3a2.supstitucijaUnaprijed(_3P*_3b);
	std::cout<<"Y vektor nakon supstitucije unaprijed: "<<std::endl;
	_3yb.printMatrix();
	Matrica _3xb=_3a2.supstitucijaUnazad(_3yb);
	std::cout<<"X vektor nakon supstitucije unazad: "<<std::endl;
	_3xb.printMatrix();
	
	//ZADATAK 4
	std::cout<<std::endl;
	std::cout<<"Zadatak 4: "<<std::endl;
	Matrica _4a("4a.txt");
	Matrica _4b("4b.txt");
	//LU	
	_4a.LUdekompozicija();
	std::cout<<"Matrica A nakon LU dekompozicije: "<<std::endl;
	_4a.printMatrix();
	Matrica _4ya=_4a.supstitucijaUnaprijed(_4b);
	std::cout<<"Y vektor nakon supstitucije unaprijed: "<<std::endl;
	_4ya.printMatrix();
	Matrica _4xa=_4a.supstitucijaUnazad(_4ya);
	std::cout<<"X vektor nakon supstitucije unazad: "<<std::endl;
	_4xa.printMatrix();
	//LUP
	Matrica _4a2("4a.txt");
	Matrica _4P=_4a2.LUPdekompozicija();
	std::cout<<"Matrica A nakon LUP dekompozicije: "<<std::endl;
	_4a2.printMatrix();
	std::cout<<"Jedinična matrica nakon LUP dekompozicije: "<<std::endl;	
	_4P.printMatrix();
	Matrica _4yb=_4a2.supstitucijaUnaprijed(_4P*_4b);
	std::cout<<"Y vektor nakon supstitucije unaprijed: "<<std::endl;
	_4yb.printMatrix();
	Matrica _4xb=_4a2.supstitucijaUnazad(_4yb);
	std::cout<<"X vektor nakon supstitucije unazad: "<<std::endl;
	_4xb.printMatrix();

	//ZADATAK 5
	std::cout<<std::endl;
	std::cout<<"Zadatak 5: "<<std::endl;
	Matrica _5a("5a.txt");
	Matrica _5b("5b.txt");
	//LU	
	_5a.LUdekompozicija();
	std::cout<<"Matrica A nakon LU dekompozicije: "<<std::endl;
	_5a.printMatrix();
	Matrica _5ya=_5a.supstitucijaUnaprijed(_5b);
	std::cout<<"Y vektor nakon supstitucije unaprijed: "<<std::endl;
	_5ya.printMatrix();
	Matrica _5xa=_5a.supstitucijaUnazad(_5ya);
	std::cout<<"X vektor nakon supstitucije unazad: "<<std::endl;
	_5xa.printMatrix();
	//LUP
	Matrica _5a2("5a.txt");
	Matrica _5P=_5a2.LUPdekompozicija();
	std::cout<<"Matrica A nakon LUP dekompozicije: "<<std::endl;
	_5a2.printMatrix();
	std::cout<<"Jedinična matrica nakon LUP dekompozicije: "<<std::endl;	
	_5P.printMatrix();
	Matrica _5yb=_5a2.supstitucijaUnaprijed(_5P*_5b);
	std::cout<<"Y vektor nakon supstitucije unaprijed: "<<std::endl;
	_5yb.printMatrix();
	Matrica _5xb=_5a2.supstitucijaUnazad(_5yb);
	std::cout<<"X vektor nakon supstitucije unazad: "<<std::endl;
	_5xb.printMatrix();

	//ZADATAK 6
	std::cout<<std::endl;
	std::cout<<"Zadatak 6: "<<std::endl;
	Matrica _6a("6a.txt");
	Matrica _6b("6b.txt");
	//LU	
	_6a.LUdekompozicija();
	std::cout<<"Matrica A nakon LU dekompozicije: "<<std::endl;
	_6a.printMatrix();
	Matrica _6ya=_6a.supstitucijaUnaprijed(_6b);
	std::cout<<"Y vektor nakon supstitucije unaprijed: "<<std::endl;
	_6ya.printMatrix();
	Matrica _6xa=_6a.supstitucijaUnazad(_6ya);
	std::cout<<"X vektor nakon supstitucije unazad: "<<std::endl;
	_6xa.printMatrix();
	//LUP	
	Matrica _6a2("6a.txt");
	Matrica _6P=_6a2.LUPdekompozicija();
	std::cout<<"Matrica A nakon LUP dekompozicije: "<<std::endl;
	_6a2.printMatrix();
	std::cout<<"Jedinična matrica nakon LUP dekompozicije: "<<std::endl;	
	_6P.printMatrix();
	Matrica _6yb=_6a2.supstitucijaUnaprijed(_6P*_6b);
	std::cout<<"Y vektor nakon supstitucije unaprijed: "<<std::endl;
	_6yb.printMatrix();
	Matrica _6xb=_6a2.supstitucijaUnazad(_6yb);
	std::cout<<"X vektor nakon supstitucije unazad: "<<std::endl;
	_6xb.printMatrix();

	return 0;	
}

