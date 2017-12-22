#include "lab1.h"

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
	//std::cout<<"Zadatak 1: "<<std::endl;	
	//Matrica A("1.txt");
	//Matrica B(A);
	//std::cout<<"Matrica A: "<<std::endl;
	//A.printMatrix();
	///std::cout<<"Matrica B: "<<std::endl;
	//B.printMatrix();
	//std::cout<<"Jednakost matrica A i B: "<<std::endl;
	//std::cout<<(A==B)<<std::endl;
	//B.divideScalar(15015.7777777);
	//B.multiplyScalar(15015.7777777);
	//std::cout<<"Jednakost matrica A i B: "<<std::endl;
	//std::cout<<(A==B)<<std::endl;
	
	//ZADATAK 1
	std::cout<<"Zadatak 1: "<<std::endl;
	Matrica _1("zad1.txt");
	//LU
	_1.printMatrix();
	std::vector<std::vector<double>> array;	
	_1.getInverse(array);
	//std::cout<<"Inverse:"<<std::endl;
	//inv.printMatrix();
	for(int i=0;i<array.size();i++)
		{
			for(int j=0;j<array[0].size();j++)
			{
				std::cout<<array[i][j]<<" ";
			}
			std::cout<<std::endl;
		}
	//ZADATAK 2
	std::cout<<std::endl;
	std::cout<<"Zadatak 2: "<<std::endl;
	Matrica _2("zad2.txt");
	//LU
	_2.printMatrix();
	array.clear();
	_2.getInverse(array);
	//std::cout<<"Inverse:"<<std::endl;
	//inv.printMatrix();
	std::cout<<"Inverz: "<<std::endl;
	for(int i=0;i<array.size();i++)
		{
			for(int j=0;j<array[0].size();j++)
			{
				std::cout<<std::setw(6);
				std::cout<<array[i][j]<<" ";
			}
			std::cout<<std::endl;
		}

	return 0;
}
