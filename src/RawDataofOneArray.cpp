#include <cmath>

#include "RawDataofOneArray.h"

#define MAX_ITERATION 30
#define	CONVERGENCE_ABS 0.001
#define GOLDEN_SECTION 0.618

// private member function

double RawDataofOneArray::P22Function(double p22,struct my_f_params_P22 *paramsP22)
{
	double pp11 = (paramsP22->p11);
	double pp12 = (paramsP22->p12);
	double pp21 = (paramsP22->p21);

	vector <double> pinten1 = (paramsP22->inten1);
	vector <double> pinten2 = (paramsP22->inten2);
	
	vector <double>::iterator ppinten1 = pinten1.begin();
	vector <double>::iterator ppinten2 = pinten2.begin();

	double costofone = 0;
	double costofall = 0;

	while (ppinten1!=pinten1.end())  // sum of the cost function.
	{
		costofone = log(*ppinten2) - log(pp12 * (*ppinten1 - pp21)/pp11 + p22);
		costofall = costofall + costofone * costofone;

		ppinten1++;
		ppinten2++;
	}

	return costofall;
}

double RawDataofOneArray::P21Function(double p21,struct my_f_params_P21 *paramsP21)
{
	double pp11 = (paramsP21->p11);
	double pp12 = (paramsP21->p12);
	double pp22 = (paramsP21->p22);

	vector <double> pinten1 = (paramsP21->inten1);
	vector <double> pinten2 = (paramsP21->inten2);
	
	vector <double>::iterator ppinten1 = pinten1.begin();
	vector <double>::iterator ppinten2 = pinten2.begin();

	double costofone = 0;
	double costofall = 0;

	while (ppinten1!=pinten1.end())  // sum of the cost function.
	{
		costofone = log(*ppinten1) - log(pp11 * (*ppinten2 - pp22)/pp12 + p21);
		costofall = costofall + costofone * costofone;

		ppinten1++;
		ppinten2++;
	}

	return costofall;
}

// public member funtion of class RawDataofOneArray

void RawDataofOneArray::setArrayID(int id) 
{
	ArrayID = id;
}

void RawDataofOneArray::setCloneID(vector <string> id)
{
	Data.CloneID = id;
}

void RawDataofOneArray::setRawData(vector<double> inten1, vector<double> inten2) 
{
	Data.Col1Inten = inten1;
	Data.Col2Inten = inten2;
}

int RawDataofOneArray::getArrayID() 
{
	return ArrayID;
}

vector <double> RawDataofOneArray::getCol1Inten()
{
	return Data.Col1Inten; 
}

vector <double> RawDataofOneArray::getCol2Inten()
{
	return Data.Col2Inten; 
}

vector <string> RawDataofOneArray::getCloneID()
{
	return Data.CloneID; 
}

vector <RawDataofOneArray> RawDataofOneArray::selectBlockData(vector <RawDataofOneArray> rawdata,Design des)
{
	vector <int> a;
	a = des.getArray();

	vector <int> a_unique(a.begin(),a.end());
	vector <int>::iterator pa_unique = unique(a_unique.begin(),a_unique.end());
	a_unique.erase(pa_unique,a_unique.end());

	vector <RawDataofOneArray>::iterator praw = rawdata.begin();

	int arrayid;
	RawDataofOneArray rawdataofone;
	vector <RawDataofOneArray> rawdataofoneblock;

	bool found = false;

	while (praw != rawdata.end())
	{
		rawdataofone = *praw;
		arrayid = rawdataofone.getArrayID();

		vector <int> tmp;
		tmp.push_back(arrayid);

		found = includes(a_unique.begin(),a_unique.end(),tmp.begin(),tmp.end());

		if (found)
		{
			rawdataofoneblock.push_back(rawdataofone); 

		}

		praw++;
	}

	return (rawdataofoneblock);
}

double RawDataofOneArray::adjustCy5(ParameterofOneArray parameter) 
{
	vector <double> inten1; 
	vector <double> inten2;

	vector <double>::iterator pinten1 = Data.Col1Inten.begin();

	double p11 = parameter.getP1Col1();
	double p12 = parameter.getP1Col2(); 
	double p21 = parameter.getP2Col1(); 

	int i = 0;

	// pick out points whose yred > p21.

	while (pinten1!=Data.Col1Inten.end())
	{
		if (*pinten1 > p21)
		{
			inten1.push_back(*pinten1);
			inten2.push_back(Data.Col2Inten[i]);
		}

		i++;
		pinten1++;
	}

	// estimation new p22.

	struct my_f_params_P22 paramsP22 = {inten1,inten2,p11,p12,p21};

	double lowerb = 0;
	double upperb = 10 * parameter.getP2Col2();

	double iter = 0;
	double a = lowerb;
	double b = upperb;
	double x1 = a + (1 - GOLDEN_SECTION) * (b - a);
	double x2 = a + GOLDEN_SECTION * (b - a);
	double f1 = P22Function(x1,&paramsP22);
	double f2 = P22Function(x2,&paramsP22);

	while ((abs(b-a) > CONVERGENCE_ABS) && (iter < MAX_ITERATION))
	{
		if (f1 < f2)
		{
			b = x2;
			x2 = x1;
			f2 = f1;
			x1 = a + (1 - GOLDEN_SECTION) * (b - a);
			f1 = P22Function(x1,&paramsP22);

		}else{
			if (f1 == f2)
			{
				a = x1;
				b = x2;
				x1 = a + (1 - GOLDEN_SECTION) * (b - a);
				x2 = a + GOLDEN_SECTION * (b - a);			
				f1 = P22Function(x1,&paramsP22);
				f2 = P22Function(x2,&paramsP22);

			}else{
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + GOLDEN_SECTION * (b - a);
				f2 = P22Function(x2,&paramsP22);
			}
		}

		iter++;

	}

	double min = (a+b)/2;
		
	return 	min;		

}

double RawDataofOneArray::adjustCy3(ParameterofOneArray parameter) 
{
	vector <double> inten1; 
	vector <double> inten2;

	vector <double>::iterator pinten2 = Data.Col2Inten.begin();

	double p11 = parameter.getP1Col1();
	double p12 = parameter.getP1Col2(); 
	double p22 = parameter.getP2Col2(); 

	int i = 0;

	// pick out points whose yred > p21.

	while (pinten2!=Data.Col2Inten.end())
	{
		if (*pinten2 > p22)
		{
			inten2.push_back(*pinten2);
			inten1.push_back(Data.Col1Inten[i]);
		}

		i++;
		pinten2++;
	}

	// estimation new p21.

	struct my_f_params_P21 paramsP21 = {inten1,inten2,p11,p12,p22};

	double lowerb = 0;
	double upperb = 10 * parameter.getP2Col1();

	double iter = 0;
	double a = lowerb;
	double b = upperb;
	double x1 = a + (1 - GOLDEN_SECTION) * (b - a);
	double x2 = a + GOLDEN_SECTION * (b - a);
	double f1 = P21Function(x1,&paramsP21);
	double f2 = P21Function(x2,&paramsP21);

	while ((abs(b-a) > CONVERGENCE_ABS) && (iter < MAX_ITERATION))
	{
		if (f1 < f2)
		{
			b = x2;
			x2 = x1;
			f2 = f1;
			x1 = a + (1 - GOLDEN_SECTION) * (b - a);
			f1 = P21Function(x1,&paramsP21);

		}else{
			if (f1 == f2)
			{
				a = x1;
				b = x2;
				x1 = a + (1 - GOLDEN_SECTION) * (b - a);
				x2 = a + GOLDEN_SECTION * (b - a);			
				f1 = P21Function(x1,&paramsP21);
				f2 = P21Function(x2,&paramsP21);

			}else{
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + GOLDEN_SECTION * (b - a);
				f2 = P21Function(x2,&paramsP21);
			}
		}

		iter++;

	}

	double min = (a+b)/2;
		
	return 	min;		

}
