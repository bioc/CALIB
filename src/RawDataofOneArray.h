#include <vector>
#include <string>

#include "Design.h"
#include "ParameterofOneArray.h"

using namespace std;

struct RawArrayData
{
	vector <double> Col1Inten;
	vector <double> Col2Inten;
	vector <string> CloneID;
};

class RawDataofOneArray
{
	int ArrayID;
	RawArrayData Data;

private:
	double P22Function(double p22,struct my_f_params_P22 *paramsP22);
	double P21Function(double p21,struct my_f_params_P21 *paramsP21);

public:
	void setArrayID(int id);
	void setCloneID(vector <string> id);
	void setRawData(vector <double> inten1,vector <double> inten2);
	int getArrayID();
	vector <double> getCol1Inten();
	vector <double> getCol2Inten();
	vector <string> getCloneID();
	vector <RawDataofOneArray> selectBlockData(vector <RawDataofOneArray> rawdata, Design des);
	double adjustCy5(ParameterofOneArray parameter);
	double adjustCy3(ParameterofOneArray parameter);
};

struct my_f_params_P22
{
	vector <double> inten1;
	vector <double> inten2;

	double p11;
	double p12;
	double p21;
};

struct my_f_params_P21
{
	vector <double> inten1;
	vector <double> inten2;

	double p11;
	double p12;
	double p22;
};
