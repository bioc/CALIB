#include <vector>
#include <string>

#include "RawDataofOneArray.h"

using namespace std;

struct NDataofOneGene
{
	vector <double> Conc;	
};

class NormalizedData
{
	vector <NDataofOneGene> NData;
	vector <string> CloneID;
	vector <string> initialCloneID;
	vector <int> condition;
	

private:
	double sum(vector <double> a);
	double mean(vector <double> a);
	vector <RawDataofOneArray> pickOneClone(string cloneid,vector <RawDataofOneArray>::iterator prawdata,size_t size);
	void calculateInitialX0ofOne(double inten1,double inten2, ParameterofOneArray par,double &x0Col1, double &x0Col2);
	void calculateIntialX0(vector <RawDataofOneArray>::iterator prawofone,vector <int> raw_id,vector <int> a_unique, vector <int> c_unique,vector <int> c,vector <ParameterofOneArray> par,vector <int> par_id,vector <double> &conc,vector <double> &upper,vector <double> &lower);
	double costFunction(double x,char emodel,struct my_f_params_A *paramsA);
	double calculateCostFunction(double ka,double mus,double inten,double conc1,double conc2,double p1,double p2,double sigmam,double sigmaa,double spoterror,char emodel);
	double spoterrorFunction(double serror,char emodel,struct my_f_params_S *paramsS);
	int calculateSpotErrorofOneSpot(double ka,double mus,double inten1,double conc1,double p11,double p21,double sigmam1,double sigmaa1,double inten2,double conc2,double p12,double p22,double sigmam2,double sigmaa2,double &Qestim,double &spoterror,char emodel);
	double calculateQestim(double inten1,double inten2,double conc1,double conc2,ParameterofOneArray parofone,char emodel);
	double calculateAllCost(vector <RawDataofOneArray>::iterator prawofone,size_t size,vector <int> c_unique,vector <int> a,vector <int> c,vector <ParameterofOneArray> par,vector <int> par_id,vector <double> x0,char emodel);
	vector <double> normalizeOneClone(vector <RawDataofOneArray>::iterator prawofone,size_t sizeofone,vector <int> raw_id,vector <int> a_unique, vector <int> c_unique,vector <int> a,vector <int> c,vector <ParameterofOneArray> par,vector <int> par_id,char emodel);

public:
	void setCloneID(RawDataofOneArray rawdata);
	vector <NDataofOneGene> getData();
	vector <string> getCloneID();
	vector <int> getCondition();
	void normalizeAllClone(vector <RawDataofOneArray>::iterator prawdata,size_t size,Design des,vector <ParameterofOneArray> par,char emodel);
};
