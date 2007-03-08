#include <vector>
#include <string>

using namespace std;

class ParameterofOneArray
{
	int ArrayID;
	double Ka;
	double MuSpot;
	double SpotErrorofOneSpot;
	double P1Col1;
	double P1Col2;
	double P2Col1;
	double P2Col2;
	double MuAddCol1;
	double MuAddCol2;
	double SigmaAddCol1;
	double SigmaAddCol2;
	double SigmaMulCol1;
	double SigmaMulCol2;
	double SigmaSpot;
	vector <double> SpotError;
	
private:	
	double mean(vector <double> a);
	double standarddeviation(vector <double> a);
	void calculateMeanPoint(vector <double> &IntenCol1,vector <double> &IntenCol2,vector <double> &ConcCol1,vector <double> &ConcCol2);
	double costFunction(double x,char emodel,struct my_f_params_A *paramsA);
	double calculateCostFunction(double ka,double mus,double inten,double conc1,double conc2,double p1,double p2,double sigmam,double sigmaa,double spoterror,char emodel);
	double spoterrorFunction(double serror,char emodel,struct my_f_params_S *paramsS);
	double calculateSpotErrorofOneSpot(double ka,double mus,double inten1,double conc1,double p11,double p21,double sigmam1,double sigmaa1,double inten2,double conc2,double p12,double p22,double sigmam2,double sigmaa2,char emodel);
	double fineKaFunction(double k,char emodel,struct my_f_params_K *paramsK);
	void calculateXs(vector <double> &conc1,vector <double> &conc2,char emodel);
	double costFunctionGivenXs(double x,struct my_f_params_AXs *paramsAXs);
	double calculateCostFunctionGivenXs(double ka,double mus,double inten,double xs,double p1,double p2,double sigmam,double sigmaa);
	double P1Function(double p1,struct my_f_params_P1 *paramsP1);
	double estimateP1(double ka,double mus,double p2, double sigmam,double sigmaa,vector<double> inten,vector<double> xs,double initialp1);


public:
	void setArrayID(int id);
	void setMuSpot(double MuS);
	void setKa(double ka);
	void setP1Col1(double p11);
	void setP1Col2(double p12);
	void setP2Col1(double p21);
	void setP2Col2(double p22);
	void setSigmaAddCol1(double SigmaA1);
	void setSigmaAddCol2(double SigmaA2);
	void setSigmaMulCol1(double SigmaM1);
	void setSigmaMulCol2(double SigmaM2);
	void setSigmaSpot(double s);
	void setInitialP1(double maxInten);
	void setFineKa(vector <double> IntenCol1,vector <double> IntenCol2,vector <double> ConcCol1,vector <double> ConcCol2,char emodel);
	void calculateP1(vector <double> IntenCol1,vector <double> IntenCol2, vector <double> ConcCol1,vector <double> ConcCol2,char emodel);
	int getArrayID();
	double getKa();
	double getMuSpot();
	double getP1Col1();
	double getP1Col2();
	double getP2Col1();
	double getP2Col2();
	double getMuAddCol1();
	double getMuAddCol2();
	double getSigmaAddCol1();
	double getSigmaAddCol2();
	double getSigmaMulCol1();
	double getSigmaMulCol2();
	vector <double> getSpotError();
	void calculateSigmaSpot();
	double getSigmaSpot();	
};

struct my_f_params_S
{
	// parameters for col1;
	double ka;
	double mus;

	double inten1;
	double conc1;
	double p11;
	double p21;
	double sigmam1;
	double sigmaa1;

	// parameters for col2;
	double inten2;
	double conc2;
	double p12;
	double p22;
	double sigmam2;
	double sigmaa2;
};

struct my_f_params_A
{
	// just one color;
	double ka;
	double mus;

	double inten;
	double conc1;
	double conc2;
	double p1;
	double p2;
	double sigmam;
	double sigmaa;
	double spoterror;
};

struct my_f_params_K
{
	// parameters for col1;
	double mus;

	double p11;
	double p21;
	double sigmam1;
	double sigmaa1;

	// parameters for col2;
	double p12;
	double p22;
	double sigmam2;
	double sigmaa2;

	vector <double> inten1;
	vector <double> inten2;
	vector <double> conc1;
	vector <double> conc2;
};

struct my_f_params_AXs
{
	double ka;
	double mus;

	double inten;
	double xs;
	double p1;
	double p2;
	double sigmam;
	double sigmaa;
};

struct my_f_params_P1
{
	double ka;
	double mus;

	double p2;
	double sigmam;
	double sigmaa;

	vector <double> inten;
	vector <double> xs;
};	
