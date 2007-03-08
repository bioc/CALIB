#include <vector>
#include <string>

using namespace std;

struct Spike
{
	vector <double> Col1Conc;
	vector <double> Col2Conc;
	vector <double> Col1Inten;
	vector <double> Col2Inten;
	vector <string> SpotID;
	vector <string> Type;
};

class SpikeofOneArray
{
	int ArrayID;
	Spike spikes;

private:
	double mean(vector <double> a);
	double standarddeviation(vector <double> a);
	double calculatestd(vector <double> X,vector <double> Y,double meanX,double meanY);

public:
	void setArrayID(int id);
	void setSpikes(vector <double> SCol1Conc,vector <double> SCol2Conc,vector <double> SCol1Inten,vector <double> SCol2Inten, vector <string> SpikeType);
	Spike selectSpikes(string selecttype);
	Spike kickoutNegatives(string kickouttype);
	double calculateAdditiveMean(vector <double> a);
	double calculateAdditiveVariance(vector <double> a);
	double calculateMultiplicativeVariance(vector <double> a,vector <double> b);
	int getArrayID();
	vector <double> getCol1Conc();
	vector <double> getCol1Conc(Spike s);
	vector <double> getCol2Conc();
	vector <double> getCol2Conc(Spike s);
	vector <double> getCol1Inten();
	vector <double> getCol1Inten(Spike s);
	vector <double> getCol2Inten();
    vector <double> getCol2Inten(Spike s);
	vector <string> getSpotID();
	vector <string> getSpotID(Spike s);
	vector <string> getType();
	vector <string> getType(Spike s);	
};
