#include <vector>
#include <algorithm>

using namespace std;

class Design
{
	vector <int> Array;
	vector <int> Cond;
	vector <int> Dye;

public:
	void setArray(vector <int> a);
	void setCond(vector <int> c);
	void setDye(vector <int> d);
	vector <int> getArray();
	vector <int> getCond();
	vector <int> getDye();
	vector <Design> splitBlock(Design d);
};  
