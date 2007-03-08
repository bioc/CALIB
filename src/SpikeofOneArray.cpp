#include <cmath>
#include <fstream>

#include "SpikeofOneArray.h"

#define PI 3.14159265

// member function of class SpikeofOneArray

double SpikeofOneArray::mean(vector <double> a)
{
	double result = 0;
	double sum = 0;
	int size = 0;
	vector <double>::iterator p = a.begin();

	while(p!= a.end())
	{
		sum += *p;
		size++;
		p++;
	}

	result = sum/size;

	return result;
}

double SpikeofOneArray::standarddeviation(vector <double> a)
{
	double result = 0;
	double avg = mean(a);
	double sumsquare = 0;
	double tmp = 0;
	int size = 0;
	vector <double>::iterator p=a.begin();

	while(p!= a.end())
	{
		tmp = (*p-avg)*(*p-avg);
		sumsquare += tmp;
		size++;
		p++;
	}

	result = sqrt(sumsquare/(size-1));

	return result;
}
			
double SpikeofOneArray::calculatestd(vector <double> X,vector <double> Y,double meanX,double meanY) 
{
	// regression of Y~aX+b; and also X~a'Y+b'

	vector <double>::iterator px = X.begin();
	vector <double>::iterator py = Y.begin();

	double aYX,aXY,bYX,bXY;
	double numerator, denomYX,denomXY;

	aYX = aXY = bYX = bXY = 0;
	numerator = denomYX = denomXY = 0;


	while (px != X.end())
	{
		numerator += (*px - meanX)*(*py - meanY); 
		denomYX += (*px - meanX)*(*px - meanX);
		denomXY += (*py - meanY)*(*py - meanY);

		px++;
		py++;
	}

	aYX = numerator/denomYX;
	bYX = meanY - aYX * meanX;

	aXY = numerator/denomXY;
	bXY = meanX - aXY * meanY;

	// approximation of orthognal regression

	double a,b;
	a = b = 0;

	a = (aYX*aXY*bXY + aXY*bYX + bXY + aYX*aXY*aXY*bYX)/(2*aXY*aXY*bYX + 2*aXY*bXY);
	b = (aXY*bYX - bXY)/(2*aXY);

	double sse = 0;
	px = X.begin();
	py = Y.begin();

	while(px != X.end())
	{
		sse += ((*px)*(*px)+(*py - b)*(*py - b))*(sin(PI/2 - atan(a) - atan((*px)/(*py - b))))*(sin(PI/2 - atan(a) - atan((*px)/(*py - b))));

		px++;
		py++;

	}

	double standd = sqrt(sse);

	return standd;
}

void SpikeofOneArray::setArrayID(int id)
{
	ArrayID = id;
}

void SpikeofOneArray::setSpikes(vector <double> SCol1Conc,vector <double> SCol2Conc,vector <double> SCol1Inten,vector <double> SCol2Inten,vector <string> SpikeType)
{
	spikes.Col1Conc = SCol1Conc;
	spikes.Col2Conc = SCol2Conc;
	spikes.Col1Inten = SCol1Inten;
	spikes.Col2Inten = SCol2Inten;
	//spikes.SpotID = SpotID;
	spikes.Type = SpikeType;
}

int SpikeofOneArray::getArrayID()
{
	return ArrayID;
}

vector <double> SpikeofOneArray::getCol1Conc()
{
	return spikes.Col1Conc;
}

vector <double> SpikeofOneArray::getCol1Conc(Spike s)
{
	return s.Col1Conc; 
}

vector <double> SpikeofOneArray::getCol2Conc()
{
	return spikes.Col2Conc;
}

vector <double> SpikeofOneArray::getCol2Conc(Spike s)
{
	return s.Col2Conc;
}

vector <double> SpikeofOneArray::getCol1Inten()
{
	return spikes.Col1Inten;
}

vector <double> SpikeofOneArray::getCol1Inten(Spike s)
{
	return s.Col1Inten;
}

vector <double> SpikeofOneArray::getCol2Inten()
{
	return spikes.Col2Inten;
}

vector <double> SpikeofOneArray::getCol2Inten(Spike s)
{
	return s.Col2Inten;
}

vector <string> SpikeofOneArray::getSpotID()
{
	return spikes.SpotID;
}

vector <string> SpikeofOneArray::getSpotID(Spike s)
{
	return s.SpotID;
}

vector <string> SpikeofOneArray::getType()
{
	return spikes.Type;
}

vector <string> SpikeofOneArray::getType(Spike s)
{
	return s.Type;
}

Spike SpikeofOneArray::selectSpikes(string selecttype) 
{
	Spike selected;
		
	vector <double>::iterator p1conc = spikes.Col1Conc.begin();
	vector <double>::iterator p2conc = spikes.Col2Conc.begin();
	vector <double>::iterator p1inten = spikes.Col1Inten.begin();
	vector <double>::iterator p2inten = spikes.Col2Inten.begin();
	//vector <string>::iterator pid = spikes.SpotID.begin();
	vector <string>::iterator ptype = spikes.Type.begin();

	while(ptype!= spikes.Type.end())  //didn't consider missing value here
	{
		if(*ptype==selecttype)
		{
			selected.Col1Conc.push_back(*p1conc);
			selected.Col2Conc.push_back(*p2conc);
			selected.Col1Inten.push_back(*p1inten);
			selected.Col2Inten.push_back(*p2inten);
			//selected.SpotID.push_back(*pid);
			selected.Type.push_back(*ptype);
		}

		p1conc++;
		p2conc++;
		p1inten++;
		p2inten++;
		//pid++;
		ptype++;
	}

	return selected;
}

Spike SpikeofOneArray::kickoutNegatives(string kickouttype)
{
	Spike selected;
		
	vector <double>::iterator p1conc = spikes.Col1Conc.begin();
	vector <double>::iterator p2conc = spikes.Col2Conc.begin();
	vector <double>::iterator p1inten = spikes.Col1Inten.begin();
	vector <double>::iterator p2inten = spikes.Col2Inten.begin();
	//vector <string>::iterator pid = spikes.SpotID.begin();
	vector <string>::iterator ptype = spikes.Type.begin();

	while(ptype!= spikes.Type.end())  //didn't consider missing value here
	{                                 // missing value could lead unequal-size vectors 
		if(*ptype!=kickouttype)
		{
			selected.Col1Conc.push_back(*p1conc);
			selected.Col2Conc.push_back(*p2conc);
			selected.Col1Inten.push_back(*p1inten);
			selected.Col2Inten.push_back(*p2inten);
			//selected.SpotID.push_back(*pid);
			selected.Type.push_back(*ptype);
		}

		p1conc++;
		p2conc++;
		p1inten++;
		p2inten++;
		//pid++;
		ptype++;
	}

	return selected;
}

double SpikeofOneArray::calculateAdditiveMean(vector <double> a)
{
	double m = 0;

	m = mean(a);

	return m;
}

double SpikeofOneArray::calculateAdditiveVariance(vector<double> a) 
{
	double var = 0;

	var=standarddeviation(a);
	
	return var;
}

double SpikeofOneArray::calculateMultiplicativeVariance(vector <double> a,vector <double> b)
{
	double var = 0;         // the multiplicative variance to return.

	vector <double> inta;     // get the value from vector a, take the log;
	vector <double> intb;     // get the value from vector b, take the log;
	vector <double> distoo;    // calculate the distance to the orign;

	vector <double>::iterator pa = a.begin();
	vector <double>::iterator pb = b.begin();

	double tmpa = 0;
	double tmpb = 0;
	double tmpdis = 0;

	while(pa!=a.end())
	{
		tmpa = log(*pa);
		tmpb = log(*pb);
		tmpdis = sqrt (tmpa * tmpa + tmpb * tmpb);

		inta.push_back(tmpa);
		intb.push_back(tmpb);
		distoo.push_back(tmpdis); 

		pa++;
		pb++;
	}
	
	/* sort (from max to min) inta and intb according to the distance to the orign,
	the descending order is for using pop_back.*/

	unsigned int i = 0;
	unsigned int j = 0;

	vector <double> sortedinta;
	vector <double> sortedintb;

	vector <double>::iterator pinta = inta.begin();
	vector <double>::iterator pintb = intb.begin();
	vector <double>:: iterator pdist = distoo.begin();

	size_t size = distoo.size();

	for (i=0; i<size; i++) /* in each iteration, we find one value: in 1st iteration, 
								   find the largest value; in 2nd iteration, find the second large one...
								   the i has nothing to do with the ith element in the vector*/

	{
		double max = 0;
		unsigned int index = 0;

		for (j=0;j<distoo.size();j++)  // in each ith iteration, find out the index of the largest value in current distoo vector
		{
			if (distoo.at(j) > max)
			{
				max = distoo.at(j);
				index = j;
			}
		}

		sortedinta.push_back(inta.at(index));  // sort the inta and intb according to the index
		sortedintb.push_back(intb.at(index));

		pdist = distoo.begin();
		pdist += index;
		distoo.erase(pdist); // delete the largest value, go to the next iteration.	

		pinta = inta.begin();
		pinta += index;
		inta.erase(pinta);

		pintb = intb.begin();
		pintb += index;
		intb.erase(pintb);
	}

	double standd = 0;
	vector <double> vstandd; /* vector to store the std, the first element of this vector is
							 the std after dropping the first point, the second one is 
							 the std after dropping the first point and the second point...*/

	double Stable_Factor = 0.05;
	double Variance_Factor = 0.05;

	unsigned int Stable_Size = int(Stable_Factor * a.size()); // need to discuss about 0.05
	bool flag = true;
	
	/* calculate the standard deviation by adding the points one by one
	   from the begining, we take three points in order to make svd thing work.*/

	vector <double> usefulinta;
	vector <double> usefulintb;

	for (i=0;i<2;i++)
	{
		usefulinta.push_back(sortedinta[i]);
		usefulintb.push_back(sortedintb[i]);
	}

	i = 2;

	// pick out the points from the largest distance, find out the stable value.

	
	vector <double> scaledinta;
	vector <double> scaledintb;
	
	while (flag)
	{
		usefulinta.push_back(sortedinta[i]);
		usefulintb.push_back(sortedintb[i]);

		double meaninta = mean(usefulinta);
		double meanintb = mean(usefulintb);
		standd = calculatestd(usefulinta,usefulintb,meaninta,meanintb);
				
		vstandd.push_back(standd); 
		
		/* to decide whether stop dropping points and stop the calculation, it means that this standd is the SigmaMul */


		if (vstandd.size() > Stable_Size)
		{
			vector <double> vtmp(vstandd.end() - Stable_Size - 1,vstandd.end() - 1); /* get the last 20 std (20 indicated by Stable_Size), take the mean of these 20 std,
								  if the difference between the current std --- standd 
								  and the mean is smaller than the threshold, stop calculation.*/
			
			double tmpmean = 0;
			tmpmean = mean(vtmp);

			if (fabs(tmpmean - standd) < Variance_Factor * standd )  // need to discuss about 0.005
			{ 
				flag = false;
				var = tmpmean;
			}
		}

		i++;
	}

	//var = var/10;

	return var;
}

