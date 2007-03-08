#include <cmath>
#include <algorithm>

#include "NormalizedData.h"

#define SQRT2 1.414
#define MAX_ITERATION_COST 10
#define MAX_ITERATION_SPOT 5
#define	CONVERGENCE_ABS 0.01
#define Widen_Factor1 1.1
#define Widen_Factor2 0.9
#define GOLDEN_SECTION 0.618


// =================================== private member function =========================================

double NormalizedData::sum(vector <double> a)
{
	double result = 0;

	vector <double>::iterator p = a.begin();

	while(p!= a.end())
	{
		result += *p;
		p++;
	}

	return result;
}

double NormalizedData::mean(vector <double> a)
{
	double result = 0;
	double sumofa = sum(a);
	double size = a.size();

	result = sumofa/size;

	return result;
}

vector <RawDataofOneArray> NormalizedData::pickOneClone(string cloneid,vector <RawDataofOneArray>::iterator prawdata,size_t size) 
{	
	// the index of the "cloneid"

	unsigned int i = 0;	

	vector <int> index;
	vector <int>::iterator pindex; 

	for (i=0;i<initialCloneID.size();i++)
	{
		if (initialCloneID[i] == cloneid)
		{
			index.push_back(i); 
		}
	}

	RawDataofOneArray rawdataofone;
	int arrayid;

	vector <double> inten1;
	vector <double> inten2;
	
	vector <RawDataofOneArray> rawdata_pick;

	for (i=0;i<size;i++)
	{
		rawdataofone = *prawdata;
		inten1 = rawdataofone.getCol1Inten(); 
		inten2 = rawdataofone.getCol2Inten();
		arrayid = rawdataofone.getArrayID(); 

		vector <double> inten1_pick;
		vector <double> inten2_pick;
		
		RawDataofOneArray rawdataofone_pick;

		pindex = index.begin();

		while (pindex != index.end())
		{
			inten1_pick.push_back(inten1[*pindex]);
			inten2_pick.push_back(inten2[*pindex]);

			pindex++;
		}

		rawdataofone_pick.setArrayID(arrayid);
		rawdataofone_pick.setRawData(inten1_pick,inten2_pick); 
		rawdata_pick.push_back(rawdataofone_pick);

		prawdata++;
	}

	return (rawdata_pick);
}

void NormalizedData::calculateInitialX0ofOne(double inten1,double inten2, ParameterofOneArray par,double &x0Col1, double &x0Col2)
{
	double ka = par.getKa();
	double mus = par.getMuSpot();
	double p11 = par.getP1Col1();
	double p12 = par.getP1Col2();
	double p21 = par.getP2Col1();
	double p22 = par.getP2Col2();

	x0Col1 = ((inten1 - p21) * 1/ka)/(p11 * mus - (inten2 - p22) * p11/p12 - (inten1 - p21));
	x0Col2 = ((inten2 - p22) * 1/ka)/(p12 * mus - (inten1 - p21) * p12/p11 - (inten2 - p22));
}

void NormalizedData::calculateIntialX0(vector<RawDataofOneArray>::iterator prawofone,vector <int> raw_id,vector <int> a_unique,vector <int> c_unique, vector <int> c, vector<ParameterofOneArray> par,vector <int> par_id,vector <double> &conc,vector <double> &upper,vector <double> &lower) 
{
	vector <int>::iterator pa_unique = a_unique.begin();

	ParameterofOneArray parameterofone;

	RawDataofOneArray rawdataofone;
	
	unsigned int i = 0;
	int index = 0;
	vector <double> inten1;
	vector <double> inten2;
	vector <double>::iterator pinten1;
	vector <double>::iterator pinten2;
	double x0col1 = 0;
	double x0col2 = 0;
	vector <double> allconc;
	vector <double> allupper;
	vector <double> alllower;

	while (pa_unique != a_unique.end())
	{
		for (i=0;i<par_id.size();i++)
		{
			if (par_id[i] == *pa_unique)
			{
				index = i;
			}
		}

		parameterofone = par[index];

		for (i=0;i<raw_id.size();i++)
		{
			if (raw_id[i] == *pa_unique)
			{
				index = i;
			}
		}

		rawdataofone = *(prawofone + index);

		inten1 = rawdataofone.getCol1Inten();
		inten2 = rawdataofone.getCol2Inten();

		pinten1 = inten1.begin();
		pinten2 = inten2.begin();

		vector <double> x0col1_rep;
		vector <double> x0col2_rep;

		while (pinten1!=inten1.end())
		{
			calculateInitialX0ofOne(*pinten1,*pinten2,parameterofone,x0col1,x0col2);
			
			// concentration should be bigger than 0.

			if (x0col1 < 0)
			{
				x0col1_rep.push_back(0);
			}else
			{
				x0col1_rep.push_back(x0col1);
			};

			
			if (x0col2 < 0)
			{
				x0col2_rep.push_back(0);
			}else
			{
				x0col2_rep.push_back(x0col2);
			};

			pinten1++;
			pinten2++;
		}

		allconc.push_back(mean(x0col1_rep));
		allconc.push_back(mean(x0col2_rep));

		vector <double>::iterator pupper1 = max_element(x0col1_rep.begin(),x0col1_rep.end()); // can be changed here
		vector <double>::iterator pupper2 = max_element(x0col2_rep.begin(),x0col2_rep.end());
	
		allupper.push_back(*pupper1);
		allupper.push_back(*pupper2);

		vector <double>::iterator plower1 = min_element(x0col1_rep.begin(),x0col1_rep.end());
		vector <double>::iterator plower2 = min_element(x0col2_rep.begin(),x0col2_rep.end());

		alllower.push_back(*plower1);
		alllower.push_back(*plower2); 

		pa_unique++;
	}

	vector <int>::iterator pc_unique = c_unique.begin();	
	vector <int> cond;
	
	while (pc_unique != c_unique.end())
	{
		vector <double> tmpx0;
		vector <double> tmpupper;
		vector <double> tmplower;

		for (i=0;i<c.size();i++)
		{
			if(c[i] == *pc_unique)
			{
				tmpx0.push_back(allconc[i]); 
				tmpupper.push_back(allupper[i]);
				tmplower.push_back(alllower[i]); 
			}
		}

		conc.push_back(mean(tmpx0));

		upper.push_back(*max_element(tmpupper.begin(),tmpupper.end()));
		lower.push_back(*min_element(tmplower.begin(),tmplower.end()));

		pc_unique++;
	}

}

double NormalizedData::costFunction(double x,char emodel,struct my_f_params_A *paramsA)
{
	double pka = (paramsA->ka);
	double pmus = (paramsA->mus);
	double pinten = (paramsA->inten);
	double pconc1 = (paramsA->conc1);
	double pconc2 = (paramsA->conc2);
	double pp1 = (paramsA->p1);
	double pp2 = (paramsA->p2);
	double psigmam = (paramsA->sigmam);
	double psigmaa = (paramsA->sigmaa);
	double pspoterror = (paramsA->spoterror);

	double xs = 0;
	
	switch (emodel){
		case 'A':
			xs = (pconc1 * (pmus + pspoterror))/(1/pka + pconc1 + pconc2);
			break;
		case 'M':
			xs = (pconc1 * pmus * exp(pspoterror))/(1/pka + pconc1 + pconc2);
			break;
	}
	
	double result = (((log(pinten - x) - log(pp1 * xs + pp2))/(psigmam * SQRT2)) * ((log(pinten - x) - log(pp1 * xs + pp2))/(psigmam * SQRT2))) + (x/(psigmaa * SQRT2))*(x/(psigmaa * SQRT2));

	return result;

}

double NormalizedData::calculateCostFunction(double ka,double mus,double inten,double conc1,double conc2,double p1,double p2,double sigmam,double sigmaa,double spoterror,char emodel)
{
	struct my_f_params_A paramsA = {ka,mus,inten,conc1,conc2,p1,p2,sigmam,sigmaa,spoterror};
	
	/* calculate the boundaries for minimization */

	double 	lowerb = 0;
	double	upperb = 0;

	switch (emodel){
		case 'A':
			upperb = inten - p1 * (conc1 * (mus + spoterror)/(1/ka + conc1 + conc2)) - p2;
			break;
		case 'M':
			upperb = inten - p1 * (conc1 * (mus * exp(spoterror))/(1/ka + conc1 + conc2)) - p2;
			break;
	}		 
		
	if (lowerb > upperb)
	{
		double tmp = 0;

		tmp = lowerb;
		lowerb = upperb;
		upperb = tmp;
	}

	double iter = 0;
	double a = lowerb;
	double b = upperb;
	double x1 = a + (1 - GOLDEN_SECTION) * (b - a);
	double x2 = a + GOLDEN_SECTION * (b - a);
	double f1 = costFunction(x1,emodel,&paramsA);
	double f2 = costFunction(x2,emodel,&paramsA);

	while ((abs(b-a) > CONVERGENCE_ABS) && (iter < MAX_ITERATION_COST))
	{
		if (f1 < f2)
		{
			b = x2;
			x2 = x1;
			f2 = f1;
			x1 = a + (1 - GOLDEN_SECTION) * (b - a);
			f1 = costFunction(x1,emodel,&paramsA);

		}else{
			if (f1 == f2)
			{
				a = x1;
				b = x2;
				x1 = a + (1 - GOLDEN_SECTION) * (b - a);
				x2 = a + GOLDEN_SECTION * (b - a);			
				f1 = costFunction(x1,emodel,&paramsA);
				f2 = costFunction(x2,emodel,&paramsA);
			}else{
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + GOLDEN_SECTION * (b - a);
				f2 = costFunction(x2,emodel,&paramsA);
			}
		}

		iter++;

	}

	double min = (a+b)/2;
	double fmin = costFunction(min,emodel,&paramsA);
	
	return 	fmin;
}

double NormalizedData::spoterrorFunction(double serror,char emodel,struct my_f_params_S *paramsS)
{
	double pka = (paramsS->ka);
	double pmus = (paramsS->mus);

	double pinten1 = (paramsS->inten1);
	double pconc1 = (paramsS->conc1);
	double pconc2 = (paramsS->conc2);
	double pp11 = (paramsS->p11);
	double pp21 = (paramsS->p21);
	double psigmam1 = (paramsS->sigmam1);
	double psigmaa1 = (paramsS->sigmaa1);

	double pinten2 = (paramsS->inten2);
	double pp12 = (paramsS->p12);
	double pp22 = (paramsS->p22);
	double psigmam2 = (paramsS->sigmam2);
	double psigmaa2 = (paramsS->sigmaa2);

	double result = calculateCostFunction(pka,pmus,pinten1,pconc1,pconc2,pp11,pp21,psigmam1,psigmaa1,serror,emodel) + calculateCostFunction(pka,pmus,pinten2,pconc2,pconc1,pp12,pp22,psigmam2,psigmaa2,serror,emodel);

	return result;
}


int NormalizedData::calculateSpotErrorofOneSpot(double ka,double mus,double inten1,double conc1,double p11,double p21,double sigmam1,double sigmaa1,double inten2,double conc2,double p12,double p22,double sigmam2,double sigmaa2,double &Qestim,double &spoterror,char emodel) //? the paremeters you wen ti?
{
	/* give all the parameters to the function we want to minimize */
	
	struct my_f_params_S paramsS = {ka,mus,inten1,conc1,p11,p21,sigmam1,sigmaa1,inten2,conc2,p12,p22,sigmam2,sigmaa2}; 
		
	/* calculate the boundaries for minimization */ 
	
	if (conc1 == 0 && conc2 == 0) 
	{
		spoterror = 0;
		Qestim = spoterrorFunction(spoterror,emodel,&paramsS);
		return(0);
	}

	double lowerb = 0;
	double upperb = 0;

	double xs1 = 0;
	double xs2 = 0;

	double s1 = 0;
	double s2 = 0;

	if (conc1 != 0)
	{
		xs1 = (inten1 - p21)/p11;
		
		switch(emodel){
			case 'A':
				if (xs1 <= 0) 
				{
					s1 = -mus;
				}else
				{
					s1 = ((xs1 * (1/ka + conc1 + conc2))/conc1 - mus)/mus;
				}
				break;
			case 'M':
				if (xs1 <= 0) 
				{
					double RATIO_MINS0_MUS =10000000;
					s1 = -log(RATIO_MINS0_MUS);
				}else 
				{
					s1 = log(xs1 * (1/ka + conc1 + conc2)/(conc1 * mus));
				}
				break;
		}
	}

	if (conc2 != 0)
	{
		xs2 = (inten2 - p22)/p12;
		switch(emodel){
			case 'A':
				if (xs2 <= 0)
				{
					s2 = -mus;
				}else
				{
					s2 = ((xs2 * (1/ka + conc1 + conc2))/conc2 - mus)/mus;
				}
				break;
			case 'M':
				if (xs2 <= 0)
				{
					double RATIO_MINS0_MUS =10000000;
					s2 = -log(RATIO_MINS0_MUS);
				}else
				{
					s2 = log(xs2 * (1/ka + conc1 + conc2)/(conc2 * mus));
				}
				break;
		}
	}


	// widen the boundaries, for minimizing the Q_a1 and Q_a2.

	if (s1 > s2)
	{
		if (s2 < 0)
		{
			//lowerb = Widen_Factor1 * s2;
			lowerb = s2;
		}else
		{
			//lowerb = Widen_Factor2 * s2;
			lowerb = 0;
		}

		if (s1 < 0)
		{
			//upperb = Widen_Factor2 * s1;
			upperb = 0;
		}else
		{
			//upperb = Widen_Factor1 * s1;
			upperb = s1;
		}
	}else
	{
		if (s1 < 0)
		{
			//lowerb = Widen_Factor1 * s1;
			lowerb = s1;
		}else
		{
			//lowerb = Widen_Factor2 * s1;
			lowerb = 0;
		}

		if (s2 < 0)
		{
			//upperb = Widen_Factor2 * s2;
			upperb = 0;
		}else
		{
			//upperb = Widen_Factor1 * s2;
			upperb = s2;
		}
		
	}

	///* make a vague guess of the initial minimal value */

	double iter = 0;
	double a = lowerb;
	double b = upperb;
	double x1 = a + (1 - GOLDEN_SECTION) * (b - a);
	double x2 = a + GOLDEN_SECTION * (b - a);
	double f1 = spoterrorFunction(x1,emodel,&paramsS);
	double f2 = spoterrorFunction(x2,emodel,&paramsS);

	while ((abs(b-a) > CONVERGENCE_ABS) && (iter < MAX_ITERATION_SPOT))
	{
		if (f1 < f2)
		{
			b = x2;
			x2 = x1;
			f2 = f1;
			x1 = a + (1 - GOLDEN_SECTION) * (b - a);
			f1 = spoterrorFunction(x1,emodel,&paramsS);

		}else{
			if (f1 == f2)
			{
				a = x1;
				b = x2;
				x1 = a + (1 - GOLDEN_SECTION) * (b - a);
				x2 = a + GOLDEN_SECTION * (b - a);			
				f1 = spoterrorFunction(x1,emodel,&paramsS);
				f2 = spoterrorFunction(x2,emodel,&paramsS);
			}else{
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + GOLDEN_SECTION * (b - a);
				f2 = spoterrorFunction(x2,emodel,&paramsS);
			}
		}

		iter++;

	}

	double min = (a+b)/2;
	double fmin = spoterrorFunction(min,emodel,&paramsS);

	spoterror = min;
	Qestim = fmin;

	return(0);

}

double NormalizedData::calculateQestim(double inten1,double inten2,double conc1,double conc2,ParameterofOneArray parofone,char emodel)
{
	double ka = parofone.getKa();
	double mus = parofone.getMuSpot();
	double p11 = parofone.getP1Col1();
	double p12 = parofone.getP1Col2();
	double p21 = parofone.getP2Col1();
	double p22 = parofone.getP2Col2();
	double sigmaa1 = parofone.getSigmaAddCol1();
	double sigmaa2 = parofone.getSigmaAddCol2();
	double sigmam1 = parofone.getSigmaMulCol1();
	double sigmam2 = parofone.getSigmaMulCol2();
	double sigmas = parofone.getSigmaSpot();

	double Qestim = 0;
	double spoterror = 0;

	calculateSpotErrorofOneSpot(ka,mus,inten1,conc1,p11,p21,sigmam1,sigmaa1,inten2,conc2,p12,p22,sigmam2,sigmaa2,Qestim,spoterror,emodel);

	double Q = 0;

	//Q = Qestim + 2 * ((spoterror/(sigmas * SQRT2)) * (spoterror/(sigmas * SQRT2))); // simply version as the following;
	Q = Qestim + (spoterror/sigmas) * (spoterror/sigmas);

	return Q;	
}


double NormalizedData::calculateAllCost(vector<RawDataofOneArray>::iterator prawofone,size_t size, vector <int> c_unique, vector <int> a, vector <int> c,vector<ParameterofOneArray> par,vector <int> par_id,vector<double> x0,char emodel) 
{
	// get the array id involved in this block.
	// the order of array id is corresponded with the order of x0.

	vector <int>::iterator pc_unique = c_unique.begin();

	vector <double>::iterator pconc = x0.begin();
	vector <double> allconc(a.size());

	unsigned int i = 0;

	while (pc_unique != c_unique.end())
	{
		for (i=0;i<c.size();i++)
		{
			if(c[i] == *pc_unique)
			{
				allconc[i] = *pconc;
			}
		}

		pc_unique++;
		pconc++;
	}

	// variables definition.

	RawDataofOneArray rawdataofone;
	ParameterofOneArray parameterofone;
	
	int arrayid = 0;

	vector <double> vinten1;
	vector <double> vinten2;
	double inten1 = 0;
	double inten2 = 0;
	vector <double> conc;
	double conc1 = 0;
	double conc2 = 0;

	double Qestim = 0;
	vector <double> Q;
	double Qnorm = 0;

	unsigned int index = 0;

	// in each iteration (prawdata), calculate the Q value of one array.

	for (index=0;index<size;index++)
	{
		rawdataofone = *prawofone;

		vinten1 = rawdataofone.getCol1Inten();
		vinten2 = rawdataofone.getCol2Inten();
		arrayid = rawdataofone.getArrayID(); 

		// get parameter.
		for (i=0;i<par_id.size();i++)
		{
			if (par_id[i] == arrayid)
			{
				parameterofone = par[i];
			}
		}

		// get concentration.
		for (i=0;i<a.size();i++)
		{
			if (a[i] == arrayid)
			{
				conc.push_back(allconc[i]); 
			}
		}

		conc1 = conc[0];
		conc2 = conc[1];	

		// get intensity. 
		// each iteration i, get value for one replicate.

		for (i=0;i<vinten1.size();i++)
		{
			inten1 = vinten1[i];
			inten2 = vinten2[i];
			Qestim = calculateQestim(inten1,inten2,conc1,conc2,parameterofone,emodel);
			Q.push_back(Qestim);
		}

		prawofone++;
		conc.clear();
	}

	Qnorm = sum(Q);

	return (Qnorm);
}

vector <double> NormalizedData::normalizeOneClone(vector<RawDataofOneArray>::iterator prawofone,size_t sizeofone, vector <int> raw_id, vector <int> a_unique,vector <int> c_unique,vector <int> a,vector <int> c, vector<ParameterofOneArray> par,vector <int> par_id,char emodel)
{
	vector <double> conc;
	vector <double> upper;
	vector <double> lower;

	// initial concentration value
	calculateIntialX0(prawofone,raw_id,a_unique,c_unique,c,par,par_id,conc,upper,lower);
	
	unsigned int i = 0;
	unsigned int numberofzero = 0;

	vector <double>::iterator pconc = conc.begin();

	while (pconc != conc.end())
	{
		if (*pconc == 0)
		{
			numberofzero++;
		}

		pconc++;
	}

	if (numberofzero == conc.size())
	{
		return (conc);

	}else
	{
		double Q_initial = 0;
		double Q = 0;
		double Q_direction = 0;

		// initial cost function value
		Q_initial = calculateAllCost(prawofone,sizeofone,c_unique,a,c,par,par_id,conc,emodel);
		Q = Q_initial;

		double Q_pre = 2 * Q;

		double CONVERGE_FACTOR = 0.005;
		
		vector <double> conc_new;

		int iteration = 0;
		int maxiter = 10;
	
		//while ((Q_pre - Q) > (Q_initial * CONVERGE_FACTOR)) 
		while ((abs(Q_pre - Q) > (Q_initial * CONVERGE_FACTOR)) && (iteration < maxiter))
		//while (abs(Q_pre - Q) > (Q_initial * CONVERGE_FACTOR))
		{
			Q_pre = Q;
	
			// decide the direction of moving:

			for (i=0;i<conc.size();i++)
			{
				double DELTA = (upper[i] - lower[i])/100000;
				conc_new = conc;

				if (conc_new[i] != 0)
				{
					conc_new[i] = conc[i] * (1 + DELTA);
					Q_direction = calculateAllCost(prawofone,sizeofone,c_unique,a,c,par,par_id,conc_new,emodel);
										
					if (Q_direction < Q_pre)
					{
						lower[i] = conc[i];
					}else
					{
						upper[i] = conc[i];
					}
				}
			}

			// after for, we have new upper bound and lower bound.
			// now take the mean of the upper and lower as the new concentration.

			vector <double>::iterator pupper = upper.begin();
			vector <double>::iterator plower = lower.begin();
			double tmp_mean = 0;

			conc_new.clear();

			while (pupper != upper.end())
			{
				tmp_mean = ((*pupper) + (*plower))/2;
				conc_new.push_back(tmp_mean); 

				pupper++;
				plower++;
			}

			Q = calculateAllCost(prawofone,sizeofone,c_unique,a,c,par,par_id,conc_new,emodel);
	
			conc = conc_new;

			iteration++;
		}

		return (conc);
	}
}

// ========================================= public member function ==================================

void NormalizedData::normalizeAllClone(vector<RawDataofOneArray>::iterator prawdata,size_t size, Design des, vector <ParameterofOneArray> par,char emodel) 
{
	// design vectors:
	
	vector <int> a;
	vector <int> c;

	a = des.getArray();
	c = des.getCond();

	vector <int> a_unique(a.begin(),a.end());
	sort(a_unique.begin(),a_unique.end());
	vector <int>::iterator pa_unique = unique(a_unique.begin(),a_unique.end());
	a_unique.erase(pa_unique,a_unique.end());

	vector <int> c_unique(c.begin(),c.end());
	sort(c_unique.begin(),c_unique.end());
	vector <int>::iterator pc_unique = unique(c_unique.begin(),c_unique.end());
	c_unique.erase(pc_unique,c_unique.end());

	// end of design vectors

	// pick out the parameter involved in this block (parameter id)

	ParameterofOneArray parameterofone;

	vector <ParameterofOneArray>::iterator ppar = par.begin();
	vector <int> par_id;

	while (ppar != par.end())
	{
		parameterofone = *ppar;
		par_id.push_back(parameterofone.getArrayID());

		ppar++;
	}

	// pick out raw data id

	RawDataofOneArray rawdataofone;

	vector <RawDataofOneArray>::iterator praw = prawdata;
	vector <int> raw_id;

	unsigned int index = 0;

	for (index=0;index<size;index++)
	{
		rawdataofone = *praw;
		raw_id.push_back(rawdataofone.getArrayID());

		praw++;
	}

	rawdataofone = *prawdata;

	vector <string>::iterator pid = CloneID.begin();

	vector <RawDataofOneArray> rawdataofoneclone;
	vector <double> conc;

	//int i = 0;

	//for(i=0;i<10;i++)
	
	while (pid != CloneID.end())
	{

		//rawdataofoneclone = pickOneClone(CloneID[i],prawdata,size);
		rawdataofoneclone = pickOneClone(*pid,prawdata,size);
		
		vector <RawDataofOneArray>::iterator prawofone = rawdataofoneclone.begin();
		size_t sizeofone = rawdataofoneclone.size();
	
		conc = normalizeOneClone(prawofone,sizeofone,raw_id,a_unique,c_unique,a,c,par,par_id,emodel);
				
		NDataofOneGene ndataofone;

		ndataofone.Conc = conc;

		NData.push_back(ndataofone); 

		pid++;
	}

	condition = c_unique;
}

void NormalizedData::setCloneID(RawDataofOneArray rawdata)
{
	initialCloneID = rawdata.getCloneID();
	CloneID = initialCloneID;

	sort(CloneID.begin(),CloneID.end());
	vector <string>::iterator pid = unique(CloneID.begin(),CloneID.end());
	CloneID.erase(pid,CloneID.end());
}

vector <NDataofOneGene> NormalizedData::getData() 
{
	return NData;
}

vector <string> NormalizedData::getCloneID()
{
	return CloneID;
}

vector <int> NormalizedData::getCondition()
{
	return condition;
}
