#include <cmath>
#include <algorithm>

#include "ParameterofOneArray.h"

#define SQRT2 1.414
#define Widen_Factor1 1.1
#define Widen_Factor2 0.9
#define MAX_ITERATION 30
#define	CONVERGENCE_ABS 0.001
#define GOLDEN_SECTION 0.618


/* ============================= private functions ================================ */

double ParameterofOneArray::mean(vector <double> a)
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

double ParameterofOneArray::standarddeviation(vector <double> a)
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

void ParameterofOneArray::calculateMeanPoint(vector <double> &IntenCol1, vector <double> &IntenCol2, vector <double> &ConcCol1,  vector <double> &ConcCol2)
{
	vector <double> inten1(IntenCol1);
	vector <double> inten2(IntenCol2);
	vector <double> conc1(ConcCol1);
	vector <double> conc2(ConcCol2);

	// get the distincive concentration values.

	vector <double>::iterator i1;

	sort(ConcCol1.begin(),ConcCol1.end());
	i1 = unique(ConcCol1.begin(),ConcCol1.end());
	ConcCol1.erase(i1,ConcCol1.end());

	vector <double>::iterator i2;

	sort(ConcCol2.begin(),ConcCol2.end());	
	i2 = unique(ConcCol2.begin(),ConcCol2.end());
	ConcCol2.erase(i2,ConcCol2.end());

	// calculate the intensity mean of each concentration.

	IntenCol1.clear();
	IntenCol2.clear();
	
	unsigned int i = 0;
	unsigned int j = 0;
	
	vector <double> inten1tmp;
	
	for (i=0;i<ConcCol1.size();i++)
	{
		for (j=0;j<conc1.size();j++)
		{
			if (conc1[j] == ConcCol1[i])
			{
				inten1tmp.push_back(inten1[j]);
			}
		}

		IntenCol1.push_back(mean(inten1tmp));
	}

	vector <double> inten2tmp;
	
	for (i=0;i<ConcCol2.size();i++)
	{
		for (j=0;j<conc2.size();j++)
		{
			if (conc2[j] == ConcCol2[i])
			{
				inten2tmp.push_back(inten2[j]);
			}
		}

		IntenCol2.push_back(mean(inten2tmp));
	}
}

double ParameterofOneArray::costFunction (double x,char emodel,struct my_f_params_A *paramsA)
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

double ParameterofOneArray::calculateCostFunction(double ka,double mus,double inten,double conc1,double conc2,double p1,double p2,double sigmam,double sigmaa,double spoterror,char emodel)
{
	
	struct my_f_params_A paramsA = {ka,mus,inten,conc1,conc2,p1,p2,sigmam,sigmaa,spoterror};
		
	/* calculate the boundaries for minimization */

	
	double lowerb = 0;
	double upperb = 0;

	if (inten > 0)
	{
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
	}else{
		
		upperb = inten;

		switch (emodel){
			case 'A':
				lowerb = inten - p1 * (conc1 * (mus + spoterror)/(1/ka + conc1 + conc2)) - p2;
				break;
			case 'M':
				lowerb = inten - p1 * (conc1 * (mus * exp(spoterror))/(1/ka + conc1 + conc2)) - p2;
				break;
		}	
	}

	double iter = 0;
	double a = lowerb;
	double b = upperb;
	double x1 = a + (1 - GOLDEN_SECTION) * (b - a);
	double x2 = a + GOLDEN_SECTION * (b - a);
	double f1 = costFunction(x1,emodel,&paramsA);
	double f2 = costFunction(x2,emodel,&paramsA);

	while ((abs(b-a) > CONVERGENCE_ABS) && (iter < MAX_ITERATION))
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

double ParameterofOneArray::spoterrorFunction(double serror,char emodel,struct my_f_params_S *paramsS)
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

double ParameterofOneArray::calculateSpotErrorofOneSpot(double ka,double mus,double inten1,double conc1,double p11,double p21,double sigmam1,double sigmaa1,double inten2,double conc2,double p12,double p22,double sigmam2,double sigmaa2,char emodel) //? the paremeters you wen ti?
{
	/* give all the parameters to the function we want to minimize */
	
	struct my_f_params_S paramsS = {ka,mus,inten1,conc1,p11,p21,sigmam1,sigmaa1,inten2,conc2,p12,p22,sigmam2,sigmaa2}; 
		
	/* calculate the boundaries for minimization */ 
	
	double xs1 = (inten1 - p21)/p11;
	double xs2 = (inten2 - p22)/p12;

	double s1 = 0;
	double s2 = 0;

	double lowerb = 0;
	double upperb = 0;

	switch (emodel){
		case 'A':
			if (xs1 <= 0 && xs2 <= 0)
			{
				return 0;
				//xs1 = 0;
				//xs2 = 0;
			}else 
			{
				if (xs1 > 0 && xs2 <= 0 )
				{
					xs2 = 0;				
				}else if (xs1 <= 0 && xs2 > 0)
				{
					xs1 = 0;				
				}
			}
			s1 = ((xs1 * (1/ka + conc1 + conc2))/conc1 - mus)/mus;
			s2 = ((xs2 * (1/ka + conc1 + conc2))/conc2 - mus)/mus;
			break;

		case 'M':
			if (xs1 <= 0 && xs2 <=0)
			{
				return 0;
			}else
			{
				if (xs1 > 0 && xs2 <= 0)
				{
					xs2 = 1/mus; 
				}else if (xs1 <= 0 && xs2 > 0)
				{
					xs1 = 1/mus;
				}
			}
			s1 = log(xs1 * (1/ka + conc1 + conc2)/(conc1 * mus));
			s2 = log(xs2 * (1/ka + conc1 + conc2)/(conc2 * mus));
			break;
	}

	
	// widen the boundaries, for minimizing the Q_a1 and Q_a2.

	if (s1 > s2)
	{
		if (s2 < 0)
		{
			lowerb = Widen_Factor1 * s2;
		}else
		{
			lowerb = Widen_Factor2 * s2;
		}

		if (s1 < 0)
		{
			upperb = Widen_Factor2 * s1;
		}else
		{
			upperb = Widen_Factor1 * s1;
		}
	}else
	{
		if (s1 < 0)
		{
			lowerb = Widen_Factor1 * s1;
		}else
		{
			lowerb = Widen_Factor2 * s1;
		}

		if (s2 < 0)
		{
			upperb = Widen_Factor2 * s2;
		}else
		{
			upperb = Widen_Factor1 * s2;
		}
		
	}

	double iter = 0;
	double a = lowerb;
	double b = upperb;
	double x1 = a + (1 - GOLDEN_SECTION) * (b - a);
	double x2 = a + GOLDEN_SECTION * (b - a);
	double f1 = spoterrorFunction(x1,emodel,&paramsS);
	double f2 = spoterrorFunction(x2,emodel,&paramsS);

	while ((abs(b-a) > CONVERGENCE_ABS) && (iter < MAX_ITERATION))
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
		
	return 	min;	
}

double ParameterofOneArray::fineKaFunction(double k,char emodel,struct my_f_params_K *paramsK) 
{
	double pmus = (paramsK->mus);

	double pp11 = (paramsK->p11);
	double pp21 = (paramsK->p21);
	double psigmam1 = (paramsK->sigmam1);
	double psigmaa1 = (paramsK->sigmaa1);

	double pp12 = (paramsK->p12);
	double pp22 = (paramsK->p22);
	double psigmam2 = (paramsK->sigmam2);
	double psigmaa2 = (paramsK->sigmaa2);

	vector <double> inten1 = (paramsK->inten1);
	vector <double> inten2 = (paramsK->inten2);
	vector <double> conc1 = (paramsK->conc1);
	vector <double> conc2 = (paramsK->conc2);


	vector <double>::iterator pinten1 = inten1.begin();
	vector <double>::iterator pinten2 = inten2.begin();
	vector <double>::iterator pconc1 = conc1.begin();
	vector <double>::iterator pconc2 = conc2.begin();

	double ka = exp(k);

	double e = 0;
	double sse = 0;
	
	SpotError.clear();

	while (pinten1!=inten1.end())
	{
		e =calculateSpotErrorofOneSpot(ka,pmus,*pinten1,*pconc1,pp11,pp21,psigmam1,psigmaa1,*pinten2,*pconc2,pp12,pp22,psigmam2,psigmaa2,emodel);
		sse += (e * e);

		SpotError.push_back(e);

		pinten1++;
		pinten2++;
		pconc1++;
		pconc2++;
	}

	return sse;
}

void ParameterofOneArray::calculateXs(vector <double> &conc1, vector <double> &conc2,char emodel)
{
	vector <double>::iterator pconc1 = conc1.begin();
	vector <double>::iterator pconc2 = conc2.begin();
	vector <double>::iterator pspoterror = SpotError.begin();

	vector <double> xs1;
	vector <double> xs2;

	double tmp1 = 0;
	double tmp2 = 0;

	while (pconc1!=conc1.end())
	{
		switch (emodel){
			case 'A':
				tmp1 = (*pconc1 * (MuSpot + *pspoterror))/(1/Ka + *pconc1 + *pconc2);
				tmp2 = (*pconc2 * (MuSpot + *pspoterror))/(1/Ka + *pconc1 + *pconc2);
				break;
			case 'M':
				tmp1 = (*pconc1 * (MuSpot * exp(*pspoterror)))/(1/Ka + *pconc1 + *pconc2);
				tmp2 = (*pconc2 * (MuSpot * exp(*pspoterror)))/(1/Ka + *pconc1 + *pconc2);
		}

		xs1.push_back(tmp1);
		xs2.push_back(tmp2);

		pconc1++;
		pconc2++;
		pspoterror++;
	}

	conc1 = xs1;    // need to return two values, so can't use "return" anymore.
	conc2 = xs2;
}

double ParameterofOneArray::costFunctionGivenXs(double x,struct my_f_params_AXs *paramsAXs)
{
	//double pka = (paramsAXs->ka);
	//double pmus = (paramsAXs->mus);
	double pinten = (paramsAXs->inten);
	double pxs = (paramsAXs->xs);
	double pp1 = (paramsAXs->p1);
	double pp2 = (paramsAXs->p2);
	double psigmam = (paramsAXs->sigmam);
	double psigmaa = (paramsAXs->sigmaa);

	double result = (((log(pinten - x) - log(pp1 * pxs + pp2))/(psigmam * SQRT2)) * ((log(pinten - x) - log(pp1 * pxs + pp2))/(psigmam * SQRT2))) + (x/(psigmaa * SQRT2)) * (x/(psigmaa * SQRT2));

	return result;
}

double ParameterofOneArray::calculateCostFunctionGivenXs(double ka,double mus,double inten,double xs,double p1,double p2,double sigmam,double sigmaa)
{
	struct my_f_params_AXs paramsAXs = {ka,mus,inten,xs,p1,p2,sigmam,sigmaa};
	
	/* calculate the boundaries for minimization */

	double lowerb = 0;
	double upperb = 0;

	if (inten > 0)
	{
		upperb = inten - p1 * xs - p2;

		if (lowerb > upperb)
		{
			double tmp = 0;

			tmp = lowerb;
			lowerb = upperb;
			upperb = tmp;
		}
	}else{
		upperb = inten;
		lowerb = inten - p1*xs - p2;
	}

	double iter = 0;
	double a = lowerb;
	double b = upperb;
	double x1 = a + (1 - GOLDEN_SECTION) * (b - a);
	double x2 = a + GOLDEN_SECTION * (b - a);
	double f1 = costFunctionGivenXs(x1,&paramsAXs);
	double f2 = costFunctionGivenXs(x2,&paramsAXs);

	while ((abs(b-a) > CONVERGENCE_ABS) && (iter < MAX_ITERATION))
	{
		if (f1 < f2)
		{
			b = x2;
			x2 = x1;
			f2 = f1;
			x1 = a + (1 - GOLDEN_SECTION) * (b - a);
			f1 = costFunctionGivenXs(x1,&paramsAXs);

		}else{
			if (f1 == f2)
			{
				a = x1;
				b = x2;
				x1 = a + (1 - GOLDEN_SECTION) * (b - a);
				x2 = a + GOLDEN_SECTION * (b - a);			
				f1 = costFunctionGivenXs(x1,&paramsAXs);
				f2 = costFunctionGivenXs(x2,&paramsAXs);

			}else{ 
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + GOLDEN_SECTION * (b - a);
				f2 = costFunctionGivenXs(x2,&paramsAXs);
			}
		}

		iter++;

	}

	double min = (a+b)/2;
	double fmin = costFunctionGivenXs(min,&paramsAXs);
	
	return 	fmin;	
}

double ParameterofOneArray::P1Function(double p1,struct my_f_params_P1 *paramsP1)
{
	double pka = (paramsP1->ka);
	double pmus = (paramsP1->mus);

	double pp2 = (paramsP1->p2);
	double psigmam = (paramsP1->sigmam);
	double psigmaa = (paramsP1->sigmaa);

	vector <double> inten = (paramsP1->inten);
	vector <double> xs = (paramsP1->xs);
	
	vector <double>::iterator ppinten = inten.begin();
	vector <double>::iterator ppxs = xs.begin();

	double costofone = 0;
	double costofall = 0;

	while (ppinten!=inten.end())  // sum of the cost function.
	{
		costofone = calculateCostFunctionGivenXs(pka,pmus,*ppinten,*ppxs,p1,pp2,psigmam,psigmaa);
		costofall += costofone;

		ppinten++;
		ppxs++;
	}

	return costofall;
}

double ParameterofOneArray::estimateP1(double ka,double mus,double p2,double sigmam,double sigmaa,vector <double> inten,vector<double> xs,double initialp1)
{
	struct my_f_params_P1 paramsP1 = {ka,mus,p2,sigmam,sigmaa,inten,xs};
		
	/* set the boundaries for minimization */ 

	double lowerb = 0;
	double upperb = 10 * P1Col1;

	double iter = 0;
	double a = lowerb;
	double b = upperb;
	double x1 = a + (1 - GOLDEN_SECTION) * (b - a);
	double x2 = a + GOLDEN_SECTION * (b - a);
	double f1 = P1Function(x1,&paramsP1);
	double f2 = P1Function(x2,&paramsP1);

	while ((abs(b-a) > 0.0001)) //&& (iter < 5))
	{
		if (f1 < f2)
		{
			b = x2;
			x2 = x1;
			f2 = f1;
			x1 = a + (1 - GOLDEN_SECTION) * (b - a);
			f1 = P1Function(x1,&paramsP1);

		}else{
			if (f1 == f2)
			{
				a = x1;
				b = x2;
				x1 = a + (1 - GOLDEN_SECTION) * (b - a);
				x2 = a + GOLDEN_SECTION * (b - a);			
				f1 = P1Function(x1,&paramsP1);
				f2 = P1Function(x2,&paramsP1);

			}else{
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + GOLDEN_SECTION * (b - a);
				f2 = P1Function(x2,&paramsP1);
			}
		}

		iter++;

	}

	double min = (a+b)/2;
		
	return 	min;		
}

/* =========================== public functions ======================== */

void ParameterofOneArray::setArrayID(int id)
{
	ArrayID = id;
}

void ParameterofOneArray::setMuSpot(double MuS)
{
	MuSpot = MuS;
}

void ParameterofOneArray::setKa(double ka)
{
	Ka = ka;
}

void ParameterofOneArray::setP1Col1(double p11)
{
	P1Col1 = p11;
}

void ParameterofOneArray::setP1Col2(double p12)
{
	P1Col2 = p12;
}

void ParameterofOneArray::setP2Col1(double p21)
{
	P2Col1 = p21;
}

void ParameterofOneArray::setP2Col2(double p22)
{
	P2Col2 = p22;
}

void ParameterofOneArray::setSigmaAddCol1(double SigmaA1)
{
	SigmaAddCol1 = SigmaA1;
}

void ParameterofOneArray::setSigmaAddCol2(double SigmaA2)
{
	SigmaAddCol2 = SigmaA2;
}

void ParameterofOneArray::setSigmaMulCol1(double SigmaM1)
{
	SigmaMulCol1 = SigmaM1;
}

void ParameterofOneArray::setSigmaMulCol2(double SigmaM2)
{
	SigmaMulCol2 = SigmaM2;
}

void ParameterofOneArray::setSigmaSpot(double s)
{
	SigmaSpot = s;
}

void ParameterofOneArray::setInitialP1(double maxInten)
{
	// estimate initial p1.

	double Max_Factor = 10;

	P1Col1 = (Max_Factor * maxInten - P2Col1)/MuSpot;
	P1Col2 = (Max_Factor * maxInten - P2Col2)/MuSpot;	
}

void ParameterofOneArray::setFineKa(vector <double> IntenCol1,vector <double> IntenCol2,vector <double> ConcCol1,vector <double> ConcCol2,char emodel)
{
	struct my_f_params_K paramsFK = {MuSpot,P1Col1,P2Col1,SigmaMulCol1,SigmaAddCol1,P1Col2,P2Col2,SigmaMulCol2,SigmaAddCol2,IntenCol1,IntenCol2,ConcCol1,ConcCol2};
		
	/* set the boundaries for minimization */ 

	double lowerb = log(0.0001/MuSpot);
	double upperb = 0;

	double iter = 0;
	double a = lowerb;
	double b = upperb;
	double x1 = a + (1 - GOLDEN_SECTION) * (b - a);
	double x2 = a + GOLDEN_SECTION * (b - a);
	double f1 = fineKaFunction(x1,emodel,&paramsFK);
	double f2 = fineKaFunction(x2,emodel,&paramsFK);

	while ((abs(b-a) > CONVERGENCE_ABS) && (iter < MAX_ITERATION))
	{
		if (f1 < f2)
		{
			b = x2;
			x2 = x1;
			f2 = f1;
			x1 = a + (1 - GOLDEN_SECTION) * (b - a);
			f1 = fineKaFunction(x1,emodel,&paramsFK);

		}else{
			if (f1 == f2)
			{
				a = x1;
				b = x2;
				x1 = a + (1 - GOLDEN_SECTION) * (b - a);
				x2 = a + GOLDEN_SECTION * (b - a);			
				f1 = fineKaFunction(x1,emodel,&paramsFK);
				f2 = fineKaFunction(x2,emodel,&paramsFK);

			}else{
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + GOLDEN_SECTION * (b - a);
				f2 = fineKaFunction(x2,emodel,&paramsFK);
			}
		}

		iter++;
	}

	double min = (a+b)/2;
	
	Ka = exp(min);  // fine ka.
	 
}

void ParameterofOneArray::calculateP1(vector <double> IntenCol1,vector <double> IntenCol2,vector <double> ConcCol1,vector <double> ConcCol2,char emodel)
{
	vector <double> xs1(ConcCol1);
	vector <double> xs2(ConcCol2);

	// calculate xs, the two parameters are changed after function calculateXs.

	calculateXs(xs1,xs2,emodel);

	double initialp11 = P1Col1;
	double initialp12 = P1Col2;

	P1Col1 = estimateP1(Ka,MuSpot,P2Col1,SigmaMulCol1,SigmaAddCol1,IntenCol1,xs1,initialp11);
	P1Col2 = estimateP1(Ka,MuSpot,P2Col2,SigmaMulCol2,SigmaAddCol2,IntenCol2,xs2,initialp12);
}

int ParameterofOneArray::getArrayID()
{
	return ArrayID;
}

double ParameterofOneArray::getKa()
{
	return Ka;
}

double ParameterofOneArray::getMuSpot() 
{
	return MuSpot;
}

double ParameterofOneArray::getP1Col1() 
{
	return P1Col1;
}

double ParameterofOneArray::getP1Col2() 
{
	return P1Col2;
}

double ParameterofOneArray::getP2Col1() 
{
	return P2Col1;
}

double ParameterofOneArray::getP2Col2()
{
	return P2Col2;
}

double ParameterofOneArray::getMuAddCol1() 
{
	return MuAddCol1;
}

double ParameterofOneArray::getMuAddCol2()
{
	return MuAddCol2;
}

double ParameterofOneArray::getSigmaAddCol1() 
{
	return SigmaAddCol1;
}

double ParameterofOneArray::getSigmaAddCol2()
{
	return SigmaAddCol2;
}

double ParameterofOneArray::getSigmaMulCol1() 
{
	return SigmaMulCol1;
}

double ParameterofOneArray::getSigmaMulCol2() 
{
	return SigmaMulCol2;
}

vector <double> ParameterofOneArray::getSpotError()
{
	return SpotError;
}

void ParameterofOneArray::calculateSigmaSpot()
{
	SigmaSpot = standarddeviation(SpotError);
}

double ParameterofOneArray::getSigmaSpot()
{
	return SigmaSpot;
}

