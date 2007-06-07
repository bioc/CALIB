#include "SpikeofOneArray.h"
#include "NormalizedData.h"

#include <cmath>
#include <algorithm>
#include <iterator>

using namespace std;

//#ifdef __cplusplus
extern "C" {
//#endif
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>


SEXP estimation(SEXP int1,SEXP int2,SEXP conc1,SEXP conc2,SEXP type,SEXP maxinten,SEXP errormodel,SEXP mus)
{
	double *xint1,*xint2,*xconc1,*xconc2,*xmax,*xgivenmus;
	const char *xerrormodel;
	
	xint1 = REAL(int1);
	xint2 = REAL(int2);
	xconc1 = REAL(conc1);
	xconc2 = REAL(conc2);
	xmax = REAL(maxinten);
	xerrormodel = CHAR(STRING_ELT(errormodel,0));
	xgivenmus = REAL(mus);

	int len = 0;
	len = length(int1);

	vector <double> inta(xint1,xint1 + len);
	vector <double> intb(xint2,xint2 + len);
	vector <double> conca(xconc1,xconc1 + len);
	vector <double> concb(xconc2,xconc2 + len);
	vector <string> spottype;
	SpikeofOneArray spikeofone;

	int i = 0;
	char *tmp; 

	// spot type: string vector
	for (i=0;i<len;i++)
	{	
		tmp = R_alloc(strlen(CHAR(STRING_ELT(type,i))),sizeof(char));
		strcpy(tmp,CHAR(STRING_ELT(type,i)));
		spottype.push_back(tmp);
	}

	// set all the values.
	
	spikeofone.setSpikes(conca,concb,inta,intb,spottype); 
	
	//define the variables used in estimating procedure.
	ParameterofOneArray parameterofone;
	Spike negative;
	vector <double> negcol1inten;
	vector <double> negcol2inten;
	double addvar1 = 0;
	double addvar2 = 0;
	double addmean1 = 0;
	double addmean2 = 0;
	Spike calibration;
	vector <double> calcol1inten;
	vector <double> calcol2inten;
	vector <double> calcol1conc;
	vector <double> calcol2conc;
	double mulvar1 = 0;
	double mulvar2 = 0;
	double mulvar = 0;
	Spike dataset;
	vector <double> d1inten;
	vector <double> d2inten;
	vector <double> d1conc;
	vector <double> d2conc;
	int sid = 0;

	// estimation begins:
	// select negative spikes of one array
	negative = spikeofone.selectSpikes("Negative");
		
	//get intensities of each color of negative spikes
	negcol1inten = spikeofone.getCol1Inten(negative);
	negcol2inten = spikeofone.getCol2Inten(negative);

	//calculate additive error for each color by using negative spikes
	addvar1 = spikeofone.calculateAdditiveVariance(negcol1inten);
	addvar2 = spikeofone.calculateAdditiveVariance(negcol2inten);

	// calculate p2 = MuAddtive 
	addmean1 = spikeofone.calculateAdditiveMean(negcol1inten); 
	addmean2 = spikeofone.calculateAdditiveMean(negcol2inten); 

	//select calibration spikes of one array
	calibration = spikeofone.selectSpikes("Calibration");

	//get intensities of each color of calibration spikes
	calcol1inten = spikeofone.getCol1Inten(calibration);
	calcol2inten = spikeofone.getCol2Inten(calibration);
	calcol1conc = spikeofone.getCol1Conc(calibration);
	calcol2conc = spikeofone.getCol2Conc(calibration);
	
	//calculate multiplicative error for each color by using calibration spikes
	mulvar = spikeofone.calculateMultiplicativeVariance(calcol1inten,calcol2inten);
	mulvar1 = 1/1.414214 * mulvar;
	mulvar2 = mulvar1;

	//====================== get the data set we use to estimate parameters ===================

	// data set: all spikes except for negative controls
	dataset = spikeofone.kickoutNegatives("Negative");

	// get the intensities and concentrations of each color 
	d1inten = spikeofone.getCol1Inten(dataset);
	d2inten = spikeofone.getCol2Inten(dataset);
	d1conc = spikeofone.getCol1Conc(dataset);
	d2conc = spikeofone.getCol2Conc(dataset);
	sid = spikeofone.getArrayID(); 

	//================================= parameter estimation ===================================

	parameterofone.setArrayID(sid); 
	parameterofone.setMuSpot(*xgivenmus);
	parameterofone.setP2Col1(addmean1);
	parameterofone.setP2Col2(addmean2);
	parameterofone.setSigmaAddCol1(addvar1);
	parameterofone.setSigmaAddCol2(addvar2);
	parameterofone.setSigmaMulCol1(mulvar1);
	parameterofone.setSigmaMulCol2(mulvar2);
	parameterofone.setInitialP1(*xmax);
	//parameterofone.setCoarseKa(calcol1inten,calcol2inten,calcol1conc,calcol2conc);
	parameterofone.setFineKa(d1inten,d2inten,d1conc,d2conc,*xerrormodel);
	parameterofone.calculateP1(d1inten,d2inten,d1conc,d2conc,*xerrormodel);
	parameterofone.calculateSigmaSpot(); 

	// estimation ends.

	// store all the parameters in a list and pass back to R.

	SEXP parlist, parlistname,parmus,parka,parp1,parp2,parsigmaa,parsigmam,parsigmas,parspoterror;
	char *name[8] = {"MuS","Ka","P1","P2","SigmaA","SigmaM","SigmaS","SpotError"};
	double *xmus,*xka,*xp1,*xp2,*xsigmaa,*xsigmam,*xsigmas,*xspoterror;
	

	// element "MuS"(numeric) in the list:
	PROTECT(parmus = allocVector(REALSXP,1));
	xmus = REAL(parmus);
	xmus[0] = parameterofone.getMuSpot();

	// element "Ka"(numeric) in the list:
	PROTECT(parka = allocVector(REALSXP,1));
	xka = REAL(parka);
	xka[0] = parameterofone.getKa();

	// element "P1"(two numerics) in the list:
	PROTECT(parp1 = allocVector(REALSXP,2));
	xp1 = REAL(parp1);
	xp1[0] = parameterofone.getP1Col1();
	xp1[1] = parameterofone.getP1Col2();

	// element "P2" (two numerics) in the list:

	PROTECT(parp2 = allocVector(REALSXP,2));
	xp2 = REAL(parp2);
	xp2[0] = parameterofone.getP2Col1();
	xp2[1] = parameterofone.getP2Col2();

	// element "SigmaA" (two numerics) in the list:
	PROTECT(parsigmaa = allocVector(REALSXP,2));
	xsigmaa = REAL(parsigmaa);
	xsigmaa[0] = parameterofone.getSigmaAddCol1(); 
	xsigmaa[1] = parameterofone.getSigmaAddCol2();

	// element "SigmaM" (two numerics) in the list:
	PROTECT(parsigmam = allocVector(REALSXP,2));
	xsigmam = REAL(parsigmam);
	xsigmam[0] = parameterofone.getSigmaMulCol1();
	xsigmam[1] = parameterofone.getSigmaMulCol2();

	// element "SigmaS" (two numerics) in the list:
	PROTECT(parsigmas = allocVector(REALSXP,1));
	xsigmas = REAL(parsigmas);
	xsigmas[0] = parameterofone.getSigmaSpot(); 

	// element "SpotError" (n numerics) in the list:
	vector <double> error;
	R_len_t size = 0;
	error = parameterofone.getSpotError();
	size = error.size();
	PROTECT(parspoterror = allocVector(REALSXP,size));
	xspoterror = REAL(parspoterror);
	for (i=0;i<size;i++)
	{
		xspoterror[i] = error[i];
	}
	
	// names of each list element:
	PROTECT(parlistname = allocVector(STRSXP,8));
	for(i=0;i<8;i++)
	{
		SET_STRING_ELT(parlistname,i,mkChar(name[i])); 
	}

	// set each element into the list:
	PROTECT(parlist = allocVector(VECSXP, 8)); 
    SET_VECTOR_ELT(parlist, 0, parmus);
	SET_VECTOR_ELT(parlist, 1, parka);
	SET_VECTOR_ELT(parlist, 2, parp1);
	SET_VECTOR_ELT(parlist, 3, parp2);
	SET_VECTOR_ELT(parlist, 4, parsigmaa);
	SET_VECTOR_ELT(parlist, 5, parsigmam);
	SET_VECTOR_ELT(parlist, 6, parsigmas);
	SET_VECTOR_ELT(parlist, 7, parspoterror);	
	setAttrib(parlist, R_NamesSymbol, parlistname); 

	UNPROTECT(10);
	return(parlist);	
}

vector <ParameterofOneArray> readinparlist(SEXP parlist, int ncol)
{
	int len_par = length(parlist);

	SEXP names = getAttrib(parlist,R_NamesSymbol);

	SEXP parmus = NULL,parka = NULL,parp1 = NULL,parp2 = NULL,parsigmaa = NULL,parsigmam = NULL,parsigmas = NULL;

	int i;

	for (i=0;i<len_par;i++)
	{
		if (strcmp(CHAR(STRING_ELT(names,i)),"MuS") == 0)
		{
			parmus = VECTOR_ELT(parlist,i);		
		}

		if (strcmp(CHAR(STRING_ELT(names,i)),"Ka") == 0)
		{
			parka = VECTOR_ELT(parlist,i);

		}

		if (strcmp(CHAR(STRING_ELT(names,i)),"P1") == 0)
		{
			parp1 = VECTOR_ELT(parlist,i);

		}

		if (strcmp(CHAR(STRING_ELT(names,i)),"P2") == 0)
		{
			parp2 = VECTOR_ELT(parlist,i);

		}

		if (strcmp(CHAR(STRING_ELT(names,i)),"SigmaA") == 0)
		{
			parsigmaa = VECTOR_ELT(parlist,i);

		}

		if (strcmp(CHAR(STRING_ELT(names,i)),"SigmaM") == 0)
		{
			parsigmam = VECTOR_ELT(parlist,i);

		}

		if (strcmp(CHAR(STRING_ELT(names,i)),"SigmaS") == 0)
		{
			parsigmas = VECTOR_ELT(parlist,i);

		}
	}

	double *xmus,*xka,*xp1,*xp2,*xsigmaa,*xsigmam,*xsigmas;

	xmus = REAL(parmus);
	xka = REAL(parka);
	xp1 = REAL(parp1);
	xp2 = REAL(parp2);
	xsigmaa = REAL(parsigmaa);
	xsigmam = REAL(parsigmam);
	xsigmas = REAL(parsigmas);

	vector <ParameterofOneArray> parameters;
	ParameterofOneArray parameterofone;

	for (i=0;i<ncol;i++)
	{
		parameterofone.setArrayID(i+1);
		parameterofone.setMuSpot(xmus[i]);
		parameterofone.setKa(xka[i]);
		parameterofone.setP1Col1(xp1[2*i]);
		parameterofone.setP1Col2(xp1[2*i+1]);
		parameterofone.setP2Col1(xp2[2*i]);
		parameterofone.setP2Col2(xp2[2*i+1]);
		parameterofone.setSigmaAddCol1(xsigmaa[2*i]);
		parameterofone.setSigmaAddCol2(xsigmaa[2*i+1]);
		parameterofone.setSigmaMulCol1(xsigmam[2*i]);
		parameterofone.setSigmaMulCol2(xsigmam[2*i+1]);
		parameterofone.setSigmaSpot(xsigmas[i]);

		parameters.push_back(parameterofone);
	}

	return (parameters);

}

SEXP adjustment(SEXP int1,SEXP int2,SEXP parlist,SEXP whichcolor)
{
	//read intensities.

	double *xint1,*xint2;
	int *xwhich;
	
	xint1 = REAL(int1);
	xint2 = REAL(int2);
	xwhich = INTEGER(whichcolor);

	int len_inten = length(int1);

	vector <double> inten1(xint1,xint1 + len_inten);
	vector <double> inten2(xint2,xint2 + len_inten);

	RawDataofOneArray rawofone;

	rawofone.setArrayID(1);
	rawofone.setRawData(inten1,inten2);

	// read parameters.

	vector <ParameterofOneArray> parameters;
	ParameterofOneArray parameterofone;

	parameters = readinparlist(parlist,1);
	parameterofone = parameters[0];

	// adjustment.

	double newp2 = 0;

	if (*xwhich == 1)
	{
		newp2 = rawofone.adjustCy3(parameterofone); 
	}else
	{
		if (*xwhich == 2)
		{
			newp2 = rawofone.adjustCy5(parameterofone); 
		}
	}
	
    SEXP result;

	PROTECT(result = allocVector(REALSXP, 1));

	REAL(result)[0] = newp2;

	UNPROTECT(1);

	return(result);
	
}	

SEXP normalization(SEXP int1,SEXP int2,SEXP dim,SEXP id,SEXP designa,SEXP designc,SEXP designd,SEXP parlist,SEXP errormodel)
{
	// read intensities of all arrays. matrix, so need extra parameter "dim".

	double *xint1;
	double *xint2;
	int *xdim;
	const char *xerrormodel;

	xint1 = REAL(int1);
	xint2 = REAL(int2);
	xdim = INTEGER(dim);
	xerrormodel = CHAR(STRING_ELT(errormodel,0));

	int nrow = xdim[0]; // number of clone
	int ncol = xdim[1]; // number of array

	int i = 0; // i and j are index.
	int j = 0;

	// cloneid : string vector

	int len_cloneid = nrow;

	vector <string> cloneid; 
	char* tmp;  // temperary variable for clone id vector.

	for (i=0;i<len_cloneid;i++)
	{	
		tmp = R_alloc(strlen(CHAR(STRING_ELT(id,i))),sizeof(char));
		strcpy(tmp,CHAR(STRING_ELT(id,i)));
		cloneid.push_back(tmp);
	}

    // read intensities.

	vector <RawDataofOneArray> rawdata;
	RawDataofOneArray rawofone;

	for (i=0;i<ncol;i++)
	{
		vector <double> inten1(xint1,xint1 + nrow);
		vector <double> inten2(xint2,xint2 + nrow);

		rawofone.setArrayID(i+1);
		rawofone.setCloneID(cloneid);
		rawofone.setRawData(inten1,inten2);

		rawdata.push_back(rawofone); 

		xint1 += nrow;
		xint2 += nrow;
	}

	// read parameters. parameters are given by a list.
	
	vector <ParameterofOneArray> parameters;

	parameters = readinparlist(parlist,ncol);

	// read design matrix.

	int *xdesigna;
	int *xdesignc;
	int *xdesignd;

	xdesigna = INTEGER(designa);
    xdesignc = INTEGER(designc);
	xdesignd = INTEGER(designd);

	int len_design = length(designa);

	vector <int> a(xdesigna,xdesigna + len_design);
	vector <int> c(xdesignc,xdesignc + len_design);
	vector <int> d(xdesignd,xdesignd + len_design);

	Design des;
	des.setArray(a);
	des.setCond(c);
	des.setDye(d);

	// find out design block from design matrix.

	vector <Design> des_block;
	des_block = des.splitBlock(des); 

	// normalize data block by block.

	vector <Design>::iterator pdes_block = des_block.begin();

	Design block;
	NormalizedData ndata;
	
	vector <RawDataofOneArray> rawdataofoneblock;
	vector <NDataofOneGene> normalizeddata;
	vector <int> condition;
	
	while (pdes_block!=des_block.end())
	{
		block = *pdes_block;

		rawdataofoneblock = rawofone.selectBlockData(rawdata,block); 

		vector <RawDataofOneArray>::iterator prawofblock = rawdataofoneblock.begin(); 
		size_t size = rawdataofoneblock.size();

		ndata.setCloneID(*prawofblock);
		ndata.normalizeAllClone(prawofblock,size,block,parameters,*xerrormodel);

		vector <NDataofOneGene> tmpndata;
		vector <int> tmpcond;

		tmpndata = ndata.getData(); 
		tmpcond = ndata.getCondition();

		copy(tmpndata.begin(),tmpndata.end(),back_inserter(normalizeddata));
		copy(tmpcond.begin(),tmpcond.end(),back_inserter(condition));

		pdes_block++;
	}
	
	// output the result:

	int rnrow = 0;
	int rncol = 0;

	vector <string> result_cloneid;

	result_cloneid = ndata.getCloneID();

	rnrow = normalizeddata.size();
	rncol = condition.size();

	NDataofOneGene dataofoneclone;

	SEXP result,dims,dimnames;

	// the result
	PROTECT(result = allocVector(REALSXP, rnrow * rncol));

	for (i=0;i<rnrow;i++)
	{
		dataofoneclone = normalizeddata[i];
		
		for (j=0;j<rncol;j++)
		{
			REAL(result)[i+rnrow*j] = dataofoneclone.Conc[j];
		}
	}

	// dimension. rnrow(row number of result.) and rncol(column number of result.)

	PROTECT(dims = allocVector(INTSXP,2));
	INTEGER(dims)[0] = rnrow;
	INTEGER(dims)[1] = rncol;
    setAttrib(result,R_DimSymbol,dims);

	// dimension name.
	SEXP dimnames_row;
	SEXP dimnames_col;
	
	vector <string>::iterator pcloneid = result_cloneid.begin();
	string tmp_id;

	PROTECT(dimnames_row = allocVector(STRSXP,rnrow));
	for (i=0;i<rnrow;i++)
	{
		tmp_id = *pcloneid;
		char *name = strdup(tmp_id.c_str());;
		SET_STRING_ELT(dimnames_row,i,mkChar(name));
		
		pcloneid++;		
	}

	int *xdimnamescol;

	PROTECT(dimnames_col = allocVector(INTSXP,rncol));
	xdimnamescol = INTEGER(dimnames_col);
	for(i=0;i<rncol;i++)
	{
		xdimnamescol[i] = condition[i];
	}

	PROTECT(dimnames = allocVector(VECSXP,2));
	SET_VECTOR_ELT(dimnames,0,dimnames_row);
	SET_VECTOR_ELT(dimnames,1,dimnames_col);
	setAttrib(result,R_DimNamesSymbol,dimnames);

	UNPROTECT(5);
	return(result); 
}

R_CallMethodDef callMethods[] = {
	{"estimation",(DL_FUNC)&estimation,8},
	{"adjustment",(DL_FUNC)&adjustment,4},
	{"normalization",(DL_FUNC)&normalization,9},
	{NULL,NULL,0}
};

void R_init_CALIB(DllInfo *info)
{
	R_registerRoutines(info,NULL,callMethods,NULL,NULL);
}

#ifdef __cplusplus
} // extern "C"
#endif
