* BIOS 509 Poisson Regression Notes;
* Spring 2023;

* Read in the data;
proc import datafile = 'C:\Users\JRENGLE\OneDrive - Emory University\Documents\Teaching\BIOS 509 Spring 2023\Guest Lecture\mls.csv' out = mls DBMS = CSV replace;
run;

* Filter and transform some variables;
data mls;
	set mls;
	where Pos ^= 'GK' and Min > 500;
	GC = Gls + Ast;
	lWage = log(Wage);* / log(1.1); * Dividing by log(1.1) give the interpretation: for a 10% increase in wage, ...;
	logeMP = log(Min / 90);
run;

/*
By default, PROC GENOMD outputs the deviance and pearson X2 statistic. To conduct
goodness-of-fit tests you can manually compare to a chi-squared distribution with
the appropriate degrees of freedom to get a p-value.
You can also manually subtract these and compare to a chi-squared distribution 
with the appropriate df to conduct the likelihood ratio tests to compare nested
models.
*/

/**== Poisson Regression: Wage-only ==**/
proc genmod data = mls;
	model GC = lWage / dist = poisson link = log offset = logeMP;
	ods output ParameterEstimates = p1_parms ModelFit = p1_gof;
run;

* If you want p-values for the GOF tests, you need to calculate them manually;
data p1_gof_p;
	set p1_gof;
	pval = 1 - CDF('CHISQUARE', value, df);
run;
title 'GOF tests';
proc print data = p1_gof_p;
	where Criterion in ('Scaled Deviance','Scaled Pearson X2');
run;

* If you want rate ratios, you need to calculate them manually;
data p1_RR;
	set p1_parms;
	array ests(4) Estimate--UpperWaldCL;
	do i = 1 to 4;
		ests(i) = exp(ests(i));
	end;
run;
title 'RR Parameter Estimates';
proc print data = p1_RR;
	where Parameter ^= 'Scale';
	var Parameter Estimate LowerWaldCL UpperWaldCL;
run;

/**== Poisson Regression: wage and position ==**/
proc genmod data = mls;
	class Pos(ref = 'DF');
	model GC = lWage Pos / dist = poisson link = log offset = logeMP;
	ods output ParameterEstimates = p2_parms ModelFit = p2_gof;
run;

data p2_gof_p;
	set p2_gof;
	pval = 1 - CDF('CHISQUARE', value, df);
run;
title 'GOF tests';
proc print data = p2_gof_p;
	where Criterion in ('Scaled Deviance','Scaled Pearson X2');
run;

data p2_RR;
	set p2_parms;
	array ests(4) Estimate--UpperWaldCL;
	do i = 1 to 4;
		ests(i) = exp(ests(i));
	end;
run;
title 'RR Parameter Estimates';
proc print data = p2_RR;
	where Parameter ^= 'Scale' and ChiSq ^= .;
	var Parameter Level1 Estimate LowerWaldCL UpperWaldCL;
run;


/**== Poisson Regression: wage and position interaction ==**/
proc genmod data = mls;
	class Pos(ref = 'DF');
	model GC = lWage Pos lWage | Pos/ dist = poisson link = log offset = logeMP;
	ods output ParameterEstimates = p3_parms ModelFit = p3_gof;
run;

data p3_gof_p;
	set p3_gof;
	pval = 1 - CDF('CHISQUARE', value, df);
run;
title 'GOF tests';
proc print data = p3_gof_p;
	where Criterion in ('Scaled Deviance','Scaled Pearson X2');
run;

data p3_RR;
	set p3_parms;
	array ests(4) Estimate--UpperWaldCL;
	do i = 1 to 4;
		ests(i) = exp(ests(i));
	end;
run;
title 'RR Parameter Estimates';
proc print data = p3_RR;
	where Parameter ^= 'Scale' and ChiSq ^= .;
	var Parameter Level1 Estimate LowerWaldCL UpperWaldCL;
run;

* Table of Deviances for all 3 Poisson model fits;
data pdev;
	set p1_gof p2_gof p3_gof;
	where Criterion = 'Scaled Deviance';
run;
proc print;
run;
* You can manually subtract these deviances to conduct drop-in-deviance (LR) tests;


/**== Quasi-Poisson Regression: main effects model ==**/
proc genmod data = mls;
	class Pos(ref = 'DF');
	model GC = lWage Pos / dist = poisson link = log offset = logeMP scale = Pearson;
	ods output ParameterEstimates = qp2_parms ModelFit = qp2_gof;
run;
/*
Note the scale parameter reported by SAS is the square root of the dispersion parameter.
This is the value that the Poisson reg. SEs are multiplied by to get the new SEs.
There are GOF statistics output, but these may not be reliable since they assume the
scale parameter is known;
*/

data qp2_RR;
	set qp2_parms;
	array ests(4) Estimate--UpperWaldCL;
	do i = 1 to 4;
		ests(i) = exp(ests(i));
	end;
run;

title 'RR Parameter Estimates';
proc print data = qp2_RR;
	where Parameter ^= 'Scale' and ChiSq ^= .;
	var Parameter Level1 Estimate LowerWaldCL UpperWaldCL;
run;

/**== Quasi-Poisson Regression: interaction model ==**/
proc genmod data = mls;
	class Pos(ref = 'DF');
	model GC = lWage Pos lWage | Pos/ dist = poisson link = log offset = logeMP scale = Pearson;
	ods output ParameterEstimates = qp3_parms ModelFit = qp3_gof;
run;

data qp3_RR;
	set qp3_parms;
	array ests(4) Estimate--UpperWaldCL;
	do i = 1 to 4;
		ests(i) = exp(ests(i));
	end;
run;

title 'RR Parameter Estimates';
proc print data = qp3_RR;
	where Parameter ^= 'Scale' and ChiSq ^= .;
	var Parameter Level1 Estimate LowerWaldCL UpperWaldCL;
run;

/*
We can use an F-test to compare these "nested" quasi-poisson models. I'm not
sure how to do this automatically in SAS, but here is how to do it manually.
*/

proc sql noprint;
	select Value into :qp2_dev from qp2_gof where Criterion = 'Deviance';
	select Value into :qp3_dev from qp3_gof where Criterion = 'Deviance';
	select Estimate**2 into :qp3_phi from qp3_parms where Parameter = 'Scale';
	select DF into :qp3_df from qp3_gof where Criterion = 'Deviance';
quit;
data qp_ftest;
	F = (&qp2_dev. - &qp3_dev.) / 2 / &qp3_phi.;
	pval = 1 - CDF('F', F, 2, &qp3_df.);
run;
proc print;
run;


/**== Negative Binomial Regression: Main effects model ==**/
proc genmod data = mls;
	class Pos(ref = 'DF');
	model GC = lWage Pos / dist = NB link = log offset = logeMP;
	ods output ParameterEstimates = nb2_parms ModelFit = nb2_gof;
run;

data nb2_RR;
	set nb2_parms;
	array ests(4) Estimate--UpperWaldCL;
	do i = 1 to 4;
		ests(i) = exp(ests(i));
	end;
run;
title 'RR Parameter Estimates';
proc print data = nb2_RR;
	where Parameter ^= 'Scale' and ChiSq ^= .;
	var Parameter Level1 Estimate LowerWaldCL UpperWaldCL;
run;

/**== Negative Binomial Regression: interaction model ==**/
proc genmod data = mls;
	class Pos(ref = 'DF');
	model GC = lWage Pos lWage | Pos / dist = NB link = log offset = logeMP;
	ods output ParameterEstimates = nb3_parms ModelFit = nb3_gof;
run;

data nb3_RR;
	set nb3_parms;
	array ests(4) Estimate--UpperWaldCL;
	do i = 1 to 4;
		ests(i) = exp(ests(i));
	end;
run;
title 'RR Parameter Estimates';
proc print data = nb3_RR;
	where Parameter ^= 'Scale' and ChiSq ^= .;
	var Parameter Level1 Estimate LowerWaldCL UpperWaldCL;
run;

/*
There is no formal GOF test for NB. But we can still compare these two
models using a LRT
*/

proc sql noprint;
	select Value into :nb2_ll from nb2_gof where Criterion = 'Full Log Likelihood';
	select Value into :nb3_ll from nb3_gof where Criterion = 'Full Log Likelihood';
quit;
data nb_lrt;
	LRT = -2*(&nb2_ll. - &nb3_ll.);
	pval = 1 - CDF('Chisquare', LRT, 2);
run;
proc print;
run;


/*
Overdispersion test via negative binomial.
*/

proc sql noprint;
	select Value into :p3_ll from p3_gof where Criterion = 'Full Log Likelihood';
quit;
data nb_lrt;
	LRT = -2*(&p3_ll. - &nb3_ll.);
	pval = 1 - CDF('Chisquare', LRT, 1);
run;
proc print;
run;


/**== Residual diagnostics example for the NB model ==**/
proc genmod data = mls plots = (STDRESCHI(index xbeta) STDRESDEV(index xbeta)); * These plots help assess linearity and independence;
	class Pos(ref = 'DF');
	model GC = lWage Pos lWage | Pos / dist = NB link = log offset = logeMP;
	output out = nb_res STDRESCHI = stdresp STDRESDEV = stdresd;
run;

* We can use a qq-plot or something similar to test normality;
proc univariate data = nb_res;
	qqplot stdresd / normal(mu = 0 sigma = 1) square;
run;
