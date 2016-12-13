* Joint Model A;
* Binary covariate;

libname data65 "c:\liulei\paper\paper11\data65";

%let datadir=c:\liulei\paper\paper11\data65;
data data65.est1_all;
set _null_;
run;
data data65.fit1_all;
set _null_;
run;

data data65.theta1_all;
set _null_;
run;
data data65.est2_all;
set _null_;
run;
data data65.fit2_all;
set _null_;
run;
data data65.theta2_all;
set _null_;
run;

data data65.est3_all;
set _null_;
run;
data data65.fit3_all;
set _null_;
run;
data data65.theta3_all;
set _null_;
run;


%macro cure();

%do ii=1 %to 200;

title "Replicate &ii";

PROC IMPORT OUT= WORK.one 
            DATAFILE= "&datadir\&ii..csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=20; 
RUN;

data one;
set one;
aa=1;
run;

data two;
set one;
if event=1;
run;

proc univariate data=two noprint;
var stoptime; 
output out=quant_r pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qr; 
run;

data quant_r;
set quant_r;
aa=1;
run;

data three;
set one;
if event=2;
run;

proc univariate data=three noprint;
var stoptime; 
output out=quant_d pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qd; 
run;

data quant_d;
set quant_d;
aa=1;
run;

* Merge data with the quantiles;

data four;
merge one quant_r quant_d;
by aa;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data five;
set four;
array quant_r {6} qr0  qr20  qr40  qr60  qr80  qr100;
array quant_d {6} qd0  qd20  qd40  qd60  qd80  qd100;

array dur_r {5} dur_r1-dur_r5;
array dur_d {5} dur_d1-dur_d5;

array event_r {5} event_r1-event_r5;
array event_d {5} event_d1-event_d5;

do i=1 to 5;
	dur_r{i}=0;
	dur_d{i}=0;
	event_r{i}=0;
	event_d{i}=0;
end;

* For recurrent event;
if event=1 then do;
	do i=2 to 6;
		if stoptime<=quant_r{i} then do;
			event_r{i-1}=1;
			i=6;
		end;
	end;
end;

else do; /* If death or censored observation */
	do i=2 to 6;
		if stoptime<=quant_r{i} then do;
			dur_r{i-1}=stoptime-quant_r{i-1};
			i=6;
		end;
		else do;
			dur_r{i-1}=quant_r{i}-quant_r{i-1};
		end;
	end;

	do i=2 to 6;
		if stoptime<=quant_d{i} then do;
			event_d{i-1}=(event=2);
			dur_d{i-1}=stoptime-quant_d{i-1};
			i=6;
		end;
		else do;
			dur_d{i-1}=quant_d{i}-quant_d{i-1};
		end;
	end;
end;

run;

proc sort data=five;
by id stoptime;
run;

data six;
set five;
by id;
retain enum;
if first.id then enum=0;
enum=enum+1;
run;

/*
proc freq data=six;
table enum;
run;

data seven;
set six;
by id;
if last.id;
run;
proc freq data=seven;
table enum;
run;
*/


* With cure fraction;


proc nlmixed data=six qpoints=5;

parms shape1=1.25 shape2=1.25 scale1=4 scale2=10
alpha0=-.5 alpha1=1 beta1=1 delta1=1 vara=1 gamma=1;
bounds shape1 shape2 scale1 scale2  >=0;

eta= alpha0 + alpha1 * covar ; 		/* for probability of no recurrent event */
prob=1/(1+exp(-eta));

base_haz_r=shape1* stoptime **(shape1-1)/(scale1 ** shape1)  ;
cum_base_haz_r=stoptime **shape1/(scale1 ** shape1)  ;
base_haz_d=shape2* stoptime **(shape2-1)/(scale2 ** shape2)  ;
cum_base_haz_d=stoptime **shape2/(scale2 ** shape2)  ;


*if stoptime >=qd100 then cum_base_haz_d=cum_base_haz_d + .02;

mu1= beta1 * covar +  a;			/* for recurrent event */
mu2= delta1 * covar + gamma * a;	/* for death event */

loglik1=-exp(mu1) * cum_base_haz_r;
loglik2=-exp(mu2) * cum_base_haz_d;


if event=1 then loglik = log(base_haz_r) + mu1; 		/*log likelihood for recurrent event */
if event=0 then do;
	if enum>1 then loglik= log(1-prob) + loglik1 + loglik2; 				/* log likelihood for subjects with at least 1 recurrent event */
	if enum=1 then loglik=log(prob + (1-prob) * exp(loglik1 + loglik2));					/*log likelihood for subjects without recurrent event */
end;
if event=2 then do;
	loglik= log(1-prob) + loglik1 + log(base_haz_d) + mu2 + loglik2; 				/* log likelihood for non-cured patients who died */
end;


model id ~ general(loglik);

random a ~ normal(0, vara) subject=id;

ods output ParameterEstimates=est2 FitStatistics=fit2 ; 

run;

* with cure fraction, set gamma=1;

proc nlmixed data=six qpoints=5;

parms shape1=1.25 shape2=1.25 scale1=4 scale2=10
alpha0=-.5 alpha1=1 beta1=1 delta1=1 vara=1 ;
bounds shape1 shape2 scale1 scale2  >=0;

eta= alpha0 + alpha1 * covar ; 		/* for probability of no recurrent event */
prob=1/(1+exp(-eta));

base_haz_r=shape1* stoptime **(shape1-1)/(scale1 ** shape1)  ;
cum_base_haz_r=stoptime **shape1/(scale1 ** shape1)  ;
base_haz_d=shape2* stoptime **(shape2-1)/(scale2 ** shape2)  ;
cum_base_haz_d=stoptime **shape2/(scale2 ** shape2)  ;

*if stoptime >=qd100 then cum_base_haz_d=cum_base_haz_d + .02;

mu1= beta1 * covar +  a;			/* for recurrent event */
mu2= delta1 * covar + a;	/* for death event */

loglik1=-exp(mu1) * cum_base_haz_r;
loglik2=-exp(mu2) * cum_base_haz_d;

if event=1 then loglik = log(base_haz_r) + mu1; 		/*log likelihood for recurrent event */
if event=0 then do;
	if enum>1 then loglik= log(1-prob) + loglik1 + loglik2; 				/* log likelihood for subjects with at least 1 recurrent event */
	if enum=1 then loglik=log(prob + (1-prob) * exp(loglik1 + loglik2));					/*log likelihood for subjects without recurrent event */
end;
if event=2 then do;
	loglik= log(1-prob) + loglik1 + log(base_haz_d) + mu2 + loglik2; 				/* log likelihood for non-cured patients who died */
end;


model id ~ general(loglik);

random a ~ normal(0, vara) subject=id;


ods output ParameterEstimates=est3 FitStatistics=fit3; 

run;



data data65.est2_all;
set data65.est2_all est2;
run;

data data65.fit2_all;
set data65.fit2_all fit2;
run;


data data65.est3_all;
set data65.est3_all est3;
run;

data data65.fit3_all;
set data65.fit3_all fit3;
run;



%end;

%mend;

%cure();
data result;
set _null_;
run;
%macro organize_result(varname);
data two;
set data65.est2_all;
if parameter="&varname";
&varname=estimate;
se_&varname=standarderror;
keep &varname se_&varname;
run;

data result;
merge result two;
run;
%mend;
%organize_result(beta1);
%organize_result(delta1);
%organize_result(alpha0);
%organize_result(alpha1);
%organize_result(vara);
%organize_result(gamma);

%organize_result(shape1);
%organize_result(scale1);
%organize_result(shape2);
%organize_result(scale2);

data _null_;
set result;
file "c:\liulei\paper\paper11\data65\AGQresult1.txt";
put alpha0 alpha1 beta1 delta1 gamma  vara shape1 scale1 shape2 scale2 se_alpha0 se_alpha1 se_beta1 se_delta1 se_gamma se_vara se_shape1 se_scale1 se_shape2 se_scale2; 
run;
data result;
set _null_;
run;
%macro organize_result(varname);
data two;
set data65.est3_all;
if parameter="&varname";
&varname=estimate;
se_&varname=standarderror;
keep &varname se_&varname;
run;

data result;
merge result two;
run;
%mend;
%organize_result(beta1);
%organize_result(delta1);
%organize_result(alpha0);
%organize_result(alpha1);
%organize_result(vara);

%organize_result(shape1);
%organize_result(scale1);
%organize_result(shape2);
%organize_result(scale2);

data _null_;
set result;
file "c:\liulei\paper\paper11\data65\AGQresult2.txt";
put alpha0 alpha1 beta1 delta1 vara  shape1 scale1 shape2 scale2  se_alpha0 se_alpha1 se_beta1 se_delta1 se_vara se_shape1 se_scale1 se_shape2 se_scale2; 
run;

**** Summarizing fit statistics;

data result2;
set _null_;
run;
%macro organize_result(criterion, name);
data two;
set data65.fit2_all;
if descr="&criterion";
&name=value;
keep &name;
run;

data result2;
merge result2 two;
run;
%mend;
%organize_result(-2 Log Likelihood, loglik2);
%organize_result(AIC (smaller is better), aic2);
%organize_result(AICC (smaller is better), aicc2);
%organize_result(BIC (smaller is better), bic2);

data result3;
set _null_;
run;
%macro organize_result(criterion, name);
data two;
set data65.fit3_all;
if descr="&criterion";
&name=value;
keep &name;
run;

data result3;
merge result3 two;
run;
%mend;
%organize_result(-2 Log Likelihood, loglik3);
%organize_result(AIC (smaller is better), aic3);
%organize_result(AICC (smaller is better), aicc3);
%organize_result(BIC (smaller is better), bic3);

data result;
merge result2 result3;
run;


data _null_;
set result;
file "c:\liulei\paper\paper11\data65\fitstat.txt";
put loglik2 loglik3  aic2 aic3  aicc2 aicc3  bic2 bic3 ;
run;
