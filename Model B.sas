libname data51 "c:\liulei\paper\";

%let datadir=c:\liulei\paper\;


PROC IMPORT OUT= WORK.one 
            DATAFILE= "&datadir\1.csv" 
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

* No cure fraction;
proc nlmixed data=six qpoints=5 MAXITER=5000;

parms r01=0.25 r02=0.25 r03=0.24 r04=0.25 r05=0.25
	 	h01=0.1 h02=0.1 h03=0.1 h04=0.1 h05=0.1
		beta1=1 delta1=1 vara=.4;
bounds r01 r02 r03 r04 r05 h01 h02 h03 h04 h05 vara >=0;

base_haz_r=r01 * event_r1 + r02 * event_r2 + r03 * event_r3 + r04 * event_r4 + r05 * event_r5 ;
cum_base_haz_r=r01 * dur_r1 + r02 * dur_r2 + r03 * dur_r3 + r04 * dur_r4 + r05 * dur_r5 ;

base_haz_d=h01 * event_d1 + h02 * event_d2 + h03 * event_d3 + h04 * event_d4 + h05 * event_d5 ;
cum_base_haz_d=h01 * dur_d1 + h02 * dur_d2 + h03 * dur_d3 + h04 * dur_d4 + h05 * dur_d5 ;

mu1= beta1 * covar +  a;			/* for recurrent event */

mu2= delta1 * covar + a;	/* for death event */

loglik1=-exp(mu1) * cum_base_haz_r;

loglik2=-exp(mu2) * cum_base_haz_d;
if event=1 then loglik= log(base_haz_r) + mu1 ; 					/*log likelihood for recurrent event */
if event=2 then loglik=loglik1 + log(base_haz_d) + mu2 + loglik2;	/*log likelihood for death */
if event=0 then loglik=loglik1 + loglik2;							/*log likelihood for censoring */

model stoptime ~ general(loglik);
random a ~ normal(0, vara) subject=id;

ods output ParameterEstimates=est1 FitStatistics=fit1 ; 

run;

* with cure fraction, set gamma=1;

proc nlmixed data=six qpoints=5;

parms r01=0.25 r02=0.25 r03=0.24 r04=0.25 r05=0.25
	 	h01=0.1 h02=0.1 h03=0.1 h04=0.1 h05=0.1
alpha0=0 alpha1=-1 beta1=1.5 delta1=1 vara=.3 ;
bounds r01 r02 r03 r04 r05 h01 h02 h03 h04 h05 vara>=0;

eta= alpha0 + alpha1 * covar ; 		/* for probability of no recurrent event */
prob=1/(1+exp(-eta));

base_haz_r=r01 * event_r1 + r02 * event_r2 + r03 * event_r3 + r04 * event_r4 + r05 * event_r5 ;
cum_base_haz_r=r01 * dur_r1 + r02 * dur_r2 + r03 * dur_r3 + r04 * dur_r4 + r05 * dur_r5 ;

base_haz_d=h01 * event_d1 + h02 * event_d2 + h03 * event_d3 + h04 * event_d4 + h05 * event_d5 ;
cum_base_haz_d=h01 * dur_d1 + h02 * dur_d2 + h03 * dur_d3 + h04 * dur_d4 + h05 * dur_d5 ;

mu1= beta1 * covar +  a;			/* for recurrent event */
mu2= delta1 * covar + a;	/* for death event */

loglik1=-exp(mu1) * cum_base_haz_r;
loglik2=-exp(mu2) * cum_base_haz_d;


if event=1 then loglik = log(base_haz_r) + mu1; 		/*log likelihood for recurrent event */
if event=0 then do;
	if enum>1 then loglik= log(1-prob) + loglik1; 				/* log likelihood for subjects with at least 1 recurrent event */
	if enum=1 then loglik=log(prob + (1-prob) * exp(loglik1));					/*log likelihood for subjects without recurrent event */
	loglik=loglik + loglik2;
end;
if event=2 then do;
	if enum>1 then loglik= log(1-prob) + loglik1; 				/* log likelihood for subjects with at least 1 recurrent event */
	if enum=1 then loglik=log(prob + (1-prob) * exp(loglik1));					/*log likelihood for subjects without recurrent event */
	loglik=loglik + log(base_haz_d) + mu2 + loglik2;
end;


model id ~ general(loglik);

random a ~ normal(0, vara) subject=id;


ods output ParameterEstimates=est3 FitStatistics=fit3; 

run;
