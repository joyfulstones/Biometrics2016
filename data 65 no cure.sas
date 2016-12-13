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



%macro cure();

%do ii=599 %to 600;

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

* No cure fraction;


proc nlmixed data=six qpoints=5 MAXITER=5000 instep=1e-2;

parms r01=0.04 r02=0.05 r03=0.07 r04=0.08 r05=0.09
	 	h01=0.01 h02=0.02 h03=0.02 h04=0.02 h05=0.03
		beta1=-.24 delta1=-.16 vara=4.4;
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


data data65.est1_all;
set data65.est1_all est1;
run;

data data65.fit1_all;
set data65.fit1_all fit1;
run;


%end;

%mend;

%cure();
data result;
set _null_;
run;
%macro organize_result(varname);
data two;
set data65.est1_all;
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
%organize_result(vara);

%organize_result(r01);
%organize_result(r02);
%organize_result(r03);
%organize_result(r04);
%organize_result(r05);
%organize_result(h01);
%organize_result(h02);
%organize_result(h03);
%organize_result(h04);
%organize_result(h05);


data _null_;
set result;
file "c:\liulei\paper\paper11\data65\AGQresult1.txt";
put beta1 delta1  vara r01 r02 r03 r04 r05 h01 h02 h03 h04 h05 se_beta1 se_delta1 se_vara se_r01 se_r02 se_r03 se_r04 se_r05 se_h01 se_h02 se_h03 se_h04 se_h05; 
run;
