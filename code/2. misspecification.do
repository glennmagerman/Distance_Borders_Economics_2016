********************************************************************************
*** Borders and Distance in the Gravity model and estimation specifications  ***
********************************************************************************

* Author: Glenn Magerman, University of Leuven, Naamsestraat 69, 3000 Leuven (Belgium). 
* email: glenn.magerman@kuleuven.be.

* First version: Oct 27, 2014.
* This version: April 15, 2015.

/* 
1. does sample size affect estimators? - random subsamples
2. Non-linearities
	- residual vs fitted plots -> we find non-linearities
	- residual vs regressor plots -> non-linearities are in distance
	- run model with polynomial in distance
3. Omitted variables
	- how to check?
4. Multicollinearity
	- VIF test on distance vs border
	- VIF test on rest of model
	- if problematic, run model without these vars (but try everything else first before dropping them)
*/
********************************************************************************

clear*
version 12.0
set matsize 11000
set maxvar 32767, perm
cd "~/Dropbox/PhD HUB-KUL/Papers/Degrees and Clustering/Clean Data"
use "Total_panel1998-2011", clear
cd "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/graphs"

*************************
*** 1. Random Samples ***
*************************
cd "~/Dropbox/PhD HUB-KUL/Papers/Degrees and Clustering/Clean Data"
use "Total_panel1998-2011", clear
	set seed 205
	sample 75
	capture drop i_* j_*
	qui tab i, gen(i_)
	qui tab j, gen(j_)
* generate Taylor approximations 
	foreach x in lndist contig lang_off colony { 			
	egen temp1 = mean(`x'), by(i)
	egen temp2 = mean(`x'), by(j)
	egen temp3 = mean(`x')
	gen tay`x' = `x' - temp1 - temp2 + temp3
	drop temp1 temp2 temp3
	}
*** Regressions
forvalues t = 1998(1)2011 {
	eststo:   reg lnv lngdp_i lngdp_j lndist contig lang_off colony rta wto_i wto_j if t==`t', robust cluster(cp)
	eststo:   reg lnv lndist contig lang_off colony rta i_* j_* if t==`t', robust cluster(cp)
	eststo:   reg lnv lngdp_i lngdp_j tay* rta wto_i wto_j if t==`t', robust cluster(cp)
	eststo:   ppml v lndist contig lang_off colony rta  i_* j_* if t==`t', technique(nr 10 dfp 10 nr) cluster(cp)
	eststo:   glm v lndist contig lang_off colony rta i_* j_* if t==`t', family(gamma) link(log) vce(robust) difficult technique(nr 10 dfp 10 nr)
}
	cd "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/tables"
	esttab using "random_75%.csv", not nostar plain r2 ar2 pr2 bic obslast nogaps replace drop(*i_* *j_*) 
	eststo clear
*** TSlines 
	cd "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/tables/robustness"
	import excel "sample_distance", clear sheet("25%") firstrow
*** distance
	sort year
	line ols lsdv bb ppml gpml year, xlabel(1998(1)2011) ylabel(-1.6(0.2)-0.4) ///
		legend(label(1 "OLS") label(2 "LSDV") label(3 "BB") label(4 "PPML") label(5 "GPML")) ///
		xtitle("Year") ytitle("Coefficient")
	cd "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/graphs"
	graph export "robustness_75%.pdf", replace
clear 
exit 
exit


*************************
*** 2. Residual plots ***
*************************
cd "~/Dropbox/PhD HUB-KUL/Papers/Degrees and Clustering/Clean Data"
use "Total_panel1998-2011", clear
cd "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/graphs"
	capture drop i_* j_*
	qui tab i, gen(i_)
	qui tab j, gen(j_)
* generate Taylor approximations 
	foreach x in lndist contig lang_off colony { 			
		egen temp1 = mean(`x'), by(i)
		egen temp2 = mean(`x'), by(j)
		egen temp3 = mean(`x')
		gen tay`x' = `x' - temp1 - temp2 + temp3
		drop temp1 temp2 temp3
	}

*** 2.1 Residuals vs fitted plots
	* do all on LSDV, since residual plots look identical across methods
* crosssection
	reg lnv lndist contig lang_off colony rta i_* j_* if t==2005, robust cluster(cp)
	predict err1, resid
	predict fit1, xb
	lpoly err1 fit1, noscatter ci legend(label(2 "Polynomial Smoother")) yline(0)
	graph export "lpoly_cross_section.pdf", replace

* pooled
	reg lnv lndist contig lang_off colony rta i_* j_* t_*, robust cluster(cp)
	predict err2, resid
	predict fit2, xb
	lpoly err2 fit2, noscatter ci legend(label(2 "Polynomial Smoother")) yline(0)
	graph export "lpoly_pooled.pdf", replace
	
* panel
	egen it = group(i t)
	qui tab it, gen(it_)
	egen jt = group(j t)
	qui tab jt, gen(jt_)
	reghdfe lnv lndist contig lang_off colony rta, absorb(fe_i=i.it fe_j=i.jt) vce(cluster cp) 
	predict err4, resid
	predict fit4, xb
	lpoly err4 fit4, noscatter ci legend(label(2 "Polynomial Smoother")) yline(0)
	graph export "lpoly_panel.pdf", replace

* robust regression (WLS) (since clustering does not show up in residual analysis)
	rreg lnv lndist contig lang_off colony rta i_* j_* t_*
	predict err3, resid
	predict fit3, xb
	lpoly err3 lnv, noscatter ci legend(label(2 "Polynomial Smoother")) yline(0)
	graph export "lpoly_test3.pdf", replace
// polynomial relationship stays in error across all specifications + a lot of noise in lower values

*** 2.2 Residual vs regressor plots
	* crosssection
	reg lnv lndist contig lang_off colony rta i_* j_* if t==2005, robust cluster(cp)
	predict err4, resid
	foreach x in lnv lndist contig lang_off colony rta {
		lpoly err4 `x', noscatter ci legend(label(2 "Polynomial Smoother")) yline(0)
		graph export temp`x'.pdf, replace
	}
// all dummies are perfectly on 0, disturbance in distance and trade value. distance converges to 0 very quickly

*** 2.3 Add polynomial to model 
	gen lndist2 = (lndist)^2
	gen lndist3 = (lndist)^3
	
	reg lnv lngdp_i lngdp_j lndist lndist2 lndist3 contig lang_off colony rta if t==2005, robust cluster(cp)
	predict err6, resid
	predict fit6, xb
	lpoly err6 fit6, noscatter ci legend(label(2 "Polynomial Smoother")) yline(0)
	graph export "lpoly_lndistpoly3.pdf", replace
// improves residuals, but destroys the model.


*** 2.4. subgroup differences
	capture drop err* fit*
	forvalues i = 1(1)5 {
		keep if continent_i=="`i'"
		keep if continent_j=="`i'"
		reg lnv lndist contig lang_off colony rta i_* j_* if t==2005, robust cluster(cp)
		predict err`i', resid
		predict fit`i', xb
		err`i' fit`i', noscatter ci legend(label(2 "Polynomial Smoother")) yline(0)
		graph export "lpoly_continent`i'.pdf", replace
	}
	// same pattern, across all continents!

****************************
*** 3. Multicollinearity ***
****************************
*** test on OLS, since dummies take too much time for VIF and cannot do subset
reg lnv lngdp_i lngdp_j lndist contig lang_off colony rta, robust cluster(cp)
estat vif
// all vif are around 1.0-1.5 so no problems at all... (<= 5-6 is ok)

clear
exit
exit
