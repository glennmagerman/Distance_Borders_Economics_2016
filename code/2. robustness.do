* Author: Glenn Magerman, KU Leuven, Naamsestraat 69, 3000 Leuven (Belgium). 
* email: glenn.magerman@kuleuven.be.
* Version: March 19, 2016.

/* 
________________________________________________________________________________
1. does sample size affect estimators? - random subsamples --> drop this
2. Non-linearities
	- residual vs fitted plots -> we find non-linearities
	- residual vs regressor plots -> non-linearities are in distance
	- run model with polynomial in distance
3. Multicollinearity
	- VIF test on distance vs border
	- VIF test on rest of model
	- if problematic, run model without these vars (try everything else first before dropping)
________________________________________________________________________________
*/

clear*
version 12.0
set matsize 11000
set scheme  s1color
global folder "~/Dropbox/Research/Papers/distance_borders"
global bilateral "lndist contig lang_off colony rta"
capture mkdir $folder/Output/Tables/Robustness

* drop lowest GDP countries to compare to older studies when using CEPII BACI data
	qui sum gdp_i, d
	keep if gdp_i > r(p25)
	qui sum gdp_j, d
	keep if gdp_j > r(p25)


*********************
*** 1. Subsamples ***
*********************
*** Keep largest GDP countries
cd $folder/Data
use gravity_98_11b, clear
foreach y in 25 50 75 {
	qui sum gdp_i, d
	keep if gdp_i > r(p`y')
	qui sum gdp_j, d
	keep if gdp_j > r(p`y')
		capture drop i_* j_* tay*
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
		* Regressions
		eststo clear
		forvalues t = 1998(1)2011 {
			capture eststo: reg  lnv lngdp_i lngdp_j $bilateral wto_i wto_j if t==`t', robust cluster(cp)
			capture eststo: reg  lnv $bilateral i_* j_* if t==`t', robust cluster(cp)
			capture eststo: reg  lnv lngdp_i lngdp_j tay* wto_i wto_j if t==`t', robust cluster(cp)
			capture eststo: ppml v   $bilateral  i_* j_* if t==`t', technique(nr 10 dfp 10 nr) cluster(cp)
			capture eststo: glm  v   $bilateral i_* j_* if t==`t', family(gamma) /// 
			link(log) vce(robust) difficult technique(nr 10 dfp 10 nr)
		}
		cd $folder/Output/Tables
		esttab using "GDPpercentiles_`y'.csv", not nostar plain r2 ar2 pr2 bic obslast nogaps replace drop(*i_* *j_*) 
	}

* Time Series distance 
// manually collected coefficients for distance in csv file

	foreach x in 25 50 75 {
		cd $folder/Output/Tables/Robustness
		import excel "subsample_GDP", clear sheet("`x'") firstrow
		tsset t
		tsline dist_*, xlabel(1998(1)2011) ylabel(-2(0.5)0) legend(label(1 "OLS") ///
			label(2 "LSDV") label(3 "BB") label(4 "PPML") label(5 "GPML")) ///
			xtitle("Year") ytitle("Coefficient")
		cd $folder/Output/Graphs
		graph export "robustness_`x'.pdf", replace
	}	
	
*************************
*** 2. Residual plots ***
*************************
cd $folder/Data
use gravity_98_11b, clear

*** 2.1 Residuals vs fitted plots
// All on LSDV, since residual plots look identical across methods
* crosssection
	reg lnv $bilateral i_* j_* if t==2005, robust cluster(cp)
	predict err1, resid
	predict fit1, xb
	lpoly err1 fit1, noscatter ci legend(label(2 "Polynomial Smoother")) yline(0)
	cd $folder/Output/Graphs
	graph export "lpoly_cross_section.eps", replace

* pooled
	reg lnv lndist contig lang_off colony rta i_* j_* t_*, robust cluster(cp)
	predict err2, resid
	predict fit2, xb
	lpoly err2 fit2, noscatter ci legend(label(2 "Polynomial Smoother")) yline(0)
	cd $folder/Output/Graphs
	graph export "lpoly_pooled.eps", replace

*** 2.2 Residual vs regressor plots
	* crosssection
	reg lnv $bilateral i_* j_* if t==2005, robust cluster(cp)
	predict err5, resid
	foreach x in $bilateral {
		lpoly err5 `x', noscatter ci legend(label(2 "Polynomial Smoother")) yline(0)
		graph export temp`x'.eps, replace
	}
// all dummies are perfectly on 0, disturbance in distance, but converges to 0 quickly

*** 2.3 Add polynomial to model 
	gen lndist2 = (lndist)^2
	gen lndist3 = (lndist)^3
	reg lnv lngdp_i lngdp_j lndist lndist2 lndist3 contig lang_off colony rta if t==2005, robust cluster(cp)
	predict err6, resid
	predict fit6, xb
	lpoly err6 fit6, noscatter ci legend(label(2 "Polynomial Smoother")) yline(0)
	graph export "lpoly_lndistpoly3.eps", replace
// no improvement, interpret joint coefficients? also, is not structural

****************************
*** 3. Multicollinearity ***
****************************
*** test on OLS, since dummies take too much time for VIF and cannot do subset
reg lnv lngdp_i lngdp_j $bilateral, robust cluster(cp)
estat vif
// all vif are around 1.0-1.5 so no problems at all... (<= 6-10 is ok)

clear
exit
exit
