* Author: Glenn Magerman, ECARES (ULB)
* email: glenn.magerman@ulb.ac.be.
* Version: March 19, 2016.
/*
________________________________________________________________________________
This do-file is part of the paper "Borders and Distance" 
by Glenn Magerman, Zuzanna Studnicka and Jan Van Hove (2016)
Install following ado's first: esttab, ppml
All remaining errors are mine. If you have suggestions or find mistakes, 
please send me an e-mail.
________________________________________________________________________________
*/

clear*
version 12.0
set matsize 11000
set scheme  s1color
global folder "~/Dropbox/Research/Published/distance_borders"
global bilateral "lndist contig lang_off colony rta"
capture mkdir $folder/Output
capture mkdir $folder/Output/Graphs
capture mkdir $folder/Output/Tables

******************
*** 0. Prelims ***
******************
// Code to explicit creation of dummies and Taylor approximated variables
cd $folder/Data
use gravity_98_11, clear
* Convergence issues PPML with data in dollars, rescale monetary values to 1000s
	foreach x in v gdp_i gdp_j {
		replace `x' = `x'/1000
	}
* generate logs	
	foreach x in v dist gdp_i gdp_j {
		gen ln`x' = ln(`x')
	}
* generate dummies
	qui tab i, gen(i_)
	qui tab j, gen(j_)
	qui tab t, gen(t_)
	egen cp = group(i j)
* generate Taylor approximations 
	foreach x in lndist contig lang_off colony { 			
		egen temp1 = mean(`x') if v>0 & v<., by(i)
		egen temp2 = mean(`x') if v>0 & v<., by(j)
		egen temp3 = mean(`x') if v>0 & v<.
		gen tay`x' = `x' - temp1 - temp2 + temp3 if v>0 & v<.
		drop temp1 temp2 temp3
	}
save gravity_98_11b, replace

****************************************
*** 1. Distance coefficients by year ***
****************************************
cd $folder/Data
use gravity_98_11b, clear
* Regressions
	eststo clear
	forvalues t = 1998(1)2011 {
		eststo: reg  lnv lngdp_i lngdp_j $bilateral wto_i wto_j if t==`t', robust cluster(cp)
		eststo: reg  lnv $bilateral i_* j_* if t==`t', robust cluster(cp)
		eststo: reg  lnv lngdp_i lngdp_j tay* rta wto_i wto_j if t==`t', robust cluster(cp)
		eststo: ppml v   $bilateral  i_* j_* if t==`t', technique(nr 10 dfp 10 nr) cluster(cp)
		eststo: glm  v   $bilateral i_* j_* if t==`t', family(gamma) link(log) /// 
		vce(robust) difficult technique(nr 10 dfp 10 nr)
	}
	cd $folder/Output/Tables
	esttab using "cross_sections.csv", not nostar plain r2 ar2 pr2 bic obslast nogaps drop(*i_* *j_*) replace 
* Time Series graph distance coefficients
// manually collected distance coefficients from cross_sections.csv
	cd $folder/Output/Tables
	insheet using "dist_coeff.csv", clear 
	tsset t
	tsline dist_*, xlabel(1998(1)2011) ylabel(-2(0.5)0) legend(label(1 "OLS") ///
		label(2 "LSDV") label(3 "BB") label(4 "PPML") label(5 "GPML")) ///
		xtitle("Year") ytitle("Coefficient")
	cd $folder/Output/Graphs
	graph export "tsline_dist.eps", replace
* Test for significant difference between estimates in 1998 and 2011. 
// equivalent to t-test
cd $folder/Data
use gravity_98_11b, clear
	gen dummy=.
	replace dummy=1 if t==2011
	replace dummy=0 if t==1998 
	reg lnv lndist contig lang_off colony rta i_* j_* dummy, robust cluster(cp)
	// dummy is significant: there is a significant difference in both years.
	
******************************************
*** 2. Border effects by world regions ***
******************************************
cd $folder/Data
use gravity_98_11b, clear
*** 1. generate dummies for regions
	ren i iso_3
	merge n:1 iso_3 using "countries_regions.dta", keepusing(region)
	drop if _merge==2
	drop _merge
	ren (iso_3 j region) (i iso_3 region_i)
	merge n:1 iso_3 using "countries_regions.dta", keepusing(region)
	drop if _merge==2
	drop _merge
	ren (iso_3 region) (j region_j)
	qui tab region_i, gen(region_i_)
	qui tab region_j, gen(region_j_)
*** 2. generate dummies for intra-regional trade	
	gen ASEAN = 0
	replace ASEAN = 1 if region_i=="ASEAN" & region_j=="ASEAN"
	gen Africa = 0
	replace Africa = 1 if region_i=="Africa (sub-sahara)" & region_j=="Africa (sub-sahara)"
	gen Asia = 0
	replace Asia = 1 if region_i=="Asia" & region_j=="Asia"
	gen Europe = 0
	replace Europe = 1 if region_i=="Europe" & region_j=="Europe"
	gen North_Africa = 0
	replace North_Africa = 1 if region_i=="North Africa" & region_j=="North Africa"
	gen North_America = 0
	replace North_America = 1 if region_i=="North America" & region_j=="North America"
	gen Pacific = 0
	replace Pacific = 1 if region_i=="Pacific" & region_j=="Pacific"
	gen South_America = 0
	replace South_America = 1 if region_i=="South America" & region_j=="South America"	
	global regions "Europe North_America South_America Asia ASEAN North_Africa Africa Pacific" 
*** 3. generate cluster var on continent for clustering SE
	gen cont_clus=0
	replace cont_clus=1 if ASEAN==1
	replace cont_clus=2 if Africa==1
	replace cont_clus=3 if Asia==1
	replace cont_clus=4 if Europe==1
	replace cont_clus=5 if North_Africa==1
	replace cont_clus=6 if North_America==1
	replace cont_clus=7 if Pacific==1
	replace cont_clus=8 if South_America==1
*** 4. pooled estimates
	eststo clear
	eststo: reg  lnv lngdp_i lngdp_j $bilateral wto_i wto_j $regions t_*, robust cluster(cont_clus)
	eststo: reg  lnv $bilateral $regions i_* j_* t_*, robust cluster(cont_clus)
	eststo: reg  lnv lngdp_i lngdp_j tay* rta wto_i wto_j $regions t_*, robust cluster(cont_clus)
	eststo: ppml v   $bilateral $regions i_* j_* t_*, technique(nr 10 dfp 10 nr) cluster(cont_clus)
	eststo: glm  v   $bilateral $regions i_* j_* t_*, family(gamma) link(log) difficult technique(nr 10 dfp 10 nr) vce(robust)
	cd $folder/Output/Tables
	esttab using "gravity_regions.csv", r2 ar2 pr2 se bic obslast nogaps replace drop(*i_* *j_* *t_*)
	
*** 5. Time series evolution of intra-continental trade (only LSDV)
	eststo clear
	forvalues t = 1998(1)2011 {
		eststo:   reg lnv $bilateral $regions i_* j_* if t==`t', robust cluster(cont_clus)
	}
	cd $folder/Output/Tables
	esttab using "continent_coefficients.csv", se nostar plain r2 ar2 pr2 bic obslast nogaps replace drop(*i_* *j_*) 
*** 6. Time Series graphs coefficients
// manually collected continent coefficients from continent_coefficients.csv
cd $folder/Output/Tables
import excel "continent_coefficients.xlsx", clear firstrow
	sort year
	line coeff_* year, xlabel(1998(1)2011) ///
		legend(label(1 "Europe") label(2 "North America") /// 
		label(3 "South America") label(4 "Asia") label(5 "ASEAN") ///
		label(6 "North Africa") label(7 "Africa") label(8 "Pacific")) ///
		xtitle("Year") ytitle("Coefficient")
	cd $folder/Output/Graphs
	graph export "tsline_cont.eps", replace	
	
**************************	
*** 3. US-Canada trade ***
**************************
cd $folder/Data
use gravity_avw, clear
/*
foreach x in exports gdp_i gdp_j {
	replace `x'=`x'*1000
	}
	drop lnv
	gen lnv =ln(exports)
*/	
	qui tab i, gen(i_)
	qui tab j, gen(j_)
	egen cp = group(i j)
	foreach x in lndist { 			
		egen temp1 = mean(`x'), by(i)
		egen temp2 = mean(`x'), by(j)
		egen temp3 = mean(`x')
		gen tay`x' = `x' - temp1 - temp2 + temp3
		drop temp1 temp2 temp3
	}
	eststo clear
	eststo: reg lnv lngdp_i lngdp_j lndist adjacency us ca, robust cluster(cp)
	eststo: reg lnv lndist adjacency us ca i_* j_*, robust cluster(cp)
	eststo: reg lnv lngdp_i lngdp_j tay* adjacency us ca, robust cluster(cp)
	eststo: ppml v  lndist adjacency us ca i_* j_*, technique(nr 10 dfp 10 nr) cluster(cp)
	eststo: glm  v lndist adjacency us ca, family(gamma)  link(log) difficult technique(nr 10 dfp 10 nr) vce(robust) 
	cd $folder/Output/Tables
	esttab using us_can_results.csv, r2 ar2 se bic obslast nogaps replace drop(i_* j_*) nogaps
	* no adjacency
	eststo clear
	eststo: reg lnv lngdp_i lngdp_j lndist us ca, robust cluster(cp)
	eststo: reg lnv lndist us ca i_* j_*, robust cluster(cp)
	eststo: reg lnv lngdp_i lngdp_j tay* us ca, robust cluster(cp)
	eststo: ppml v  lndist us ca i_* j_*, technique(nr 10 dfp 10 nr) cluster(cp)
	eststo: glm  v  lndist us ca, family(gamma)  link(log) difficult technique(nr 10 dfp 10 nr) vce(robust) 
	cd $folder/Output/Tables
	esttab using us_can_results_noadj.csv, r2 ar2 se bic obslast nogaps replace drop(i_* j_*) nogaps

clear
exit 
exit
