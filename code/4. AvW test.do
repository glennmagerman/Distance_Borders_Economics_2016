*****************************************
*** Estimating AvW with other methods ***
*****************************************

clear*
version 12.0
set matsize 11000
set maxvar 32767, perm
cd "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/Literature/AvW (2003) Data/"

/* Not a perfect match with the original data, but results are close:
- US-US trade is corrupt file, I downloaded it again from the website: http://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=831&DB_Short_Name=CFS
- Still need to correct exchange rate of flows US-CA
- calculate internal distance
- keep only 40 regions
*/

******************************************
*** 1. Import data for 2-country model ***
******************************************

* distance
insheet using distance.txt,tab clear
drop v1-v5
gen from = _n
order from
reshape long v,i(from) j(to)
replace to=to-5
ren v dist
sort from to
save "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/Literature/AvW (2003) Data/stata/distance", replace

* names states
insheet using states.txt, tab clear
ren v1 state
ren v2 name
replace name = "District of Columbia" in 51
save "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/Literature/AvW (2003) Data/stata/states", replace

* trade US-US (1993) - in mio dollars
insheet using own_shipments.txt, tab clear
ren from name
merge n:1 name using "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/Literature/AvW (2003) Data/stata/states"
drop if _merge==1
drop _merge
drop name 
ren state from
ren to name
merge n:1 name using "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/Literature/AvW (2003) Data/stata/states"
drop if _merge==1
drop _merge
drop name
ren state to
order from to value
save "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/Literature/AvW (2003) Data/stata/usus", replace

* trade US-CA and CA-CA
insheet using trade1.txt, tab clear
drop v1-v2
ren v3 from
ren v4 to
ren v5 y1988
ren v6 y1989
ren v7 y1990
ren v8 y1991
ren v9 y1992
ren v10 y1993
ren v11 y1994
ren v12 y1995
reshape long y,i(from to) j(year)
ren y value
keep if year==1993
drop year
replace value=value/1000
save "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/Literature/AvW (2003) Data/stata/trade", replace

* nom GDP for 63 prov and states
insheet using nomGDP.txt, tab clear
gen year=_n + 1987
order year
reshape long v, i(year) j(state)
ren v gdp
keep if year ==1993
drop year
save "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/Literature/AvW (2003) Data/stata/gdp", replace

* exch rate
insheet using erate.txt, tab clear
gen year=_n + 1969
order year
ren v1 erate
keep if year==1993
// erate = .7751337

* GDP*rate
use "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/Literature/AvW (2003) Data/stata/gdp", clear
replace gdp=gdp*.7751337 if state>=52
save "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/Literature/AvW (2003) Data/stata/gdp", replace

/*
erate=erate[24,1]; /* US dollars per Canadian dollar in 1993 */
trade1=erate*trade1[.,8]/1000 ; /* transforms trade to millions of US dollars */
*/

**********************************
*** 2. Merge to dataset (1993) ***
**********************************
cd "~/Dropbox/PhD HUB-KUL/Papers/distance_borders/Literature/AvW (2003) Data/stata/"
use trade, clear
append using usus
sort from to
gen us_i=0
replace us_i=1 if from <=51
gen us_j=0 
replace us_j=1 if to<= 51
gen us=0
replace us=1 if us_i==1 & us_j==1

gen ca=0
replace ca=1 if us_i==0 & us_j==0
gen border=0
replace border=1 if us_i!=us_j
merge 1:1 from to using distance
keep if _merge==3
drop _merge
ren from state
merge n:1 state using gdp
keep if _merge==3
drop _merge
ren state from
ren gdp gdp_i
ren to state
merge n:1 state using gdp
keep if _merge==3
drop _merge
ren state to
ren gdp gdp_j
ren from i 
ren to j
gen lnv = ln(value)
gen lngdp_i = ln(gdp_i)
gen lngdp_j = ln(gdp_j)
gen lndist = ln(dist)
order i j lnv lngdp_i lngdp_j lndist border
save gravity_avw, replace

**********************
*** 3. Regressions ***
**********************

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

drop if i==j

eststo clear
eststo: reg lnv lngdp_i lngdp_j lndist us ca, robust cluster(cp)
eststo: reg lnv lngdp_i lngdp_j lndist us ca i_* j_*, robust cluster(cp)
eststo: reg lnv lngdp_i lngdp_j tay* us ca, robust cluster(cp)
eststo: ppml v  lndist us ca  i_* j_*, technique(nr 10 dfp 10 nr) cluster(cp)
eststo: glm  v  lndist us ca i_* j_*, family(gamma)  link(log) vce(robust) difficult technique(nr 10 dfp 10 nr)
esttab using avw_sample.csv, replace drop(i_* j_*) nogaps
