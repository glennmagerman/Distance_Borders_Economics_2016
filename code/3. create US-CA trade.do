* Author: Glenn Magerman, KU Leuven, Naamsestraat 69, 3000 Leuven (Belgium). 
* email: glenn.magerman@kuleuven.be.
* Version: March 20, 2016.

/*
________________________________________________________________________________
This do-file is part of the paper "Borders and Distance" 
by Glenn Magerman, Zuzanna Studnicka and Jan Van Hove (2016)
This file estimates coefficients for US-Canada trade
Install following ado's first: esttab, ppml
All remaining errors are mine. If you have suggestions or find mistakes, 
please send me an e-mail.
________________________________________________________________________________
*/

clear*
version 12.0
set matsize 11000
set scheme  s1color
global folder "~/Dropbox/Research/Papers/distance_borders"
global bilateral "lndist contig lang_off colony rta"

******************************************
*** 1. Import data for 2-country model ***
******************************************
cd $folder/Data/AvW
* distance
insheet using Distance.txt,tab clear // balanced pairs 63*63 = 3969 - OK
	drop v1-v5
	gen from = _n
	order from
	reshape long v,i(from) j(to)
	replace to=to-5
	ren v dist
	sort from to
	* follow AvW: create internal distance as 1/4 distance closest region (Wei)
	forvalues x=1/63 {
	 sum(dist) if from==`x' & dist>1 
	 replace dist=r(min)*0.25 if from==to & from==`x' 
	}
	save distance, replace
	// 1-51 US, 52-63 CAN
	// checked some distances, seems to be ok

* exch rate
insheet using erate.txt, tab clear
	gen year=_n + 1969
	ren v1 erate
	keep if year==1993
	di erate // .77513373 = US$ per CAN$
	global erate .77513373
	
* trade US-CA and CA-CA // only 40 regions!!
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
	ren y v
	keep if year==1993
	drop year
	replace v = v* $erate /1000 // 1000s of CAD$ to mio US$
	save trade, replace 
/*		
* trade US-US (1993) - in mio US $
*insheet using own_shipments_matrix_transposed.txt, tab clear
insheet using own_shipments_matrix.txt, tab clear
	gen from = _n
	order from
	reshape long v,i(from) j(to)
	replace v=(3025/5846)*v // explained by AvW in appendix
	ren v v
	save usus, replace
*/	

* trade US-US (1993) - in mio US$
// original AvW file is corrupt
// use CFS data for 1997 (only I can find), and deflate to 1993 values
// http://www.transtats.bts.gov/DL_SelectFields.asp?Table_ID=831&DB_Short_Name=CFS	
// delete DC trade (wrong order and never as destination)
import excel shipments_1997.xlsx, sheet("Blad1") firstrow case(lower) clear
	reshape wide v,i(state) j(state_name) string
	gen from = _n
	order from
	reshape long
	egen to = group(state_name)
	*replace v=(3025/5846)*v // explained by AvW in appendix
	replace v = v*0.89  // deflate 1997 values to 1993 values
	drop if state=="District of Columbia"
	keep from to value
	ren v v
	save usus, replace

* nom GDP for 63 prov and states - transform to mio 1993 US$ 
insheet using nomGDP.txt, tab clear
	gen year=_n + 1987
	order year
	reshape long v, i(year) j(state)
	ren v gdp
	keep if year ==1993
	drop year
	replace gdp=gdp*0.8
	save gdp, replace // all in 1000s of CAN$ // correct, checked with AvW readme file
	// 1-51 US, 52-63 CAN
	
* adjacency info
	insheet using adjacency.csv, clear	// file "doubles" from-to and back
	keep from to
	gen adjacency = 1
	save adjacency, replace	

* state/province names and ids 1-63
insheet using states.txt, clear
	ren v1 state
	ren v2 name
	save states,replace	// 1-51 US, 52-63 CAN

**********************************
*** 2. Merge to dataset (1993) ***
**********************************
cd $folder/Data/AvW
use trade, clear
	append using usus
	fillin from to 
	drop _fillin
	replace v=0 if v==.
	* add GDP
	ren from state
	merge n:1 state using gdp, nogen keep(match)
	merge n:1 state using states, nogen keep(match)
	ren state from
	ren name state_i
	ren gdp gdp_i
	ren to state
	merge n:1 state using gdp, nogen keep(match)
	merge n:1 state using states, nogen keep(match)
	ren state to
	ren name state_j
	ren gdp gdp_j
	drop if from==51 | to==51 // DC
	drop if from==57 | to==57 // NW Territories and 
	drop if from==63 | to==63 // Yukon territory 
	* AvW 30 + 10 subsample
	drop if from==2 | to==2
	drop if from==4 | to==4
	drop if from==6 | to==6
	drop if from==7 | to==7
	drop if from==8 | to==8
	drop if from==11 | to==11
	drop if from==15 | to==15
	drop if from==16 | to==16
	drop if from==24 | to==24
	drop if from==27 | to==27
	drop if from==28 | to==28
	drop if from==31 | to==31
	drop if from==36 | to==36
	drop if from==37 | to==37
	drop if from==39 | to==39
	drop if from==40 | to==40
	drop if from==41 | to==41
	drop if from==44 | to==44
	drop if from==48 | to==48
	drop if from==50 | to==50
	sort from to
	
	* gen within US and Canada dummies
	gen us_i=0
	replace us_i=1 if from <= 50
	gen us_j=0 
	replace us_j=1 if to <= 50
	gen us=0
	replace us=1 if us_i==1 & us_j==1
	gen ca=0
	replace ca=1 if us_i==0 & us_j==0 
	* gen adjacency dummy
	merge 1:1 from to using adjacency
	drop if _merge==2
	drop _merge
	replace adjacency = 0 if adjacency==.
	* add distance
	merge 1:1 from to using distance
	keep if _merge==3
	drop _merge
	ren from i 
	ren to j
	drop if i==j
	gen lnv = ln(v)
	gen lngdp_i = ln(gdp_i)
	gen lngdp_j = ln(gdp_j)
	gen lndist = ln(dist)
	order i j lnv lngdp_i lngdp_j lndist adjacency
	cd $folder/Data
	save gravity_avw, replace
	
clear
exit 
exit	
