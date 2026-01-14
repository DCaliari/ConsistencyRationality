clear

import excel "DATI.xlsx", sheet("Variables") firstrow


gen patient = 0 if _n < 146
replace patient =1 if beta > 0.99

gen impatient = 0 if _n < 146
replace impatient =1 if beta < 0.91


gen neutral = 0 if _n < 146
replace neutral =1 if rho < 0.2

gen averse = 0 if _n < 146
replace averse =  1 if rho > 2.2


************************
** Generate variables **
************************


egen _Duration_XYZW = std(Duration_XYZW)
egen _WARP_L_MAIN = std(WARP_L_MAIN)
egen _WARP_D_MAIN = std(WARP_D_MAIN)
egen _Duration_D = std(Duration_D)

egen _Duration_ABCD = std(Duration_ABCD)
egen _Raven = std(Raven_W)
egen _RavenTime = std(RavenTime)
egen _Duration_L = std(Duration_L)

gen Mon_Risk = 0 if _n < 146
replace Mon_Risk = 1 if FOSDmoney == 1 | FOSDprob == 1

gen Mon_Risk_C = FOSDmoney + FOSDprob


*****************
**** Table 3  ***
*****************

reg WARP_D_MAIN _Raven _RavenTime  _Duration_XYZW if impatient==1 | patient==1 | del_time==1 , robust
outreg2 using c.xls, replace ctitle(WARP time)
reg WARP_D_MAIN _Raven _RavenTime  _Duration_XYZW if impatient==0 & patient==0 & del_time==0 , robust
outreg2 using c.xls, append ctitle(WARP time)


reg WARP_L_MAIN _Raven _RavenTime  _Duration_ABCD if averse==1 | neutral==1 | deliberate==1 , robust
outreg2 using c.xls, replace ctitle(WARP time)
reg WARP_L_MAIN _Raven _RavenTime  _Duration_ABCD if averse==0 & neutral==0 & deliberate==0 , robust
outreg2 using c.xls, append ctitle(WARP time)

*****************
**** Table 4  ***
*****************


reg WARP_D_MAIN Impatience, robust
outreg2 using c.xls, replace ctitle(WARP time)
reg WARP_D_MAIN Impatience  impatient patient del_time , robust
outreg2 using c.xls, append ctitle(WARP time)
reg WARP_D_MAIN Impatience  impatient patient del_time _Duration_D , robust
outreg2 using c.xls, append ctitle(WARP time)

reg WARP_D_MAIN  Mon_time , robust
outreg2 using c.xls, append ctitle(WARP time)
reg WARP_D_MAIN  Mon_time impatient patient del_time , robust
outreg2 using c.xls, append ctitle(WARP time)
reg WARP_D_MAIN  Mon_time impatient patient del_time _Duration_D , robust
outreg2 using c.xls, append ctitle(WARP time)

reg WARP_L_MAIN Mon_Risk  , robust
outreg2 using c.xls, replace ctitle(WARP risk)
reg WARP_L_MAIN Mon_Risk   averse neutral deliberate , robust
outreg2 using c.xls, append ctitle(WARP risk)
reg WARP_L_MAIN Mon_Risk    averse neutral _Duration_L deliberate , robust
outreg2 using c.xls, append ctitle(WARP risk)

reg WARP_L_MAIN  SOSD , robust
outreg2 using c.xls, append ctitle(WARP risk)
reg WARP_L_MAIN  SOSD   averse neutral deliberate , robust
outreg2 using c.xls, append ctitle(WARP risk)
reg WARP_L_MAIN  SOSD   averse neutral _Duration_L deliberate , robust
outreg2 using c.xls, append ctitle(WARP risk)


	
