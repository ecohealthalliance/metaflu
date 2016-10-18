# Parameters



###Previous Model Parameters [@Hosseini_2013]:

Direct Transmission $\beta = 0.004$  

Recovery Rate $\upsilon = (0-0.167)$

Disease Mortality Rate $\alpha = (0 - 0.4)$

Environmental Uptake Rate $\nu = 0.001$

Virion Infectiousness $\phi = 1.96 * 10^{-4}$

Virus Degradation Rate $\eta = 0.14$

Virion Shedding Rate $\sigma = 10^5$

Movement Rate $\omega = 0.30$

###New Model Paramaters

__Poultry to Poultry Direct Transmission $\beta_{PP}$__
Originally adapted from duck study [@Sturm_Ramirez_2004] ~ but unclear how exactly. Existing studies / models deal with village- or flock-level transmission, not individuals. For example, in West Bengal[@Pandit_2013]: $1.00 * 10^{-8}$, originally taken from a study in wild birds [@Vaidya_2012]. 

__Wild Bird to Poultry Direct Transmission $\beta_{WP}$__  
Transmission between wild birds and poultry *communities* ($1.00 * 10^{-8}$) was based on previous study of avian influenza dynamics in *wild birds* [@Vaidya_2012]

__Poultry to Human Direct Transmission $\beta_{PH}$__  
Fitted parameters from models predicting the seasonality of cumulative cases of HPAI H5N1 [@TUNCER_2013], the estimated range for direct poultry to human transmission was $1.9 * 10^{-11}$ to $2.3 * 10^{-11}$.

__Recovery Rate $\upsilon$__  
Recovery for HPAI infected poultry does not seem to be commonly modeled -- instead either SI models [@WARD_2008; @TUNCER_2013] or SIC (Culled) [@Pandit_2013] are used.  

__Disease Mortality Rate $\alpha$__  
Originally given as a range from studies based on ducks / chickens. When re-calculated based on chicken-specific parameter paper[@Bouma_2009], $\alpha = 0.39$.

__Environmental Uptake Rate $\nu$__  
__Virion Infectiousness $\phi$__  
__Virus Degradation Rate $\eta$__  
__Virion Shedding Rate $\sigma$__  

Environmental transmission not deemed relevant for initial model ~ more relevant for fecal-oral LPAI transmission when no epidemic occurs. In a long-term model could be interesting to consider effects of cleaning. HPAI environmental persistence has been researched under laboratory conditions [@Wood_2010]. 

__Movement Rate $\omega$__

Original paper uses values from FAO Poultry Review in Egypt [@Unknown_2016], although exact decision on how to quantify $\omega$ not immediately clear. Once the method is understood, this resource may provide more specific information for Ethiopia: [@greycite54620]. Alternatively, FAO data may be relevant here. 




## References

<div id="refs"></div>
