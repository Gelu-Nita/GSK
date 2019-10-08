;+
; NAME:
;       GSK
;
; PURPOSE:
;       This function computes the lower and upper detection thresholds
;       for the Generalized Spectral Kurtosis estimator corresponding to a given false alarm (PFA) probability level
;
; CATEGORY:
;       Numerical Analysis.
;
; CALLING SEQUENCE:
;       Result = GSK(M,N=N,d=d,PFA=PFA)
;
; INPUTS:
;      M: Off-board accumulation length
;      N: On-board accumulation length
;      d: Shape factor of Gamma distribution according to wich the raw samples being accumulated are assumed to be distributed
;         ChiSqr (Time Domain) d=0.5
;         Exponential (Spectral Domain) d=1
;         NOTE: if not povided, d=1 (exponential distribution) is assumed by default 
;
; KEYWORD PARAMETERS:
;       PFA:	false alarm probability level. If not provided the default PFA=0.0013499 (3-sigma normal PDF) is used.

; OUTPUTS:
;       This function returns the lower and upper thresholds in the array form [lower,upper]
;
; REFERENCES: Gelu Nita and Dale Gary, The Generalized Spectral Kurtosis Estimator,
;             Monthly Notices of the Royal Astronomical Society: Letters, Volume 406, Issue 1, 1 July 2010, Pages L60–L64, 
;             https://doi.org/10.1111/j.1745-3933.2010.00882.x
;             
;             Gelu M Nita, Jack Hickish, David MacMahon, Dale E Gary
;             EOVSA Implementation of a Spectral Kurtosis Correlator for Transient Detection and Classification
;             Journal of Astronomical Instrumentation, Volume 5, Issue 4, December 2016, Pages 1641009
;             https://doi.org/10.1142/S2251171716410099
;
;
; MODIFICATION HISTORY:
;           Written 2010 April by Gelu M. Nita (gnita@njit.edu).
;
;-
FUNCTION gsk_upper_root, X
 common share_gsk, m1,m2,m3,m4,p
 if (-(m3-2*m2^2)/m3+x[0])/(m3/2/m2) lt 0 then return, 0
 return,abs((1-igamma(4*(m2^3)/(m3^2),(-(m3-2*m2^2)/m3+x[0])/(m3/2/m2)))-p)
END

FUNCTION gsk_lower_root, X
 common share_gsk, m1,m2,m3,m4,p
 if (-(m3-2*m2^2)/m3+x[0])/(m3/2/m2) lt 0 then return, 0
 return,abs(igamma(4*(m2^3)/(m3^2),(-(m3-2*m2^2)/m3+x[0])/(m3/2/m2))-p)
END



function GSK,M,N=N,d=d,PFA=PFA
     common share_gsk, m1,m2,m3,m4,p
     on_error,1

     ;Define default parameters
     if n_elements(N) eq 0 then N=1
     if n_elements(d) eq 0 then d=1
 	 if n_elements(M) eq 0 then M=6104
 	 if n_elements(PFA) eq 0 then PFA=0.0013499

     ;Force double precision
	 p=double(pfa)
	 M=double(M)
	 Nd=double(N*d)
   PFA=double(PFA)

	 ;Compute GSK moments according to equation 9 (used in threshold computation)
	 m1=1d
	 m2=(2*( M^2)* Nd *(1 + Nd))/((-1 + M) *(6 + 5* M* Nd + (M^2)*( Nd^2)))
	 m3=(8*( M^3)* Nd* (1 + Nd)* (-2 + Nd* (-5 + M *(4 + Nd))))/(((-1 + M)^2)* (2 + M* Nd) *(3 +M* Nd)* (4 + M* Nd)* (5 + M* Nd))
	 m4=(12*( M^4)* Nd* (1 + Nd)* (24 + Nd *(48 + 84* Nd + M *(-32 + Nd *(-245 - 93 *Nd + M* (125 + Nd* (68 + M + (3 + M)* Nd)))))))/$
	    (((-1 + M)^3)* (2 + M* Nd)* (3 + M *Nd)* (4 + M* Nd) *(5 + M *Nd)* (6 + M* Nd)* (7 + M *Nd))

	 ;Compute normalized moments entering Pearson criterion (not used below)
	 beta1=(8* (2 + M *Nd)* (3 + M* Nd)* (-2 + Nd* (-5 + M* (4 + Nd)))^2)/((-1 + M)* Nd* (1 + Nd) *((4 + M *Nd)^2)* (5 + M *Nd)^2)
	 beta2=(3* (2 + M *Nd) *(3 + M* Nd)* (24 + Nd* (48 + 84* Nd + M *(-32 + Nd *(-245 - 93* Nd + M *(125 + Nd* (68 + M + (3 + M)* Nd)))))))/((-1 +M)* Nd *(1 + Nd)* (4 + M* Nd)* (5 + M *Nd)* (6 + M* Nd)* (7 + M *Nd))

     ;Compute Pearson criterion according to equation 11 (not used below)
     k=2*m2*(3*(m2^2)-m4)+3*(m3^2)

     ;Compute Type III parameters defined by equation 19
     delta=m1-2*(m2^2)/m3
     beta=4*(m2^3)/(m3^2)
     alpha=m3/(2*m2)

     ;Compute fourth moment error according to equation 21
     err4=abs((3* beta *(2 + beta)* (alpha^4)/m4-1)*100)

     ;Compute the thresholds according to equation 20

     upper = NEWTON([1], 'gsk_upper_root',/double,tolf=1e-8)
     lower = NEWTON([1], 'gsk_lower_root',/double,tolf=1e-8)

     return,[lower,upper]
end