;+
; NAME:
;       GSK
;
; PURPOSE:
;       This function computes the the Generalized Spectral Kurtosis estimator corresponding to a set of
;       accumulated powers (s1) and squared powers (s2) over an accumulation length M
;
; CATEGORY:
;       Numerical Analysis.
;
; CALLING SEQUENCE:
;       Result = get_sk(s1,s2,M,d=d)
;
; INPUTS:
;      s1: Accumalted power samples
;      s2: Accumulated squared power samples
;      M:  Accumulation length
;      d:  Shape factor of Gamma distribution according to wich the raw samples being accumulated are assumed to be distributed
;          ChiSqr (Time Domain) d=0.5
;          Exponential (Spectral Domain) d=1
;          N onboard accumulaions of time domain power samples d=0.5*N
;          N onboard accumulations of spectral domain power samples d=1*N
;      NOTE: if not povided, d=1 (exponential distribution) is assumed by default    
;
; KEYWORD PARAMETERS:
;     

; OUTPUTS:
;       This function returns the GSK estimator
;
; REFERENCES: Gelu Nita and Dale Gary, 
;             Monthly Notices of the Royal Astronomical Society: Letters, Volume 406, Issue 1, 1 July 2010, Pages L60–L64, 
;             https://doi.org/10.1111/j.1745-3933.2010.00882.x
;
;
; MODIFICATION HISTORY:
;           Written April 2010 by Gelu M. Nita (gnita@njit.edu)
;
;-

function get_sk,s1,s2,m,d=d
default,d,1
return,(m*d+1)*(m*s2/s1^2-1)/(m-1)
end