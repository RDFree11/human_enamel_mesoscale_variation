function [V,FWHM,N] = PseudoVoigt(xx,Io,mean,FWHMg,FWHMl)
% PsuedoVoigt returns a psuedo-voigt profile,V, with approximate FWHM 
% defined based on Thompson et al. 1987 (THC formula) 
% with 4 defining parameters:
%
% xx     vector of x values
% Io     Integrated intensity of peak
% mean   mean parameter
% FWHMg  gaussian full width half max
% FWHMl  lorentzian full width half max

gauss = @(X,sigma) (1./(sqrt(pi).*sigma)).*exp((-X.^2)./(sigma.^2));
lorentz = @(X,gamma) (1/pi).*(gamma)./(X.^2+gamma.^2);
f = @(fg,fl) ((fg.^5)+2.69269.*(fg.^4).*(fl)+2.42843.*(fg.^3).*(fl.^2)+4.47163.*(fg.^2).*(fl.^3)+0.07842.*(fg).*(fl.^4)+(fl.^5)).^(1/5);
n = @(fg,fl) (1.36603*(fl./f(fg,fl))-0.47719*(fl./f(fg,fl)).^2+0.11116*(fl./f(fg,fl)).^3);

V = Io.*(n(FWHMg,FWHMl).*lorentz(xx-mean,f(FWHMg,FWHMl)/2)+(1-n(FWHMg,FWHMl)).*gauss(xx-mean,f(FWHMg,FWHMl)/(2*sqrt(log(2)))));
FWHM = f(FWHMg,FWHMl);
N = n(FWHMg,FWHMl);
end