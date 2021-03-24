function [Class_Fecundity,Fecunds] = Func_fecunds(Class,Lengths,c,d)
% Fecundity at length at the start of age year (f)
%
% Jess Hopf
% Nov 2012
%
% Based on Chan et al 2007
%
% Inputs; Class = vector of ages/stages/lengths etc
%        Lengths = vector of lengths       
%        c = length-fecundity parameters (Samoilys 00)
%        d = length-fecundity parameters (Samoilys 00)
%
% Outputs; Fecunds = vector of fecundities
%          Class_Fecundity = 'table' of classes and fecunditites

    for i=1:length(Class)
        Fecunds(i,:) = c*(Lengths(i).^d);
    end        
    
    Fecunds;
    Class_Fecundity = [Class,Fecunds];
    
end

