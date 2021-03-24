function [Age_Length,Lengths] = Func_Length(Age,Linf,K,T0)

% Length at start of age year (L)
%
% Jess Hopf
% Nov 2012
%
% P.leopardus
%
% Based on Chan et al 2007
% von Bertalanffy growth curve
%
% Inputs; Age = vector of ages
%         Linf = avg legth at asymtote  
%         K = growth rate 
%         T0 = age at which length is zero
% Outputs; Lengths = vector of lengths
%          Age_Length = 'table' of age and lengths


    for i=1:length(Age)
        Lengths(i,:) = Linf*(1-exp(-K*(Age(i)-T0)));
    end
    
    Lengths;
    Age_Length=[Age,Lengths];

end

