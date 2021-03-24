function [Class_Weight,Weights] = Func_Weights(Class,Lengths,y,z)

% Weights (g) at start of age year (W)
%
% Jess Hopf
% Nov 2012
%
% P.leopardus
%
% Based on Chan et al 2007
%
% Inputs; Age = vector of ages
%         Lenghts = vector of lengths 
%         y = length-weight parameter
%         z = length-weight parameter
% Outputs; Weights = vector of weights
%          Age_Weights = 'table' of age and weights


     for i=1:length(Lengths)
        Weights(i,:) = y*(Lengths(i)).^z;
     end
     
    Weights;
    Class_Weight = [Class,Weights];

end
