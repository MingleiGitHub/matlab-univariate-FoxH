% the test script calculates the PDF of the Rayleigh distribution with 0 mean
% and unit variate in the form of univariate Fox H function and in its 
% simplified representation.
% The difference between the two calculations is one the level of 1e-14.
%
% For details regarding the univariate Fox H function and its related
% issues, please refer to the book of "H Function: theory and
% application"
%
% For details between other PDF/CDF/MGF of fading channels
% and their Fox H representations, please refer to the paper and reference
% therein.
%
% ****************************************************************************************************************
% Citation: If you use this software or any (modified) part of it, please cite it as:
% Minglei You, Hongjian Sun, Jing Jiang and Jiayi Zhang, ``Unified framework for
% the effective rate analysis of wireless communication systems over MISO fading channels'',
% IEEE Transactions on Communications, 65(4), pp. 1775-1785, April, 2017.
% https://ieeexplore.ieee.org/abstract/document/7782422
% https://arxiv.org/pdf/1612.03882.pdf
% ****************************************************************************************************************
% The complex gamma function credits to Paul Godfrey, see http://winnie.fit.edu/~gabdo/gamma.txt. 
% ****************************************************************************************************************

evaluation_sample_range = 0.1:0.1:3;
Rayleigh_pdf_Fox_H = zeros(1, length(evaluation_sample_range));
Rayleigh_pdf_SimpleFunction = Rayleigh_pdf_Fox_H;

the_counter=0;
for z=evaluation_sample_range
    % syms z
    the_counter=the_counter+1;
    O=[1 0 0 1];
    P={[1],[1],[],[0.5],[],[0.5]};
    
    
    Rayleigh_pdf_Fox_H(the_counter)=func_Fox_H_univariate(z, O, P);

    Rayleigh_pdf_SimpleFunction(the_counter)=2*z*exp(-z^2);
end

diff_result=Rayleigh_pdf_Fox_H-Rayleigh_pdf_SimpleFunction

figure(1)
plot(Rayleigh_pdf_Fox_H, '-*k')
hold on
plot(Rayleigh_pdf_SimpleFunction, '-or')


