%% version 20200128 by Ming
% this is a function to evaluate the univariate fox-h function using matlab
% See definiation of the Fox_H in the book of "H Function: theory and
% application"

% ****************************************************************************************************************
% Citation: If you use this software or any (modified) part of it, please cite it as:
% Minglei You, Hongjian Sun, Jing Jiang and Jiayi Zhang, ``Unified framework for
% the effective rate analysis of wireless communication systems over MISO fading channels'',
% IEEE Transactions on Communications, 65(4), pp. 1775-1785, April, 2017.
% https://ieeexplore.ieee.org/abstract/document/7782422
% https://arxiv.org/pdf/1612.03882.pdf
% ****************************************************************************************************************

function Fox_H_value=func_Fox_H_univariate(z, O, P)

% parameters defination here
length=50;
leave_from_real_axis=0.1;
% disamble the O and P
m=O(1);
n=O(2);
p=O(3);
q=O(4);

k=P{1};
c=P{2};
a=P{3};
b=P{4};
A=P{5};
B=P{6};
Fox_H_value_temp=real(integral(@(s)mnpq_part_generation(c*z,s,m,n,p,q,a,b,A,B), leave_from_real_axis-1i*length, leave_from_real_axis+1i*length,'RelTol', 1.e-7,'ArrayValued',true)/(2*pi*1i));
Fox_H_value=Fox_H_value_temp*k;

end

function mnpq_part=mnpq_part_generation(z,s,m,n,p,q,a,b,A,B)
% initialization for the integration theta(s)
m_part= 1;
n_part= 1;
p_part= 1;
q_part= 1;
% create the integration theta(s)
for iii=1:1:m
    m_part= m_part.*gammaz(b(iii)-B(iii).*s);
end

for iii=1:1:n
    n_part= n_part.*gammaz(1-a(iii)+A(iii).*s);
end

for iii=n+1:1:p
    p_part= p_part.*gammaz(a(iii)-A(iii).*s);
end

for iii=m+1:1:q
    q_part= q_part.*gammaz(1-b(iii)+B(iii).*s);
end

mnpq_part=z.^s.*m_part.*n_part./p_part./q_part;

end
