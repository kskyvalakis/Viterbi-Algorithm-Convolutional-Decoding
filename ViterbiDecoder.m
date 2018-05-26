function [m_est] = ViterbiDecoder(y,N,epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT :
% y : Coded sequence
% N : Length of uncoded message
% epsilon : Error probability of the BSC
%
%------------------------------------------------------------
% OUTPUT :
% y_est : estimated decoded sequence of y
%
%------------------------------------------------------------
% Variable Naming Convention :
%
% 1 --> 00
% 2 --> 01
% 3 --> 10
% 4 --> 11
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializations
w_e = log((1-epsilon)/epsilon);
cf = zeros(4,N);
w = zeros(4,N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward-Pass

for i=1:N-1
    coded = y(2*i-1:2*i);
    if(i==1)
        w(1,1) = double(coded(1)==0)*w_e + double(xor(coded(1),coded(2))==0)*w_e;
        w(2,1) = double(coded(1)==0)*w_e + double(xor(coded(1),coded(2))==1)*w_e;
        w(3,1) = double(coded(1)==1)*w_e + double(xor(coded(1),coded(2))==0)*w_e;
        w(4,1) = double(coded(1)==1)*w_e + double(xor(coded(1),coded(2))==1)*w_e;
        cf(1,1) = 0;
        cf(2,1) = 0;
        cf(3,1) = 0;
        cf(4,1) = 0;
    else
        [w1, cf(1,i)] = max([w(1,i-1) -inf w(3,i-1) -inf]);
        w(1,i) = double(coded(1)==0)*w_e + double(xor(coded(1),coded(2))==0)*w_e + w1;
        
        [w2, cf(2,i)] = max([w(1,i-1) -inf w(3,i-1) -inf]);
        w(2,i) = double(coded(1)==0)*w_e + double(xor(coded(1),coded(2))==1)*w_e + w2;
        
        [w3, cf(3,i)] = max([-inf w(2,i-1) -inf w(4,i-1)]);
        w(3,i) = double(coded(1)==1)*w_e + double(xor(coded(1),coded(2))==0)*w_e + w3;
        
        [w4, cf(4,i)] = max([-inf w(2,i-1) -inf w(4,i-1)]);
        w(4,i) = double(coded(1)==1)*w_e + double(xor(coded(1),coded(2))==1)*w_e + w4;
    end
end

[w1, cf(1,N)] = max([w(1,N-1) -inf w(3,N-1) -inf]);
w(1,N) = double(y(2*N-1)==0)*w_e + w1;

[w3, cf(3,N)] = max([-inf w(2,N-1) -inf w(4,N-1)]);
w(3,N) = double(y(2*N-1)==1)*w_e + w3;

[~, fn] = max([w1 -inf w3 -inf]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backward-Pass

path = zeros(1,N);
path(N) = fn;
for i=N-1:-1:1
    path(i) = cf(path(i+1),i+1);
end

y_est = zeros(1,2*N-1);
for i=1:N-1
    if(path(i)==1)
        y_est(2*i-1:2*i) = [0 0];
    elseif(path(i)==2)
        y_est(2*i-1:2*i) = [0 1];
    elseif(path(i)==3)
        y_est(2*i-1:2*i) = [1 0];
    else
        y_est(2*i-1:2*i) = [1 1];
    end
end
if(path(N)==1)
    y_est(2*N-1) = 0;
else
    y_est(2*N-1) = 1;
end

m_est = zeros(1,N);
for i=1:N
    m_est(i) = y_est(2*i-1);
end

end
