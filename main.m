close all;
clear;
clc;

N = 128;
P = 10^5;
epsilon = [1/100 1/20 1/10 1/8 1/6 1/5];
BER = zeros(1,length(epsilon));

m = rand(1,N) < 0.5;
b = zeros(1,2*N-1);
for k=1:N-1
    b(2*k-1) = m(k);
    b(2*k) = xor(m(k),m(k+1));
end
b(2*N-1) = m(N);

for i=1:length(epsilon)
    for j=1:P
        y = bsc(b,epsilon(i));
        m_est = ViterbiDecoder(y,N,epsilon(i));
        BER(i) = BER(i) + sum(m_est~=m);
    end
end

BER = BER./(N*P);

figure, semilogy(epsilon,BER), xlabel('Epsilon values'), ylabel('Bit Error Probability')