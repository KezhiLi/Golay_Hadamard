function X_mat=test_mseq_fft(powerval,odd_even)

% powerval: powervalue required to generate m sequence; 
% odd_even: 1: odd lenghth sequence, L=2^(powerval+1)+1; otherwise, 
%  
%Output parameter: X_mat: L\times L OSTM; 
% Example: X_mat=test_mseq_fft(6,1);

aa=mseq(2,powerval);

% Odd-length sequence;
if odd_even==1,
a_exp=[1;aa;aa(length(aa):-1:1)];
else 
% Even-length sequence;
a_exp=[-1;aa;1;aa(length(aa):-1:1)];
end

La=length(a_exp);
bb=ifft(a_exp);
 %bb_imag=imag(bb);
 %max(abs(bb_imag(:)))
b_mag=abs(bb);
%[norm(a_exp),norm(bb)]



 X_mat=fft(eye(La))*diag(a_exp)*ifft(eye(La));
 X_mat=real(X_mat);
 % The remaining commands aim to check the magnitude of X_mat; 
% [max(b_mag),min(b_mag)]
% 1/sqrt(length(X_mat))

