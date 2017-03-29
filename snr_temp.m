function x = snr_temp(sig, ref)
% written by Junfeng Yang, Wotao Yin, and Yin Zhang

sig=reshape(sig, 1,length(sig));
ref=reshape(ref, 1,length(ref));

dv = sig-ref;
x = 10*log10((ref*ref')/(dv*dv'));