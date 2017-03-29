close all; clear; clc;

% --- load image ---
F = phantom(128);
%F=imread('phantom256.png');
%F = rgb2gray(F);
%F=double(F)/256;
%F=(double(F)-min(min(double(F))))/(max(max(double(F)))-min(min(double(F))));

[m n] = size(F);



% --- # of subsamples ---

%
% --- PF parameters ---
% importnat parameters are .maxItr .beta3 .relchg_tol
% for REAL data, set: opts.bsymm = true, opts.bComplex = true
aTV = 1e-8;
opts = [];
    opts.maxItr = 3000; % max # of iterations
    opts.gamma = 1.618; % noiseless choice = 1.618
    opts.beta1 = 100; % TV ADM parameter
    opts.beta2 = 10; % L1 ADM parameter, not used if aL1=0
    opts.beta3 = 100*aTV; % CU = V ADM parameter
    opts.bsymm = true; % use conj_symm ifft2? If set yes, intermediate U are real
    opts.relchg_tol = 1e-5; % stopping tolerance based on relative change
    opts.real_sol = false; % return solution of real type (i.e., U = real(U))
    opts.bPrint = true; %true; % screen display option
    opts.normalize = false; % New parameter/data normalization was added to make
    opts.bComplex = true;

% --- call PC solver ---
aL1 = 0; WT = []; W = []; % turn off L1-wavelets; TV-minimization is kept
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root=1;                        % Zadoff code
zadoff_seq=zadoff(root, m*n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
loop=1; %change number of iteration here
tt=1;

Toe_deter = zeros(2,tt);
Toe_ran=zeros(2,tt);
Golay=zeros(2,tt);
FZC=zeros(2,tt);


for pp=1:tt
nSamples = round((0.105-0.005*pp)*m*n); % sample rate 10% when pp=1
result_loop1=zeros(2,loop);
result_loop2=zeros(2,loop);
result_loop3=zeros(2,loop);
result_loop4=zeros(2,loop);
for qq=1:loop;
    
% --- generate subsample mask ---
picks = randsample(m*n,nSamples);   % random sampling -kezhi
picks = sort(picks);
if (picks(1) ~= 1); picks(1) = 1; end  % make sure 1st sample is taken
pick = false(m,n); pick(picks) = true; % generate logical matrix "pick"

 picks_deter = 1:1:1*length(picks);  % deterministic sampling -kezhi
 picks_deter = picks_deter';
if (picks_deter(1) ~= 1); picks_deter(1) = 1; end  % make sure 1st sample is taken
pick_deter = false(m,n); pick_deter(picks_deter) = true;


%OTF = exp(rand(m,n)*(2*pi*1i));   % random phase (random convolution)
OTF = sign(randn(m,n));            % random binary sequence
%OTF = reshape(generate_golay(log2(m*n)),128,128);
%OTF = reshape(zadoff_seq,128,128);
OTF = conjugate_symmetrize(OTF); OTF = OTF./abs(OTF);
PSF = otf2psf(OTF,[m n]); if ~isreal(PSF); error('PSF is not real.'); end
CB = ifft2(OTF.*fft2(F),'symmetric'); CB = CB(pick_deter);  %derministic sampling
[U1,Out] = RecPF_Circ_Eq(m,n,aTV,aL1,pick_deter,PSF,CB,2,opts,WT,W,range(F(:)),F); % this uses equality constraints
U1 = abs(U1); 
result_loop1(1,qq)=snr(U1);
result_loop1(2,qq)=psnr(U1,F);

%OTF = exp(rand(m,n)*(2*pi*1i));
OTF = sign(randn(m,n));
OTF = conjugate_symmetrize(OTF); OTF = OTF./abs(OTF);
PSF = otf2psf(OTF,[m n]); if ~isreal(PSF); error('PSF is not real.'); end
CB = ifft2(OTF.*fft2(F),'symmetric'); CB = CB(pick);
[U2,Out] = RecPF_Circ_Eq(m,n,aTV,aL1,pick,PSF,CB,2,opts,WT,W,range(F(:)),F); % this uses equality constraints
U2 = abs(U2); 
result_loop2(1,qq)=snr(U2);
result_loop2(2,qq)=psnr(U2,F);


OTF = reshape(zadoff_seq,128,128);
OTF = conjugate_symmetrize(OTF); OTF = OTF./abs(OTF);
PSF = otf2psf(OTF,[m n]); if ~isreal(PSF); error('PSF is not real.'); end
CB = ifft2(OTF.*fft2(F),'symmetric'); CB = CB(pick);
[U3,Out] = RecPF_Circ_Eq(m,n,aTV,aL1,pick,PSF,CB,2,opts,WT,W,range(F(:)),F); % this uses equality constraints
U3 = abs(U3); 
result_loop3(1,qq)=snr(U3);
result_loop3(2,qq)=psnr(U3,F);

OTF = reshape(generate_golay(log2(m*n)),128,128);
OTF = conjugate_symmetrize(OTF); OTF = OTF./abs(OTF);
PSF = otf2psf(OTF,[m n]); if ~isreal(PSF); error('PSF is not real.'); end
CB = ifft2(OTF.*fft2(F),'symmetric'); CB = CB(pick);
[U4,Out] = RecPF_Circ_Eq(m,n,aTV,aL1,pick,PSF,CB,2,opts,WT,W,range(F(:)),F); % this uses equality constraints
U4 = abs(U4); 
result_loop4(1,qq)=snr(U4);
result_loop4(2,qq)=psnr(U4,F);

end
Toe_deter(1,pp)=mean(result_loop1(1,:));
Toe_ran(1,pp)=mean(result_loop2(1,:));
FZC(1,pp)=mean(result_loop3(1,:));
Golay(1,pp)=mean(result_loop4(1,:));
Toe_deter(2,pp)=mean(result_loop1(2,:));
Toe_ran(2,pp)=mean(result_loop2(2,:));
FZC(2,pp)=mean(result_loop3(2,:));
Golay(2,pp)=mean(result_loop4(2,:));
end

%figure, hold on, 
figure, imshow(U1); %title('original');
figure, imshow(U2); %title('original');
figure, imshow(U3); %title('original');
figure, imshow(U4); %title('original');

Toe_deter
Toe_ran
FZC
Golay

% RC
% Golay
% FZC 
% RB
%plot(RC,'r');hold on;plot(Golay,'g'); plot(FZC,''); plot(OSTM);




% figure;hold on;
% plot(0.1:-0.005:0.07, Toe-deter,'-b+');
% plot(0.1:-0.005:0.07, Toe-ran, '-k*');
% plot(0.1:-0.005:0.07, FZC, '-rs');
% plot(0.1:-0.005:0.07, Golay, '-mo');
% 
% grid
% %success_Gauss
% legend('Toeplitz deterministic','Toeplitz random ' ,'FZC sequence', 'Golay sequence');
% %legend('i.i.d Gaussian Operator','Partial DFT','CDM Matrix','Zadoff-Chu Code');
% %legend('i.i.d Gaussian Operator','Partial WHT ','OSTM');
% xlabel('Rate of K/N');
% ylabel('SNR');




% figure(1); clf;
% subplot(221); imshow(F,[]); title('original');
% subplot(222); imshow(pick,[]); title(sprintf('sample mask, %4.1f%%',100*nnz(pick)/numel(pick)));
% subplot(223); imshow(U,[]); title(sprintf('recovery, SNR %4.1f',snr(U)));
% subplot(224); imshow(F-U,[]); colorbar; title('error');
