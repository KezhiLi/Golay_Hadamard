

close all; clear; clc;

% --- load image ---

F = wavread('guita_note.wav');

F = F;

N=512;
[m n] = size(F);

const = N/2; % signal length
M = N/4;
k=  M/4;

F_block_num=floor(m/const);
F_block_num=F_block_num+sign(m-const*F_block_num);
F_0 = [F;zeros(1,(const*F_block_num-m))'];  % add zeros 




% --- # of subsamples ---
%nSamples = round(0.20*m*n); % sample rate 20%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gauss_opRD = randn(2*M,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gauss_op2 = randn(1,M+N-1);
Gauss_op3 = randn(1,M+N-1);
Gauss_op_row=complex(Gauss_op2(1:N),Gauss_op3(1:N));
Gauss_op_col=complex(Gauss_op2(N:M+N-1),Gauss_op3(N:M+N-1));
Gauss_op1 = toeplitz(Gauss_op_row(1:M), Gauss_op_row);

Gauss_op = [real(Gauss_op1); imag(Gauss_op1)];  % for real matrix
% for ii=1:N,
%     Gauss_op(:,ii)=Gauss_op1(:,ii);%/norm(Gauss_op1(:,i));
% end
Gauss_op2 = [];
Gauss_op3 = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thesize=N;   % Rauhut's code, convolution based on Rademacher sequence in fourier domain

para= sign(randn(1,thesize));  
      aa225=fft(para);
  %aa225=para;
    
    aa226=zeros(1,N);
    aa226(1)=aa225(1);
 for ii=2:N;
     aa226(ii)=aa225(N-ii+2);
 end
matr5=toeplitz(aa226,aa225)/thesize;
%matr5=toeplitz(aa226,aa225);


rp22=randperm(N);

matr15=matr5(rp22(1:M),:);

matr5=[];
%matr1=matr(1:1:1*M,:);
matr15=matr15*sqrt(N/M);
% Circulant Determinisic Matrix
WHT_op =matr15;

WHT_op = [real(WHT_op); imag(WHT_op)];  % for real matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root=3;                        % Zadoff code
zadoff_seq=zadoff(root, N);
fistrow_zad2=fft(zadoff_seq);

fistrow_zad1=zeros(1,N);
fistrow_zad1(1)=fistrow_zad2(1);
 for ii=2:N;
     fistrow_zad1(ii)=fistrow_zad2(N-ii+2);
 end
 fistrow_zad=toeplitz(fistrow_zad1,fistrow_zad2);
 rp3=randperm(N);
 
Zad_op=fistrow_zad(rp3(1:M),:);
Zad_op=Zad_op*sqrt(N/M);

Zad_op = [real(Zad_op); imag(Zad_op)];  % for real matrix
%for i=1:N,
%    CDM_op(:,i)=CDM_op(:,i)/norm(CDM_op(:,i));
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [a,b] = generate_golay(log2(thesize));                 % Golay sequence

fistrow_golay2=fft(a);
fistrow_golay1=zeros(1,N);
fistrow_golay1(1)=fistrow_golay2(1);
 for ii=2:N;
     fistrow_golay1(ii)=fistrow_golay2(N-ii+2);
 end
 fistrow_golay=toeplitz(fistrow_golay1,fistrow_golay2);
 
%  fistrow_golay=hadamard(thesize)*diag(a)*dftmtx(thesize)/sqrt(thesize);
 rp8=randperm(N);
 
golay_op=fistrow_golay(rp8(1:M),:);
golay_op=golay_op*sqrt(N/M);
 CDM_op=golay_op;

 CDM_op = [real( CDM_op); imag( CDM_op)];  % for real matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F_GaussRD = zeros(length(F_0),1);
F_Gauss = zeros(length(F_0),1);
F_Zad = zeros(length(F_0),1);
F_WHT = zeros(length(F_0),1);
F_CDM = zeros(length(F_0),1);

        x0=ifft(F_0);
        
        x0= x0.* (abs(x0)>0.002);
            

disp('Now starts the simulations');
for jj = 1:F_block_num

        jj
        
        ifft_F0_blk = x0((jj-1)*const+1:jj*const); 
        real_ifft_F0_blk = [real(ifft_F0_blk); imag(ifft_F0_blk)]; 
        y0 = Gauss_opRD*real_ifft_F0_blk;
        y1 = Gauss_op*real_ifft_F0_blk;
        y2 = Zad_op*real_ifft_F0_blk;
        y3 = WHT_op*real_ifft_F0_blk;
        y4 = CDM_op*real_ifft_F0_blk;
        
        % observation
%         ifft_F0_blk = ifft(F_0((jj-1)*const+1:jj*const)); %calculate the ifft of the block signal
%         real_ifft_F0_blk = [real(ifft_F0_blk); imag(ifft_F0_blk)]; % construct the real block signal with twice length
%         y0 = Gauss_opRD*real_ifft_F0_blk;
%         y1 = Gauss_op*real_ifft_F0_blk;
%         y2 = Zad_op*real_ifft_F0_blk;
%         y3 = WHT_op*real_ifft_F0_blk;
%         y4 = CDM_op*real_ifft_F0_blk;

        nn = N;
        tol0=10;
        tol1=10;
        tol2=10;
        tol3=10;
        tol4=10;
        sigma = 1e-10;
        step_size =5;
        lambda=0.0001;
        rel_tol=0.005;
         mu=0.0001;
        
        alp0 = CSRec_SP(k, Gauss_opRD, y0);
      %[alp0,Out0] = fpc_AS(N,Gauss_opRD,y0,mu,[]);
   %   [alp0 iter_num0] = SAMP(y0, Gauss_opRD, step_size, sigma);
        
   %     alp1 = omp(x, Gauss_op, y1, k, sigma);
        %[alp1,x_debias,objective,times,debias_start,mses]= GPSR_Basic1(y1,Gauss_op,tau);
      alp1 = CSRec_SP(k, Gauss_op, y1);
        %[alp1,d1]=cosamp(Gauss_op*dct(N),y1,k,tol1);
  %      [alp1 iter_num1] = SAMP(y1, Gauss_op, step_size, sigma);
         %alp1 = l1dantzig_pd(Gauss_op'*y1, Gauss_op, [], y1, 4e-3, 5e-3); % l1 minimization
        %alp1 =l1_ls(Gauss_op,y1,lambda,rel_tol);
        %[alp1,Out1] = fpc_AS(N,Gauss_op,y1,mu,[]);
        
   %     alp2 = omp(x, Zad_op, y2, k, sigma);
        %[alp2,x_debias,objective,times,debias_start,mses]= GPSR_Basic1(y2,Gold_op,tau);
      alp2 = CSRec_SP(k, Zad_op, y2);
      %[alp2,d2]=cosamp(Zad_op*dct(N),y2,k,tol1);
    %    [alp2 iter_num2] = SAMP(y2, Zad_op, step_size, sigma);
        %alp2 = l1dantzig_pd(Zad_op'*y1, Zad_op, [], y2, 4e-3, 5e-3); % l1 minimization
      %alp2 =l1_ls(Zad_op,y2,lambda,rel_tol);
      %[alp2,Out2] = fpc_AS(N,Zad_op,y2,mu,[]);
      
   %     alp3 = omp(x, WHT_op, y3, k, sigma);
        %[alp3,x_debias,objective,times,debias_start,mses]= GPSR_Basic1(y3,WHT_op,tau);
      alp3 = CSRec_SP(k, WHT_op, y3);
        %[alp3,d3]=cosamp(WHT_op*dct(N),y3,k,tol1);
     %   [alp3 iter_num3] = SAMP(y3, WHT_op, step_size, sigma);
        %alp3 = l1dantzig_pd(WHT_op'*y3, WHT_op, [], y3, 4e-3, 5e-3); % l1 minimization
       %alp3 = l1_ls(WHT_op,y3,lambda,rel_tol);
       %[alp3,Out3] = fpc_AS(N,WHT_op,y3,mu,[]);
        
   %     alp4 = omp(x, CDM_op, y4, k, sigma);
        %[alp4,x_debias,objective,times,debias_start,mses]= GPSR_Basic1(y4,CDM_op,tau);
      alp4 = CSRec_SP(k, CDM_op, y4);
        %[alp4,d4]=cosamp(CDM_op*dct(N),y4,k,tol1);
      %  [alp4 iter_num4] = SAMP(y4, CDM_op, step_size, sigma);
         %alp4 = l1dantzig_pd(CDM_op'*y4, CDM_op, [], y4, 4e-3, 5e-3); % l1 minimization
        %alp4 = l1_ls(CDM_op,y4,lambda,rel_tol);
        %[alp4,Out4] = fpc_AS(N,CDM_op,y4,mu,[]);
        




                alp0 = complex(alp0(1:const),alp0(const+1:2*const));
                alp1 = complex(alp1(1:const),alp1(const+1:2*const));
                alp2 = complex(alp2(1:const),alp2(const+1:2*const));
                alp3 = complex(alp3(1:const),alp3(const+1:2*const));
                alp4 = complex(alp4(1:const),alp4(const+1:2*const));

              F_GaussRD((jj-1)*const+1:jj*const) = alp0;
            F_Gauss((jj-1)*const+1:jj*const) = alp1;
            F_Zad((jj-1)*const+1:jj*const) = alp2;
            F_WHT((jj-1)*const+1:jj*const) = alp3;
             F_CDM((jj-1)*const+1:jj*const) = alp4;

%             F_GaussRD((jj-1)*const+1:jj*const) = fft(alp0);
%             F_Gauss((jj-1)*const+1:jj*const) = fft(alp1);
%             F_Zad((jj-1)*const+1:jj*const) = fft(alp2);
%             F_WHT((jj-1)*const+1:jj*const) = fft(alp3);
%              F_CDM((jj-1)*const+1:jj*const) = fft(alp4);

end 

F_GaussRD = fft(F_GaussRD);
F_Gauss = fft(F_Gauss);
F_Zad = fft(F_Zad);  % random sign
F_WHT = fft(F_WHT);  % FZC 
F_CDM = fft(F_CDM);  % Golay

F_GaussRD = F_GaussRD(1:length(F));
F_Gauss = F_Gauss(1:length(F));
F_Zad = F_Zad(1:length(F));  % FZC
F_WHT = F_WHT(1:length(F));  % random sign
F_CDM = F_CDM(1:length(F));  % Golay

snr(F_GaussRD, F)
snr(F_Gauss, F)
snr(F_Zad, F)
snr(F_WHT, F)
snr(F_CDM, F)


figure,
subplot(6,1,1), plot(F)
subplot(6,1,2),plot(real(F_GaussRD));
subplot(6,1,3),plot(real(F_Gauss));
subplot(6,1,4),plot(real(F_Zad));
subplot(6,1,5),plot(real(F_WHT));
subplot(6,1,6),plot(real(F_CDM));


