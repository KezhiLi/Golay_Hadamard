%function testgold();
% test gold code (size 32x1024);
% parameters
clear

N = 1024; % signal length
M = 128;

% gaussian measurement matrix
 Gauss_op0 = randn(M,N);
 Gauss_op1 = randn(M,N);
 Gauss_opRD= complex(Gauss_op0,Gauss_op1);
%Gauss_opRD = randn(2*M,N);
% normalization of the Gaussian matrix;


% Gauss_op1 = randn(1,N+M-1);
% Gauss_op2 = toeplitz(Gauss_op1(1:M),Gauss_op1(M:(M+N-1)));
% Gauss_op = Gauss_op2;

% Gauss_op=zeros(M,N);
% Gauss_op2=randn(1,N);
% Gauss_op3=randn(1,N);
% Gauss_op_row=complex(Gauss_op2(1:N),Gauss_op3(1:N));
% Gauss_op1 = toeplitz(Gauss_op_row, Gauss_op_row);
% rp00=randperm(N);
% Gauss_op=Gauss_op1(rp00(1:M),:);

Gauss_op2 = randn(1,M+N-1);
Gauss_op3 = randn(1,M+N-1);
Gauss_op_row=complex(Gauss_op2(1:N),Gauss_op3(1:N));
Gauss_op_col=complex(Gauss_op2(N:M+N-1),Gauss_op3(N:M+N-1));
Gauss_op1 = toeplitz(Gauss_op_row(1:M), Gauss_op_row);
for ii=1:N,
    Gauss_op(:,ii)=Gauss_op1(:,ii);%/norm(Gauss_op1(:,i));
end
Gauss_op2 = [];
Gauss_op3 = [];



%Gold_op = Goldcode32x1024/sqrt(M);

% % partial Hadamard with random permuatation
% %%%%%
% %WHT = hadamard(N)/sqrt(N);
% %WHT = hadamard(N);
% %%%%%
% WHT = dftmtx(N);
% for i=1:N,
%     WHT(:,i)=WHT(:,i)/norm(WHT(:,i));
% end
% %%%%%
% q2 = randperm(N-1);
% WHT_op=[WHT(1,:);WHT(q2(1:(M-1))+1,:)];
% WHT=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  thesize=N;   %random convolution
% %  
% % [a,b] = generate_golay(log2(thesize)-1);
% aaa=exp(2*rand(1,(thesize/2-1))*pi );  % random phase for random convolution
%     
%     AA=zeros(1,(thesize/2+1));
%     AA(1)=sign(randn(1));
%     AA(thesize/2+1)=sign(randn(1));
%     AA(2:thesize/2)=aaa;
%     
%     B=zeros(1,thesize);
%     B(1:(thesize/2+1))=AA(1:(thesize/2+1));
%     for count=1:(thesize/2-1);
%         B(thesize/2+1+count)=conj(AA(thesize/2+1-count));
%     end
%     para=B;
% 
% matr5=conj(func_dctmtx(thesize))*diag(para)*func_dctmtx(thesize);    
% %    aa225=fft(para);
% %matr5=toeplitz(aa225,aa225)/thesize;
% for ii=1:N,
%     matr5(ii,:)=matr5(ii,:)/norm(matr5(ii,:));
% end
% 
% rp22=randperm(N);
% 
% matr15=matr5(rp22(1:M),:);
% 
% matr5=[];
% %matr1=matr(1:1:1*M,:);
% matr15=matr15*sqrt(N/M);
% % Circulant Determinisic Matrix
% WHT_op =matr15;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% para=zeros(1,N);   %square phase code
% for pp=1:N;
%     para(pp)=exp(-(2*pi*sqrt(-1)/N)*(pp-1)^2);
% end
% aa22=fft(para);
% matr=toeplitz(aa22,aa22)*dftmtx(N);  % This is F* diag() but not F* diag() F
% 
% rp2=randperm(N);
% matr1=matr(rp2(1:M),:)/norm(matr(1,:));
% %matr1=matr(1:1:1*M,:);
% matr1=matr1*sqrt(N/M);
% % Circulant Determinisic Matrix
% CDM_op=matr1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matra=zeros(N,N);
% for ii=1:N
%     for jj=1:N
%         matra(ii,jj)=exp(-(pi*sqrt(-1)/N)*(ii-jj)^2);
%     end
% end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root=1;                        % Zadoff code
zadoff_seq=zadoff(root, N);
fistrow_zad2=fft(zadoff_seq);

fistrow_zad1=zeros(1,N);
fistrow_zad1(1)=fistrow_zad2(1);
 for ii=2:N;
     fistrow_zad1(ii)=fistrow_zad2(N-ii+2);
 end
 fistrow_zad=toeplitz(fistrow_zad1,fistrow_zad2)/N;
 rp3=randperm(N);
 
Zad_op=fistrow_zad(rp3(1:M),:);
Zad_op=Zad_op*sqrt(N/M);

%for i=1:N,
%    CDM_op(:,i)=CDM_op(:,i)/norm(CDM_op(:,i));
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [a,b] = generate_golay(log2(thesize));                 % Golay sequence

% fistrow_golay2=fft(a);
% fistrow_golay1=zeros(1,N);
% fistrow_golay1(1)=fistrow_golay2(1);
%  for ii=2:N;
%      fistrow_golay1(ii)=fistrow_golay2(N-ii+2);
%  end
%  fistrow_golay=toeplitz(fistrow_golay1,fistrow_golay2)/N;
 
 fistrow_golay=hadamard(thesize)*diag(a)*dftmtx(thesize)/sqrt(thesize);
 rp8=randperm(N);
 
golay_op=fistrow_golay(rp8(1:M),:);
golay_op=golay_op*sqrt(N/M);
 CDM_op=golay_op;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w=16;          % Generalized Chirp-like sequence
% w_s=0:(1/w):(1-1/w);
% GCL_seq=zeros(1,N);
% zadoff_seq_temp=zadoff(1, N);
% ran_phase=exp(2*w_s*pi*i ); 
% for ooo=1:N;
% GCL_seq(ooo)= zadoff_seq_temp(ooo)*ran_phase(mod(ooo,w)+1);
% end
% 
% fistrow_gcl2=fft(GCL_seq);
% 
% fistrow_gcl1=zeros(1,N);
% fistrow_gcl1(1)=fistrow_gcl2(1);
%  for i=2:N;
%      fistrow_gcl1(i)=fistrow_gcl2(N-i+2);
%  end
%  gcl_full=toeplitz(fistrow_gcl1,fistrow_gcl2)/N;
%  rp4=randperm(N);
%  
% GCL_op=gcl_full(rp4(1:M),:);
% GCL_op=GCL_op*sqrt(N/M);
% CDM_op=GCL_op;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_num = 200;
kvec=4:4: (M/2+M/16);
knum=length(kvec);
success_GaussRD = zeros(knum,1);
success_Gauss = zeros(knum,1);
success_Zad = zeros(knum,1);
success_WHT = zeros(knum,1);
success_CDM = zeros(knum,1);

disp('Now starts the simulations');
for jj = 1:knum,
    k=kvec(jj);    
    
    for i=1:trial_num
        jj
        % create a sparse signal in Psi domain
        alp = [randn(k,1); zeros(N-k,1)];

        p = randperm(N);
        x = alp(p);
        alp=x;
        
        
        % observation
        y0 = Gauss_opRD*x;
        y1 = Gauss_op*x;
        y2 = Zad_op*x;
        y3 = WHT_op*x;
        y4 = CDM_op*x;
        
           % reconstruction using OMP alg.
        %sigma = 0;
        %tau = 0.5;
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
        
        %alp0 = CSRec_SP(k, Gauss_opRD, y0);
        %[alp0,Out0] = fpc_AS(N,Gauss_opRD*idct(eye(N)),y0,mu,[]);
        [alp0,Out0] = fpc_AS(N,Gauss_opRD,y0,mu,[]);
        
   %     alp1 = omp(x, Gauss_op, y1, k, sigma);
        %[alp1,x_debias,objective,times,debias_start,mses]= GPSR_Basic1(y1,Gauss_op,tau);
      %alp1 = CSRec_SP(k, Gauss_op, y1);
        %[alp1,d1]=cosamp(Gauss_op*dct(N),y1,k,tol1);
  %      [alp1 iter_num1] = SAMP(y1, Gauss_op*dct(N), step_size, sigma);
         %alp1 = l1dantzig_pd(Gauss_op'*y1, Gauss_op, [], y1, 4e-3, 5e-3); % l1 minimization
        %alp1 =l1_ls(Gauss_op,y1,lambda,rel_tol);
        %[alp1,Out1] = fpc_AS(N,Gauss_op*idct(eye(N)),y1,mu,[]);
        [alp1,Out1] = fpc_AS(N,Gauss_op,y1,mu,[]);
        
   %     alp2 = omp(x, Zad_op, y2, k, sigma);
        %[alp2,x_debias,objective,times,debias_start,mses]= GPSR_Basic1(y2,Gold_op,tau);
      %alp2 = CSRec_SP(k, Zad_op, y2);
      %[alp2,d2]=cosamp(Zad_op*dct(N),y2,k,tol1);
  %      [alp2 iter_num1] = SAMP(y2, Zad_op*dct(N), step_size, sigma);
        %alp2 = l1dantzig_pd(Zad_op'*y1, Zad_op, [], y2, 4e-3, 5e-3); % l1 minimization
      %alp2 =l1_ls(Zad_op,y2,lambda,rel_tol);
      %[alp2,Out2] = fpc_AS(N,Zad_op*idct(eye(N)),y2,mu,[]);
      [alp2,Out2] = fpc_AS(N,Zad_op,y2,mu,[]);
      
   %     alp3 = omp(x, WHT_op, y3, k, sigma);
        %[alp3,x_debias,objective,times,debias_start,mses]= GPSR_Basic1(y3,WHT_op,tau);
     % alp3 = CSRec_SP(k, WHT_op, y3);
        %[alp3,d3]=cosamp(WHT_op*dct(N),y3,k,tol1);
   %     [alp3 iter_num3] = SAMP(y3, WHT_op*dct(N), step_size, sigma);
        %alp3 = l1dantzig_pd(WHT_op'*y3, WHT_op, [], y3, 4e-3, 5e-3); % l1 minimization
       %alp3 = l1_ls(WHT_op,y3,lambda,rel_tol);
       %[alp3,Out3] = fpc_AS(N,WHT_op*idct(eye(N)),y3,mu,[]);
       [alp3,Out3] = fpc_AS(N,WHT_op,y3,mu,[]);
        
   %     alp4 = omp(x, CDM_op, y4, k, sigma);
        %[alp4,x_debias,objective,times,debias_start,mses]= GPSR_Basic1(y4,CDM_op,tau);
     % alp4 = CSRec_SP(k, CDM_op, y4);
        %[alp4,d4]=cosamp(CDM_op*dct(N),y4,k,tol1);
   %     [alp4 iter_num4] = SAMP(y4, CDM_op*dct(N), step_size, sigma);
         %alp4 = l1dantzig_pd(CDM_op'*y4, CDM_op, [], y4, 4e-3, 5e-3); % l1 minimization
        %alp4 = l1_ls(CDM_op,y4,lambda,rel_tol);
        %[alp4,Out4] = fpc_AS(N,CDM_op*idct(eye(N)),y4,mu,[]);
        [alp4,Out4] = fpc_AS(N,CDM_op,y4,mu,[]);
        
%figure;
%plot(1:N,alp-alp1);
% if snr >50, it is considered as perfect reconstruction

      if snr(alp,alp0)>50
            success_GaussRD(jj) = success_GaussRD(jj)+1;
      end
      
      if snr(alp,alp1)>50
            success_Gauss(jj) = success_Gauss(jj)+1;
      end
        
        
         if snr(alp,real(alp2))>50
            success_Zad(jj) = success_Zad(jj)+1;
         end
        
          if snr(alp,real(alp3))>50
            success_WHT(jj) = success_WHT(jj)+1;
          end
          
           if snr(alp,alp4)>50
             success_CDM(jj) = success_CDM(jj)+1;
           end
    end
    
    success_GaussRD(jj) = success_GaussRD(jj)/trial_num;
    success_Gauss(jj) = success_Gauss(jj)/trial_num;
    success_Zad(jj) = success_Zad(jj)/trial_num;
    success_WHT(jj)=success_WHT(jj)/trial_num;
    success_CDM(jj)=success_CDM(jj)/trial_num;
end 


figure;hold on;
plot([kvec], [success_GaussRD'],'--g');
plot([kvec], [success_Gauss'],'-b+');
plot([kvec], [success_WHT'], '-k*');
plot([kvec], [success_Zad'], '-rs');
plot([kvec], [success_CDM'], '-mo');


grid
%success_Gauss
legend('Complex Random Gaussian','Complex Gaussian Toeplitz','Random binary sequence' ,'FZC squence', 'Golay sequence');
%legend('i.i.d Gaussian Operator','Partial DFT','CDM Matrix','Zadoff-Chu Code');
%legend('i.i.d Gaussian Operator','Partial WHT ','OSTM');
xlabel('No. of Non-zero coefiicients K');
ylabel('Frequency of Exact Reconstruction');

%function snr_val=snr(x0,x1)
%err=norm(x0-x1);
%sig=norm(x0);|
%snr_val=-20*log10(err/sig);

