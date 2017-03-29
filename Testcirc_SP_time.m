% Test Different Circulant operators

clear all;

N = 512; % signal length
% 32* 512 complex operators; 
M = 32;

thesize=N;   
para= sign(randn(1,thesize));  
fft_mtx=fft(eye(N));  
ifft_mtx=ifft(eye(N));
Rad_mat=ifft_mtx*diag(para)*fft_mtx;

% Rauhut's code, convolution based on Rademacher sequence in fourier domain

% Deterministic Sampling
Rad_fix_comp=Rad_mat(1:M,:);
Rad_fix=[real(Rad_fix_comp);imag(Rad_fix_comp)];

p=randperm(N);
ind=p(1:M);
Rad_rs_comp=Rad_mat(ind,:);
Rad_rs=[real(Rad_rs_comp);imag(Rad_rs_comp)];

%Zaddof_chu 
k=0:N-1;
a=exp(-j*pi*k.^2/N);

p2=randperm(N);
ind2=p2(1:M);
X0=ifft(diag(a));
X1=ifft(X0');
Zad_mat=X1'*N;
Zad_op_comp=Zad_mat(ind2,:);
Zad_op=[real(Zad_op_comp);imag(Zad_op_comp)];

%2M*N Gaussian Operator
Gauss_op=randn(2*M,N);

dct_sparse=0;

if dct_sparse==1,
    psi=dctmtx(N)';
Gauss_op=Gauss_op*psi;
Zad_op=Zad_op*psi;
Rad_rs=Rad_rs*psi;
Rad_fix=Rad_fix*psi;
end

wave_sparse=0;

if wave_sparse==1,
    load wavelet_mat;
    Gauss_op=Gauss_op*Y';
    Gauss_op=norm_mat(Gauss_op);  
    Zad_op=Zad_op*Y';
    Zad_op=norm_mat(Zad_op);
    Rad_rs=Rad_rs*Y';
    Rad_rs=norm_mat(Rad_rs);
    Rad_fix=Rad_fix*Y';
    Rad_fix=norm_mat(Rad_fix);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_num = 300;
kvec=2:2: M+2;
knum=length(kvec);

success_radfix = zeros(knum,1);
success_radrs = zeros(knum,1);
success_zad = zeros(knum,1);
success_Gauss=zeros(knum,1);

disp('Now starts the simulations');
for jj = 1:knum,
    k=kvec(jj);    
    jj
    for i=1:trial_num
        
        % create a sparse signal in Psi domain
        alp = [randn(k,1); zeros(N-k,1)];

        p = randperm(N);
        x = alp(p);
        alp=x;
        
        
        % observation
        y1 = Rad_fix*x;
        y2 = Rad_rs*x;
        y3 = Zad_op*x;
        y4=Gauss_op*x;

        
%            % reconstruction using OMP alg.
%         %sigma = 0;
%         %tau = 0.5;
%         nn = N;
%         tol1=10;
%         tol2=10;
%         tol3=10;
%         tol4=10;
%         sigma = 1e-10;
%         step_size =5;
%         lambda=0.0001;
%         rel_tol=0.005;
%         mu=0.0001;
%         
   %     alp1 = omp(x, Gauss_op, y1, k, sigma);
        %[alp1,x_debias,objective,times,debias_start,mses]= GPSR_Basic1(y1,Gauss_op,tau);
      alp1 = CSRec_SP(k, Rad_fix, y1);
        %[alp1,d1]=cosamp(Gauss_op*eye(N),y1,k,tol1);
  %      [alp1 iter_num1] = SAMP(y1, Gauss_op*eye(N), step_size, sigma);
         %alp1 = l1dantzig_pd(Gauss_op'*y1, Gauss_op, [], y1, 4e-3, 5e-3); % l1 minimization
        %alp1 =l1_ls(Gauss_op,y1,lambda,rel_tol);
        %[alp1,Out1] = fpc_AS(N,Gauss_op,y1,mu,[]);
        
   %     alp2 = omp(x, Zad_op, y2, k, sigma);
        %[alp2,x_debias,objective,times,debias_start,mses]= GPSR_Basic1(y2,Gold_op,tau);
      alp2 = CSRec_SP(k, Rad_rs, y2);
      %[alp2,d2]=cosamp(Zad_op*eye(N),y2,k,tol1);
  %      [alp2 iter_num1] = SAMP(y2, Zad_op*eye(N), step_size, sigma);
        %alp2 = l1dantzig_pd(Zad_op'*y1, Zad_op, [], y2, 4e-3, 5e-3); % l1 minimization
      %alp2 =l1_ls(Zad_op,y2,lambda,rel_tol);
      %[alp2,Out2] = fpc_AS(N,Zad_op,y2,mu,[]);
      
   %     alp3 = omp(x, WHT_op, y3, k, sigma);
        %[alp3,x_debias,objective,times,debias_start,mses]= GPSR_Basic1(y3,WHT_op,tau);
      alp3 = CSRec_SP(k, Zad_op, y3);
        %[alp3,d3]=cosamp(WHT_op*eye(N),y3,k,tol1);
   %     [alp3 iter_num3] = SAMP(y3, WHT_op*eye(N), step_size, sigma);
        %alp3 = l1dantzig_pd(WHT_op'*y3, WHT_op, [], y3, 4e-3, 5e-3); % l1 minimization
       %alp3 = l1_ls(WHT_op,y3,lambda,rel_tol);
       %[alp3,Out3] = fpc_AS(N,WHT_op,y3,mu,[]);
        
   %     alp4 = omp(x, CDM_op, y4, k, sigma);
        %[alp4,x_debias,objective,times,debias_start,mses]= GPSR_Basic1(y4,CDM_op,tau);
     % alp4 = CSRec_SP(k, CDM_op*eye(N), y4);
        %[alp4,d4]=cosamp(CDM_op*eye(N),y4,k,tol1);
   %     [alp4 iter_num4] = SAMP(y4, CDM_op*eye(N), step_size, sigma);
         %alp4 = l1dantzig_pd(CDM_op'*y4, CDM_op, [], y4, 4e-3, 5e-3); % l1 minimization
        %alp4 = l1_ls(CDM_op,y4,lambda,rel_tol);
      %  [alp4,Out4] = fpc_AS(N,CDM_op*eye(N),y4,mu,[]);
         alp4 = CSRec_SP(k, Gauss_op, y4);
%figure;
%plot(1:N,alp-alp1);
% if snr >50, it is considered as perfect reconstruction
      if snr(alp,alp1)>50
            success_radfix(jj) = success_radfix(jj)+1;
      end
        
        
         if snr(alp,real(alp2))>50
            success_radrs(jj) = success_radrs(jj)+1;
         end
        
          if snr(alp,real(alp3))>50
            success_zad(jj) =  success_zad(jj)+1;
          end
          
          if snr(alp,alp4)>50
            success_Gauss(jj) =  success_Gauss(jj)+1;
          end
           
    end
    
    success_radfix(jj) = success_radfix(jj)/trial_num;
    success_radrs(jj) = success_radrs(jj)/trial_num;
    success_zad(jj)=success_zad(jj)/trial_num;
    success_Gauss(jj)=success_Gauss(jj)/trial_num;
end 


figure;hold on;
plot([kvec], [success_Gauss'], '-go');
plot([kvec], [success_radfix'],'-b+');
plot([kvec], [success_radrs'], '-k*');
plot([kvec], [success_zad'], '-rs');


grid
%success_Gauss
legend('Gaussian','Rad-fix','Rad-rs' ,'Chirp-rs');
xlabel('No. of Non-zero coefiicients K');
ylabel('Frequency of Exact Reconstruction');
title('Result for sparse signal in the time domain');
%function snr_val=snr(x0,x1)
%err=norm(x0-x1);
%sig=norm(x0);|
%snr_val=-20*log10(err/sig);

