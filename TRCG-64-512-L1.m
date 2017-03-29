figure;hold on;
success_Gauss=[1,0.995,0.98,0.95,0.91,0.85,0.65,0.455,0.27,0.103,0.042,0.015,0.005,0.001,0,0];
success_WHT=[1,0.993,0.975,0.93,0.902,0.818,0.623,0.433,0.218,0.098,0.038,0.014,0.003,0,0,0];
success_Zad=[1,0.995,0.965,0.92,0.895,0.794,0.598,0.415,0.212,0.092,0.038,0.006,0,0,0,0];
success_CDM=[1,0.994,0.971,0.918,0.899,0.804,0.603,0.428,0.221,0.091,0.0391,0.012,0.001,0,0,0];


plot([2:2:32], [success_Gauss'],'-b+');
plot([2:2:32], [success_WHT'], '-k*');
plot([2:2:32], [success_Zad'], '-rs');
plot([2:2:32], [success_CDM'], '-mo');


grid
%success_Gauss
legend('Complex Gaussian Toeplitz','Rademacher sequence' ,'Zadoff-Chu Code', 'Golay sequence');
%legend('i.i.d Gaussian Operator','Partial DFT','CDM Matrix','Zadoff-Chu Code');
%legend('i.i.d Gaussian Operator','Partial WHT ','OSTM');
xlabel('No. of Non-zero coefiicients K');
ylabel('Frequency of Exact Reconstruction');