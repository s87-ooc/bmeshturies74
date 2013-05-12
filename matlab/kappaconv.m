

spectral = [];

for k=0:0.05:10
    I = A - k*M + B; 
    D = diag(diag(I)); 
    spectral = [spectral; k max(abs(eigs(inv(D)*(I-D))))];
end

plot(spectral(:,1), spectral(:,2));
hold on;
xlabel('kappa');
ylabel('rayon spectral');