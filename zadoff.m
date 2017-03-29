function seq = zadoff(root,N)    
for n=0:(N-1)
    seq(n+1)=exp(-j*(pi*root*n*(n+1))/N);
end
