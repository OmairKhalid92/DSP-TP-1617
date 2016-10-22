

function q = Gaussian(N)

y = randn(1,N);
figure;
u = histfit(y);
title('Histogram of Gaussian Distribution')
q=y;
   
end
