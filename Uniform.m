function q = Uniform(N)

y = rand(1,N);
figure;
u = histogram (y);
title('Histogram of Uniform Distribution')
q=y;
   
end
