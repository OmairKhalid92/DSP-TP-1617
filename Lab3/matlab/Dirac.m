function y = Dirac(n, N)

if  n >= (N-1)
    disp('Value not in range ');
   
else


arr = zeros(1,N, 'uint32');

arr(n+1) = 1; 

y = arr;

end
end