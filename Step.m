function y = Step(n, N)

if  n >= (N-1)
    disp('Value not in range ');
    
else


arr = zeros(1,N, 'uint32');


for a = n+1:N 
    arr(a) = 1;   
end
 
y = arr;
end
end