
function y = Geo(a, n, N)

if  n >= (N-1)
    disp('Value not in range ');
    
else


arr = zeros(1,N, 'uint32');


for f = n+1:N 
    arr(f) = 1;
    
end

k=0;


for e = n+1:N 

    arr(e) = (a^k).*arr(e);

    k=k+1;
    
end
 
y = arr;
 

end
end