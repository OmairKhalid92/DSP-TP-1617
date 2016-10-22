function y = Ramp(n, N)


t= n+1:1:N;
size(t);


if  n >= (N-1)
    Disp('Value not in range ');
    return
end


arr = zeros(1,N, 'double');


for a = n+1:N 
    arr(a) = 1;
    
end


for a = n+2:N 
    arr(a) = arr(a-1) + arr(a);
    
end
 
y = arr;

end
