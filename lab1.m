function lab1()
y = Dirac(10, 20);
figure(1); 
stem(y);  title('Dirac Function');

z =Step(10,20);
figure(2); 
stem(z); title ('Step Function');

e = Ramp(10,20);
figure(3); 
stem(e); title ('Step Function');




end


function y = Dirac(n, N)

if  n >= (N-1)
    disp('Value not in range ');
   
else


arr = zeros(1,N, 'uint32');

arr(n+1) = 1; 
y = arr;

end
end

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