

function q = Box(hw, n, N)

disp(hw);
disp(n);
disp(N);
disp((n-hw));

if (hw+n) > (N)
    disp(' Upper boundary greater than Range ');
    
elseif (n-hw) < 1
    disp(' Lower boundary greater than Range ');

else
    
arr = zeros(1,N, 'uint32');


for f = n-hw:n+hw 
    arr(f) = 1;
    
end

q = arr;


   
end
end
