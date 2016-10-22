

function [t,q] = Sinfn(f,n,fs)

   
%sin(2pifnTs)

t = 0:1/fs:n/f; 
    arr = sin(2*pi*f*t);
   % disp(arr(w));


q = arr;


   
end