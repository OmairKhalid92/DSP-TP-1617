function r = convFn(x, y)




lx = length(x);
ly = length(y);

if length(x)>=length(y);
    
    %ld = length(x)-length(y);
    %y = [y zeros(1,ld)];

    
    yf = fliplr(y);
    %yf = [xf zeros(ly) zeros(lx)]
    zy =zeros(1,length(y)-1);
    xe = cat(2, zy ,x, zy);
    
    yf = double(yf);
    %cr = xe.*(double(yf));
    
    
yp = yf;

cr = zeros(1,length(xe));

    for i=1:1:(lx+ly-1)
        %pad, multiply, shift and add
        
        %SLIDE & PAD
        
        yp = yf;
        yp1 = cat(2,zeros(1,i-1),yp);
        yp2 = cat(2,yp1,zeros(1,length(xe)-length(yp1))); 
        
        
        
        yp = double(yp2);
        
        
        %MULTIPLY
        
        cv = xe.*yp;
        
        %ADD
        
        cr(1,i) = sum(cv);
          
    end

end
% cr(1,1:length(cr)-1) = cr(1,2:length(cr));
% cr(1,length(cr))= 0;

r = cr;

end