function r = conv2dFn (fil, img)

 [e,t] = size(img);
 %e = no. of traversals
 %t = no. of weighted sum operations per traversal
 
 res_img = zeros(e,t);

fil_f = fliplr(fil);
pad = floor(size(fil)/2);
p = pad(1);
p_img= padarray(img, [p p], 0);
p_img = double(p_img);


start_point_col = p+1;
end_point_col = p+t;

start_point_row = p+1;
end_point_row = p+e;

 a=1;b=1;
dres = zeros(e,t);
for i=start_point_row:1:end_point_row

   
    for j=start_point_col:end_point_col
        dres = p_img(i-p:i+p,j-p:j+p).*fil_f;
        S = sum(dres);
        S = sum(S);
        res_img(a,b) = S;
  
        b=b+1;
    end
    a=a+1;
    b=1;
end

r = res_img;
end