%%OMAIR KHALID - OK

%%1.1

%ConvFn is implemented in a .m file

%%1.2

% %Defining Functions
x = [1 2 3 4];
d = double(Dirac(0,5));
s = double(Step(1,5));

st = 1:1:5;
e = exp(st/2);

sig = [-1 1];
% 
% %Applying convFn function and displaying the results
% 
dr = convFn(d,x);
 
figure();
stem(dr)
title('Convolution Result: x & Dirac');

sr = convFn(s,x);

figure();
stem(sr)
title('Convolution Result: x & Step');


er = convFn(e,x);

figure();
stem(er)
title('Convolution Result: x & Exp');

sigr = convFn(x,sig);

figure;
stem(sigr);
title('Convolution Result: x & Sig = [-1 1]');

%%1.3

%Convolution with Symmetrically Padded Input

x = [ 1 2 3 4 ];
lx=length(x);


s = double(Step(1,5));

ls=length(s);

spx = padarray(x,ls-1,'symmetric','both');


anss = zeros(1,lx+ls-1);

for i=1:1:(lx+ls-1)
        
        
        
         anss(i) = sum(spx(i:i+4).*s);
        
end
figure;
stem(anss);
title('Convolution with Symmetrically Padded Input');

%Convolution with Periodically Padded Input

x = [ 1 2 3 4 ];
lx=length(x);


s = double(Step(1,5));

ls=length(s);

spx = padarray(x,ls-1,'circular','both');


anss = zeros(1,lx+ls-1);

for i=1:1:(lx+ls-1)
        
        
        
         anss(i) = sum(spx(i:i+4).*s);
        
end
figure;
stem(anss);
title('Convolution with Periodically Padded Input');  


%Convolution with Consant-Padded Input

x = [ 1 2 3 4 ];
lx=length(x);


s = double(Step(1,5));

ls=length(s);

spx = padarray(x,[0 ls-1],1,'both');


anss = zeros(1,lx+ls-1);

for i=1:1:(lx+ls-1)
        
        
        
         anss(i) = sum(spx(i:i+4).*s);
        
end
figure;
stem(anss);
title('Convolution with Constant-Padded Input'); 

%%2.1

%Implemented in another .m file

%2.2

K = [1 4 6 4 1; 4 16 24 16 4; 6 24 36 24 6; 4 16 24 16 4; 1 4 6 4 1]./256;


lena = imread('E:\VIBOT\DSP\DSP LAB 1\DSP-TP-1617\Lab3\images\lena-grey.bmp');
imshow(lena);
img = double(lena)./256;

r = conv2dFn(K,img);

imshow(lena);
figure();imshow(r);


%2.3

sob_x = [-1 0 1;-2 0 2;-1 0 1];

rx = conv2dFn(sob_x,img);

figure;imshow(rx);title('Sobel X filter on Lena');

sob_y = [1 2 1;0 0 0;-1 -2 -1];

ry = conv2dFn(sob_y,img);

figure;imshow(ry);title('Sobel Y filter on Lena');


%%3.1, 3.2 & 3.3

close all;
clc;
clear all;


text=imread('E:\VIBOT\DSP\DSP LAB 1\DSP-TP-1617\Lab3\images\text.png');
a=imread('E:\VIBOT\DSP\DSP LAB 1\DSP-TP-1617\Lab3\images\a.png');
a= imcomplement(a); %converting the black font to white and black background to white  

%applying OTSU thresholding 
otsu_text = graythresh(text);
otsu_a = graythresh(a);

%Converting the OSTU threshold image to binary for both Text and A and
%displaying it
binary_a = im2bw(a,otsu_a);
figure('Name','a.png binarized using Otsu tresholding','NumberTitle','off');
imshow(binary_a);

binary_text = im2bw(text,otsu_text);
figure('Name','a.png binarized using Otsu tresholding','NumberTitle','off');
imshow(binary_text)

%Applying cross correlationg for both the binary images and showing it
correlated_image= xcorr2(im2double(binary_text),im2double(binary_a));
figure('Name','Correlation of binarized images of text and a','NumberTitle','off');
imshow(correlated_image,[]);

[rowMax, columnMax]=find(correlated_image == max(max(correlated_image))); %finding the indices of the maximum values of the crosscorellated image
figure('Name','Detected letter A encircled in red','NumberTitle','off');
imshow(binary_text, []); 
i=find(rowMax==(max(max(rowMax)))) %finding the indices of the 11th a
rowMax(i)=[];
columnMax(i)=[]; %removing the 11th a to show only first 10 a
hold on;
plot(columnMax-10,rowMax-10, 'ro','MarkerSize',15), %plotting red circles over the detected A.





