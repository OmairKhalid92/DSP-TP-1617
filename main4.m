%1.1
f = 5; fs = 50;
t = 0: 1/fs : 1;
figure(1);
title('Sin Wave');
xn = sin(2*pi*f*t);
N = length(xn);
fr = (-N/2 : N/2-1)*fs/N;
xf = fftshift(fft(xn));
subplot(221); plot(t, xn); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');
 
%1.2
figure(2);
title('Cos Wave');
xc = cos(2*pi*f*t);
t = 0: 1/fs : 1;
N = length(xc);
fr = (-N/2 : N/2-1)*fs/N;
xd = fftshift(fft(xc));
subplot(221); plot(t, xc); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xd)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xd)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xd)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');
 
%1.3
 
t = 0: 1/fs : 1;

xw = square(2*pi*f*t,50);

N = length(xw);
fr = (-N/2 : N/2-1)*fs/N;
xe = fftshift(fft(xw));
figure(3);

subplot(221); plot(t, xw); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xe)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xe)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xe)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');


%1.4

noise = rand(1,10000);
tt = 1:1:10000
N = length(noise);
fr = (-N/2 : N/2-1)*fs/N;
xf = fftshift(fft(noise));
figure(4);

subplot(221); plot(tt, noise); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');
 

%%2
 
%Original Signal
 
n = 0:1/1000000:1;
x = 3*cos(2*pi*f1*n) + 4*sin(2*pi*f2*n);
figure; plot(n,x);
suptitle('Original Signal');
 
%2.1 : fs  = 10
 
f1 = 5;
f2 = 20;
fs = 10;
n = 0:1/fs:1;
length(n)

x = 3*cos(2*pi*f1*n)+ 4*sin(2*pi*f2*n);

figure; plot(n,x);
subplot(221); plot(n, x); title('F = 10 Sampled Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
 
% % fs  = 20

fs = 20;
n = 0:1/fs:1;
length(n)

x = 3*cos(2*pi*f1*n)+ 4*sin(2*pi*f2*n);
subplot(222); plot(n, x); title('F = 20 Sampled Signal'); xlabel('Time(sec)'); ylabel('Amplitude');

% % fs  = 25
fs = 25;
n = 0:1/fs:1;
length(n)

x = 3*cos(2*pi*f1*n)+ 4*sin(2*pi*f2*n);
subplot(223); plot(n, x); title('F = 25 Sampled Signal'); xlabel('Time(sec)'); ylabel('Amplitude');

% % fs  = 40

fs = 40;
n = 0:1/fs:1;
 
x = 3*cos(2*pi*f1*n)+ 4*sin(2*pi*f2*n);
 subplot(224); plot(n,x); title('Fs = 40'); xlabel('Time'); ylabel('Amplitude');


% % fs  = 50

fs = 50;
n = 0:1/fs:1;

x = 3*cos(2*pi*f1*n)+ 4*sin(2*pi*f2*n);
figure;subplot(321); plot(n, x); title('Fs = 50'); xlabel('Time'); ylabel('Amplitude');

% % fs  = 100
fs = 100;
n = 0:1/fs:1;
x = 3*cos(2*pi*f1*n)+ 4*sin(2*pi*f2*n);
subplot(322); plot(n, x); title('Fs = 100'); xlabel('Time'); ylabel('Amplitude');


% % fs  = 150

fs = 150;
n = 0:1/fs:1;
 
x = 3*cos(2*pi*f1*n)+ 4*sin(2*pi*f2*n);
subplot(323); plot(n, x); title('Fs = 150'); xlabel('Time'); ylabel('Amplitude');
 

Discuss Aliasing Effects in Time Domain

fs = [10 20 25 40 50 100 150];
f1 = 5;
f2 = 20;
for k = 1:length(fs)
    
fsamp=fs(1,k);
n = 0:1/fsamp:1;
x = 3*cos(2*pi*f1*n)+ 4*sin(2*pi*f2*n);

N = 1000;
xq = fftshift(fft(x,N));
fr = (-N/2 : N/2-1)*fsamp/N;
figure;
subplot(221); plot(fr, abs(xq)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(222); plot(fr, real(xq)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(223); plot(fr, imag(xq)); title('Imaginary'); xlabel('Frequency'); ylabel('Re(X(f))');
suptitle( ['Fsamp = ' num2str(fsamp) ] );

end


%3
% %[COULD NOT COMPLETE IT. UNDER CONSTRUCTION]
% %3.1
% 
% srcFiles = dir('C:\Users\Omair Khalid\Downloads\DSP-TP-1617 (1)\DSP-TP-1617\Lab4\images\1D-DFT\*.tif');  % the folder in which ur images exists
% N = 54;
% my_im_cell = cell(N,1);
% I_norm = cell(N,1);
% I_norm_r = cell(N,1);
% I_norm_r_p = cell(N,1);
% I_norm_r_p_fft = cell(N,1);
% 
% for i = 1 : length(srcFiles)
%     filename = strcat('C:\Users\Omair Khalid\Downloads\DSP-TP-1617 (1)\DSP-TP-1617\Lab4\images\1D-DFT\',srcFiles(i).name);
%     I = imread(filename);
%     size(I);
%     r = I(:,:,1);
% 
%     my_im_cell{i} = r; % where k is the image index/number
%     
%     I_norm{i} = (r - min(min(r)))./( max( max(r) ) - min( min(r) ) );
%     
%     min_size1 = size(I_norm(1),2);
%     min_size = size(I_norm(i),2);
%     
%     if min_size <= min_size1
%          min_size1 = min_size;
%     end
%     
%    
% end
% 
% %-----------------------
% for i = 1 : length(srcFiles)
% 
%     
%     I_norm_r{i} = imresize(I_norm{i},[min_size1 (min_size1/(size(I_norm_r{i},1)/size(I_norm_r{i},2)));
%     
%     
%     I_norm_r_ptemp = I_norm_r{i}( 50 , :);
%     
%     %imshow(I_norm_r_ptemp,[]);   
%     I_norm_r_p{i} = I_norm_r_ptemp;  
%     
%    
%     
%     NN = min_size1;
%     I_norm_r_p_fft{i} = fftshift(fft(I_norm_r_p{i}));
%     fr = (-NN/2 : NN/2-1);
%    
% %subplot(221); plot(tt, I_norm_r_p_fft{1}); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
%    
% end
% 
% 
% %-----------------------
% 
% figure;
% stem(fr, abs(I_norm_r_p_fft{1})); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');

    
%Sort and Comment

image_folder='C:\Users\Tajwar Abrar Aleef\Documents\GitHub\DSP-TP-1617\Lab4\images\1D-DFT\'; 
cd (image_folder)
files=dir;
amount_of_images=length(files)-2
 for i = 1:amount_of_images
     cell{i} = imread(files(i+2).name);
 end
 
[y,x] = size(cell{1});
 for idx = 2:amount_of_images
     [a,b]=size(cell{idx})
     if a*b<x*y
     y=a;
     x=b;
     end
 end
 
 p=1;q=1;
 for i=1:amount_of_images
    I = imread(files(i+2).name);
    I = double(I);
    [a,b,c]= size(I);
        if c==4 % some of the images are four dimensional and I took only one dimentsion
            I=I(:,:,1);
        end
        I=I/(max(I(:)));  
        I_resized = imresize(I, [y x]);
        z=round(y/2);
        Profile_1D =I_resized(z,:);
   
    % 1-D DFT of the profile
    N=x;
    fr = (-N/2 : N/2-1);
    x1 = ifftshift(fft(Profile_1D));
%      figure;
%     stem(image_DFT)
    % 1-D DFT of the profile another row. if the image is barcode image the
    % DFT of the the second profile will be almost same as the dft of the
    % second profile which taken far from the center
    Profile_1D2 = I_resized(z+2,:);
    Profile_1D3 = I_resized(z-2,:);
    x2 = ifftshift(fft(Profile_1D2));
    x3 = ifftshift(fft(Profile_1D3));
%    plot(fr, abs(image_DFT)); title('Magnitude'); ...
%    xlabel('Frequency'); ylabel('|image_DFT|');

% Discussion
% For fourier transform of barcode images must have a SIGNIFICANT value when f different from zero, while
%  barcode images have DFT which concentrated at around f=0


%==============================================

% Therefore to determine wether agiven image is bar code or not we can
% see the value of the FFT of the 1-D profile of the images at high frequency(taking Threshold) or use a HIGH
% PASS FILTER. Barcodes Have high frequency components
           v1=abs(x1);
           v2=abs(x2);
           v3=abs(x3);
           threshold1=max(v1);
           threshold2=abs(2*max(v1)-max(v2)-max(v3));
           threshold(i) = threshold1*threshold2;
           
           if (threshold(i)<48)
                T(p)=i ; 
                disp('Image is Barcode')   
                p=p+1;
           else
                disp('not barcode')
                O(q)=i;
                q=q+1;
           end
 end
 figure;
 stem(threshold);
 disp('Barcode');
 disp(T);
 disp('non barcode');
 disp(O);
 Barcode=[1,2,6,44:54];
 NonBarcode = [3,4,5,7:43];
 
 x=length(Barcode);
 u=0;v=0;
 for i=1:length(Barcode)
   for j=1:length(T)
     if T(j)==Barcode(i)
         u=u+1;
     end
   end
 end
 
 for i=1:length(NonBarcode)
   for j=1:length(T)
     if T(j)==NonBarcode(i)
         v=v+1;
     end
   end
 end
 
 v=v+x-u;
 percentage_of_accuracy= ((amount_of_images-v)/amount_of_images)*100;
 disp('')
 sprintf('Percentage of correct distinction between the Barcode and NonBarcode images is: %d percent', round(percentage_of_accuracy))



 
 for i = 1 : length(srcFiles)
 
    I_norm_r{i} = imresize(I_norm{i},[min_size1 (min_size1/(size(I_norm_r{i},1)/size(I_norm_r{i},2)));

 end



