% 1.1 - Chebyshev Filter

[b,a] = cheby1(3,1,0.5,'low');
[H, W]= freqz(b,a);

subplot(2,2,1);plot(W/pi, abs(H)); title('Low Pass - Chebyshev');


[b,a] = cheby1(3,1, 0.5, 'high');
[H, W]= freqz(b,a);
subplot(2,2,2);plot(W/pi, abs(H)); title('High Pass - Chebyshev');

[b,a] = cheby1(3,1, [0.2 0.4]);
[H, W]= freqz(b,a);
subplot(2,2,3);plot(W/pi, abs(H));title('Band-Pass - Chebyshev');

[b,a] = cheby1(3,1, [0.2 0.4] , 'stop');
[H, W]= freqz(b,a);
subplot(2,2,4); plot(W/pi, abs(H)); title('Band-Stop - Chebyshev');



%1.2 - Butterworth Filter
figure;

[b,a] = butter(3,0.5,'low');
[H, W]= freqz(b,a);

subplot(2,2,1);plot(W/pi, abs(H)); title('Low Pass - Butterworth');


[b,a] = butter(3,0.5, 'high');
[H, W]= freqz(b,a);
subplot(2,2,2);plot(W/pi, abs(H)); title('High Pass - Butterworth');

[b,a] = butter(3,[0.2 0.4]);
[H, W]= freqz(b,a);
subplot(2,2,3);plot(W/pi, abs(H));title('Band-Pass - Butterworth');

[b,a] = butter(3,[0.2 0.4] , 'stop');
[H, W]= freqz(b,a);
subplot(2,2,4); plot(W/pi, abs(H)); title('Band-Stop - Butterworth');

%1.3 - Low Pass Cheyshev with varying order
figure;
[b,a] = cheby1(3,1,0.5,'low');
[H, W]= freqz(b,a);

subplot(2,2,1);plot(W/pi, abs(H)); title('Low Pass Chebyshev filter - Order:3');

[b,a] = cheby1(5,1,0.5,'low');
[H, W]= freqz(b,a);

subplot(2,2,2);plot(W/pi, abs(H)); title('Low Pass Chebyshev filter - Order:5');

[b,a] = cheby1(10,1,0.5,'low');
[H, W]= freqz(b,a);

subplot(2,2,3);plot(W/pi, abs(H)); title('Low Pass Chebyshev filter - Order:10');

[b,a] = cheby1(20,1,0.5,'low');
[H, W]= freqz(b,a);

subplot(2,2,4);plot(W/pi, abs(H)); title('Low Pass Chebyshev filter - Order:20');

% As the order increases, the number of ripples increase and The transition
% band gets steeper (close to vertical drop on n=20)


%Question 2.1
dirac = zeros(1,40);
dirac(20) = 1;

%Setting value of scaling
sc = 0.5; 
Ts = 1; 
a = sc*Ts;  
a_signal = exp(-a) ;

%Anti-Causal part of smoothing filter
dirac_aCausal_smooth = zeros(length(dirac),1);
dirac_length=length(dirac)-2:-1:1;

for i =  dirac_length 
 dirac_aCausal_smooth(i) = sc*a*a_signal*dirac(i+1)+(2*a_signal)*dirac_aCausal_smooth(i+1)-(a_signal^2)*dirac_aCausal_smooth(i+2);
end
figure()
stem (dirac_aCausal_smooth) ;
title('Anti-causal Smoothing'); 

%Causal part of smoothing
dirac_causal_smooth = zeros(length(dirac),1);
for i = 3:length(dirac);
 dirac_causal_smooth(i) = -sc*a*a_signal*dirac(i-1)+(2*a_signal)*dirac_aCausal_smooth(i-1)-(a_signal^2)*dirac_aCausal_smooth(i-2);
end    
figure()
stem (dirac_causal_smooth) ;
title('Causal Smmothing'); 

%st input Signal
st10=step(40,10);
st30=step(40,30);
st=st10-st30;
figure()
stem(st)

%Causal Derivative function
st_causal = zeros(length (st),1) ;

for i = 3 : length(st)
 st_causal(i) = st(i)+a_signal*(a-1)*st(i-1)+(2*a_signal)*st_causal(i-1)-(a_signal^2)*st_causal(i-2) ;
end
figure()
stem (st_causal) ;
title('Causal Deravative');

%Anti-causal derivative function
st_aCausal = zeros(length (st),1);
st_length = length(st)-2 : -1 : 1 ; 

for i = st_length  
 st_aCausal(i) = a_signal*(a+1)*st(i+1)-(a_signal^2)*st(i+2)+(2*a_signal)*st_aCausal(i+1)-(a_signal^2)*st_aCausal(i+2) ;  
end
figure()
stem (st_aCausal);
title('Anti-Causal Deravative');

%Question 3
barb = imread('E:\VIBOT\DSP\DSP LAB 1\DSP-TP-1617\Lab6-7\images\barbara.gif');
figure();
imshow(barb);  
case1 = zeros(size(barb)); 
case2 = zeros(size(barb)); 


for i = 1:size(barb, 2)
    image1 = barb(:,i);
    
    resp_causal = zeros(length (image1),1) ;
    for i = 3 : length(image1)
     resp_causal(i) = image1(i)+a_signal*(a-1)*image1(i - 1)+(2*a_signal)*resp_causal(i-1)-(a_signal^2)*resp_causal(i-2) ;
    end
    
    resp_aCausal = zeros(length (image1 ),1) ;
    barb_length = length(image1)-2 : -1 : 1 ;
    for i =  barb_length
     resp_aCausal(i) = a_signal*(a+1)*image1(i+1)-(a_signal^2)*image1(i+2)+(2*a_signal)*resp_aCausal(i+1)-(a_signal^2)*resp_aCausal(i+2) ;
    end
    
    resp = resp_causal + resp_aCausal;
    
    case1(:,i) = resp;            
end
figure();
imshow (case1, []);

for i = 1:size(barb, 2)

    image2 = barb(i,:);
    
    resp_causal = zeros(length (image2),1) ;
    for i = 3 : length(image2)
     resp_causal(i) = image2(i)+a_signal*(a-1)*image2(i - 1)+(2*a_signal)*resp_causal(i-1)-(a_signal^2)*resp_causal(i-2) ;
    end
    
    resp_aCausal = zeros(length (image2 ),1) ;
    barb_length = length(image2)-2 : -1 : 1 ;
    for i =  barb_length
     resp_aCausal(i) = a_signal*(a+1)*image2(i+1)-(a_signal^2)*image2(i+2)+(2*a_signal)*resp_aCausal(i+1)-(a_signal^2)*resp_aCausal(i+2) ;
    end
    resp2 = resp_causal + resp_aCausal;
    
    case2(i,:) = resp2;   
end
figure();
imshow (case2, []);













