% Histogram matching
% Lav standard normalfordelte data ud fra uniformt fordelte data.
N = 1000;
x = rand(1,N);          % N datapunkter, som er uniformt fordelte på intervallet [0,1]
Fx = unifcdf(x,0,1);    % Værdien af fordelingsfunktionen (CDF) for hvert x
                        % Bemærk, at for en uniform fordeling på
                        % intervallet [0,1], så er Fx(x) = x;
yrange = -5:0.1:5;      % Range af y-værdier, som vi skal bruge til at beregne den inverse af Fy.
Fy = normcdf(yrange,0,1); % En slags lookup table
% Invers mapping:
%   For hvert x, beregn Fx(x) og søg efter den værdi af Fy(y), som ligger
%   tættest på Fx(x).
for i = 1:N
    distance = abs(Fx(i)-Fy);
    [min_val,min_ix] = min(distance);
    y(i) = yrange(min_ix);
end
% Verificer, at de generede y-værdier følger en standard normalfordeling
step = 0.5;
yrange = -5:step:5;
h = hist(y,yrange);
fy = normpdf(yrange,0,1);
figure
bar(yrange,h/N/step)
hold on
plot(yrange,fy,'r')
hold off
legend('Data histogram','Teoretisk PDF')

% Histogram udligning i et billede
I = imread('camera.png');
x = 0:255;
num_pixels = prod(size(I));
fx = hist(I(:),x)/num_pixels;
Fx = cumsum(fx);
yrange = 0:255;
I_histeq = zeros(size(I),'uint8');
for i = 1:length(x)
    ix = find(I==x(i));
    distance = abs(Fx(i)-Fy);
    [min_val,min_ix] = min(distance);
    y = yrange(min_ix);
    I_histeq(ix) = y;
end
fy = hist(I_histeq(:),yrange)/num_pixels;
Fy = cumsum(fy);
figure,
subplot(2,3,1),plot(x,fx),title('f_x')
subplot(2,3,2),plot(x,Fx),title('F_x')
subplot(2,3,3),imshow(I)
subplot(2,3,4),plot(yrange,fy),title('f_y')
subplot(2,3,5),plot(yrange,Fy),title('F_y')
subplot(2,3,6),imshow(I_histeq)