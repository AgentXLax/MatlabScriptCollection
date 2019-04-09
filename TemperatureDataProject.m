myFile = fopen('data.txt'); 
tLine = fgetl(myFile);
while tLine(1) == '%'
    tLine = fgetl(myFile);
end
data = textscan(myFile,'%d %f %f %f %f %f %f %f %f'); 

% reads the file column by column. Each column is a cell array
% Year, Annual Anomaly, Annual Unc., Five-year Anomaly, Five-year Unc., 
% Annual Anomaly, Annual Unc., Five-year Anomaly, Five-year Unc.
% 
% Estimated Jan 1951-Dec 1980 global mean temperature (C)
%   Using air temperature above sea ice:   14.720 +/- 0.046
%   Using water temperature below sea ice: 15.304 +/- 0.046
% this point is the reference point

fclose(myFile);

syms a1 b1 a2 b2 c2 a3 b3 c3 d3
avecLin = [a1 b1];
avecQuad = [a2 b2 c2];
avecCub = [a3 b3 c3 d3];


p1 = @(a,xdata) a(1) + a(2)*(xdata - 1850);
p2 = @(a,xdata) a(1) + a(2)*(xdata - 1850) + a(3)*(xdata - 1850).^2;
p3 = @(a,xdata) a(1) + a(2)*(xdata - 1850) + a(3)*(xdata - 1850).^2 + ...
    a(4)*(xdata - 1850).^3;
pe = @(a,xdata) a(1) + a(2)*exp(a(3)*(xdata - 1850));
pmixLin = @(a,xdata) (a(1) + a(2)*(xdata - 1850))...
    .*exp(a(3)*(xdata - 1850));
pmixQuad = @(a,xdata) (a(1) + a(2)*(xdata - 1850) + ...
    a(3)*(xdata - 1850).^2).*exp(a(4)*(xdata - 1850));

tdata = double(data{1});
ydata = data{2};
%%

difLin(1) = sum(diff((p1(avecLin,tdata) - ydata).^2,avecLin(1))) == 0;
difLin(2) = sum(diff((p1(avecLin,tdata) - ydata).^2,avecLin(2))) == 0;
[A, b] = equationsToMatrix(difLin, avecLin)
avecLin = double(A\b)
p1 = @(x) avecLin(1) + avecLin(2)*(x - 1850);
 
f1 = figure;
scatter(tdata, ydata, 10,[249/256 110/256 59/256], 'filled')
hold on;
fplot(p1, [1850 2100])
title('Linear Fit: p_1(t)')
xlabel('Year')
ylabel('Temperature Difference (C)')
hold off;
print -depsc LinFit.eps
est2100Lin = p1(2100)

resErrorLin = abs(sum((p1(double(tdata)) - ydata).^2))
%%
diffQuad(1) = sum(diff((p2(avecQuad,tdata) - ydata).^2,avecQuad(1))) == 0;
diffQuad(2) = sum(diff((p2(avecQuad,tdata) - ydata).^2,avecQuad(2))) == 0;
diffQuad(3) = sum(diff((p2(avecQuad,tdata) - ydata).^2,avecQuad(3))) == 0;

[A, b] = equationsToMatrix([diffQuad(1) diffQuad(2) diffQuad(3)], avecQuad)
avecQuad = double(A\b)
p2 = @(x) avecQuad(1) + avecQuad(2)*(x-1850) + avecQuad(3)*(x-1850).^2;

f2 = figure;
scatter(tdata, ydata, 10,[249/256 110/256 59/256], 'filled')
hold on;
fplot(p2, [1850 2100])
title('Quadratic Fit: p_2(t)')
xlabel('Year')
ylabel('Temperature Difference (C)')
hold off;
print -depsc QuadFit.eps

est2100Quad = p2(2100)

resErrorQuad = abs(sum((p2(double(tdata)) - ydata).^2))
%%
diffCub(1) = sum(diff((p3(avecCub,tdata) - ydata).^2,avecCub(1))) == 0;
diffCub(2) = sum(diff((p3(avecCub,tdata) - ydata).^2,avecCub(2))) == 0;
diffCub(3) = sum(diff((p3(avecCub,tdata) - ydata).^2,avecCub(3))) == 0;
diffCub(4) = sum(diff((p3(avecCub,tdata) - ydata).^2,avecCub(4))) == 0;

[A, b] = equationsToMatrix(diffCub, avecCub)
avecCub = double(A\b)
p3 = @(x) avecCub(1) + avecCub(2)*(x - 1850) + ...
    avecCub(3)*(x - 1850).^2 + avecCub(4)*(x - 1850).^3;

f3 = figure;
scatter(tdata, ydata, 10,[249/256 110/256 59/256], 'filled')
hold on;
fplot(p3, [1850 2100])
title('3rd Degree Polynomial Fit: p_3(t)')
xlabel('Year')
ylabel('Temperature Difference (C)')
hold off;
print -depsc CubFit.eps

est2100Cub = p3(2100)

resErrorCub = abs(sum((p3(double(tdata)) - ydata).^2))

%%
[avecExp resErrorExp] = lsqcurvefit(pe,[1 1 .1],tdata,ydata)
pe = @(x)avecExp(1) + avecExp(2)*exp(avecExp(3)*(x - 1850));

fe = figure;
fplot(pe,[1850 2100])
hold on;
scatter(tdata, ydata, 10,[249/256 110/256 59/256], 'filled')
title('Exponential Fit: a + be^{ct}')
xlabel('Year')
ylabel('Temperature Difference (C)')
hold off;
print -depsc ExpFit.eps
% large error because exp can't be negative

est2100Exp = pe(2100)

%%
[avecMixLin resErrorMixLin] = lsqcurvefit(pmixLin,[1 1 .01],tdata,ydata)
pmixLin = @(x)(avecMixLin(1) + ...
    avecMixLin(2)*(x - 1850))*exp(avecMixLin(3)*(x - 1850));

fmix = figure;
fplot(pmixLin,[1850 2100])
hold on;
scatter(tdata, ydata, 10,[249/256 110/256 59/256], 'filled')
title('Mixed Linear Fit: p_1(t)e^{at}')
xlabel('Year')
ylabel('Temperature Difference (C)')
hold off;
print -depsc MixFit1.eps

est2100MixLin = pmixLin(2100)

%%
[avecMixQuad resErrorMixQuad] = ...
    lsqcurvefit(pmixQuad,[1 1 1 .01],tdata,ydata)
pmixQuad = @(x)(avecMixQuad(1) + avecMixQuad(2)*(x - 1850)...
    + avecMixQuad(3)*(x - 1850)^2)*exp(avecMixQuad(4)*(x - 1850));

fmix2 = figure;
fplot(pmixQuad,[1850 2100])
hold on;
scatter(tdata, ydata, 10,[249/256 110/256 59/256], 'filled')
title('Mixed Quadratic Fit: p_2(t)e^{at}')
xlabel('Year')
ylabel('Temperature Difference (C)')
hold off;
print -depsc MixFit2.eps

est2100MixQuad = pmixQuad(2100)

%% 
fOrigin = figure;
scatter(tdata, ydata, 10,[249/256 110/256 59/256], 'filled')
title('Temperature Anomalies Above Sea Level')
xlabel('Year')
ylabel('Temperature Difference (C)')
hold off;
print -depsc RawData.eps