data = [34 35 45 40 46 38 47 36 38 34 33 36 43 43 37 38 32 38 40 33 38 40 48 39 32 36 40 40 36 34];
% all data points are in arcminutes
sample1 = [46 36 38 43 37];
sample2 = [35 38 33 39 47];
histogram(data);
mean = avg(data);
median = mid(data);
[Mode,F] = mode(data);

%Deviation = sqrt(sum((mean - x_i)^2)/N)
SD = standDev(data, mean);
SDM = standDevMean(data, mean);
SD2 = standDev2(data, mean);
sum = standSum(data,mean);
avgSample1 = avg(sample1);
s1SD = standDev(sample1,avgSample1);
avgSample2 = avg(sample2);
s2SD = standDev(sample1,avgSample2);

function m = avg(data)
    sum = 0;
    for i = 1:length(data)
        sum = sum + data(i);
    end
    m = sum / length(data);      
end
function m = mid(data)
    data = vec2mat(sort(data), ceil(length(data) / 2));
    m = (data(1,length(data)) + data(2,1)) / 2;
end
function S = standDev(data, mean)
    sum = 0;
    for i = 1:length(data)
        sum = sum + (mean - data(i))^2;
    end
    S = sqrt(sum / (length(data) - 1));
end%4
function S = standDevMean(data, mean)
    SD = standDev(data, mean);
    S = SD / sqrt(length(data));
end
function S = standDev2(data, mean)
    SD = standDev(data, mean);
    S = SD / sqrt(2 * (length(data) - 1));
end
function S = standSum(data, mean)
    sum = 0;
    for i = 1:length(data)
        sum = sum + (mean - data(i))^2;
    end
    S = sum;
end

