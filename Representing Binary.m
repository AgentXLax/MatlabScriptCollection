%{
Nathan Wisla
Numerical Methods (AUMAT 340)
January 20, 2016
%}

binRep123 = binaryArray(123,8);    %Q1a
ternRep123 = ternArray(123,8);     %Q1b
short7Rep = shortBin(7);           %Q2a
short17Rep = shortBin(17);
sumReps = sumBinaryArray(7,17);    %Q2b
compareReps = shortBin(24);

recursive = secondOrder(3);        %Q3

%{
Q4 Using the command window, the highest value for x in a = 10^x
   is 308 where 'a' gets mapped to infinity. The smallest value is
   -323 where going beyond, 'a' gets mapped to 0. 
   This is consistent with the in-class notes where those numbers 
   would exceed 64 bits of information.
   NOTE: after 'a' goes beyond 64 bits, it is too big to be 
   stored into memory and therefore cannot be reduced or increased
   to a new number unless reassigned entirely.
%}
a1 = 10^308;% == 1.000e308
a2 = 10^309;% == Inf
a2_1 = a2 * 10^-30;% == Inf
a3 = 10^-323;% == 1.000e-323
a4 = 10^-324;% == 0
a4_1 = a4 * 10^45; % == 0

bin2Float = IEEE757('01001110110001111001011100001010'); %Q5
    
function b = binaryArray(dec,bit)      %Question 1a)
    b = mod(floor(dec * pow2(1 - bit:0)),2);%based on the dec2bin function
end
function t = ternArray(dec,digit)      %Question 1b)
    ternNum = 1 - digit:0;%creates an array from -i to 0
    i = digit;
    while i > 0
        ternNum(i) = mod(floor(dec * 3 ^ ternNum(i)),3);
        i = i - 1;
    end
    t = ternNum;
end

function s = shortBin(dec)             %Question 2a)
    [f,e] = log2(max(dec));%largest amount of digits needed (from dec2bin)
    s = mod(floor(dec * pow2(1 - max(e):0)),2);
end
function a = sumBinaryArray(dec1,dec2) %Question 2b)
    [f,i] = log2(max(dec1 + dec2));%i displays the index of the binary representation of the number
    binRep1 = mod(floor(dec1 * pow2(1 - max(i):0)),2);
    binRep2 = mod(floor(dec2 * pow2(1 - max(i):0)),2);
    sumRep = binRep1 + binRep2;
    
    while i > 0
        if sumRep(i) > 1 %if a bit is > 1
            sumRep(i-1) = sumRep(i-1) + 1; %shift the higher bit up 1 
            sumRep(i) = mod(sumRep(i),2); %and clock the lower bit back
        end
        i = i - 1;
    end
    a = sumRep;    
end

function y = secondOrder(n)            %Question 3
    if n == 0 %recursion end step
        y = 1 / 4;
    else %recursive function
        y = 1 / (4 + secondOrder(n - 1));
    end
end

function x = IEEE757(bin)              %Question 5 input binary number as char array
    array = bin == '1'; %the character array into a boolean array
    places = ones(1,32); %denotes the position of the place
    %construct the places array
    for n = 1:32
        if n < 32
            places(1,33 - n) = 2^(n - 24); %place 2^n indexed appropriately to IEEE format
        else
            places(1,33 - n) = (-1)^array(1,1); %place (-1)^n in the first entry
        end
    end
    sign = places(1,1); %create the sign
    exp = 2^((array(2:9)*places(2:9)') - 127);%matrix multiply the exponent bits
    mant = 1 + array(10:32)* places(10:32)';%matrix multiply the mantissa

    x = sign * exp * mant;

end
    