close all
clc
clear all
x = linspace(0,1/30);
n = 1;
hold on

while n < 13
    plot(x,165/pi^2*1/n^2*sin(n*pi/2)*sin(2*pi*60*n*x))
    n = n+2;
end
xlabel('Time (s)');
ylabel('Current (A)');
figure()
sum = zeros(1,100);
i = zeros(1,100);

for t = 1:100
    for n = 1:2:13
        sum(t) = sum(t) + 1/n^2 * sin(n*pi/2) * sin(2*pi*n*60*x(t));
    end
    i(t) = 165/pi^2 * sum(t);
end
plot(x,i)
xlabel('Time (s)');
ylabel('Current (A)');

amplitude = zeros(1,7);
harmonic = zeros(7,100);
k = 1;
for n = 1:2:13
    for t = 1:100
        harmonic(k,t) = 165/pi^2*1/n^2*sin(n*pi/2)*sin(2*pi*n*60*x(t));
    end
    amplitude(k)=max(harmonic(k,:));
    k = k + 1;
end

Irms = zeros(1,7);
THD = zeros(1,7);
THD(1) = 1;
for k = 1:7
    Irms(k) = amplitude(k) / sqrt(3);
end

for k = 2:7
    j = 2;
    while j <= k
        THD(k) = THD(k) + Irms(j)^2;
        j = j + 1;
    end
    THD(k) = sqrt(THD(k))/Irms(1);
end

convergence = zeros(1,6);
for k = 7:-1:2
    convergence(k-1) = abs((THD(k)-THD(k-1))/THD(k)) * 100;
end
figure()
x = 5:2:13;
plot(x,convergence(2:6))
xlabel('Harmonic');
ylabel('Convergence (%)')