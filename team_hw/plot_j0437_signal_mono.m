clear, close all

data = load('jul4_nsxhnew_obl_partial_mono.txt');
data1 = load('jul4_nsxhpi_obl_partial_mono.txt');
data2 = load('jul4_nsxhe_obl_partial_mono.txt');
ener = linspace(0.095,3.105,301);
t = linspace(0,1-1/16,16);

figure,
for k = 1:100
i = rem(k,16)+1;
loglog(ener,data(i,2:end))
hold on,
loglog(ener,data1(i,2:end))
loglog(ener,data2(i,2:end))
axis([1e-1 5 1e-2 2e1])
legend('nsxh','nsxhpi','nsxhe','location','southwest')
hold off
pause(0.1)
end
%
figure,
for i = [2:50 10]
plot(t,data(:,i))
hold on
plot(t,data1(:,i))
plot(t,data2(:,i))
plot(t,data2(:,i)-data1(:,i))
axis([0 1 0 10])
legend('nsxh','nsxhpi','nsxhe','location','best')
hold off
pause(0.03)
end
figure,
for i = [51:200 51]
plot(t,data(:,i))
hold on
plot(t,data1(:,i))
plot(t,data2(:,i))
plot(t,data2(:,i)-data1(:,i))
hold off
axis([0 1 0 5])
legend('nsxh','nsxhpi','nsxhe','location','best')
pause(0.03)
end
figure,
for i = [201:301 201]
plot(t,data(:,i))
hold on
plot(t,data1(:,i))
plot(t,data2(:,i))
plot(t,data2(:,i)-data1(:,i))
hold off
axis([0 1 0 0.1])
legend('nsxh','nsxhpi','nsxhe','location','best')
pause(0.03)
end




