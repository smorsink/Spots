clear, close all

data = load('jun6_nsxhnew_obl_j0437_mono.txt');
data1 = load('jun6_mcphacc_obl_j0437_mono.txt');

for i = 1:16
loglog(data(i,2:end))
hold on
loglog(data1(i,2:end))
hold off
axis([1 400 1 2e3])
legend('nsxhnew','mcphac','location','best')
pause(0.2)
end
pause(0.5)
for i = 2:50
plot(data(:,i))
hold on
plot(data1(:,i))
hold off
axis([0 16 0 2e3])
legend('nsxhnew','mcphac','location','best')
pause(0.05)
end
pause(0.5)
for i = 51:200
plot(data(:,i))
hold on
plot(data1(:,i))
hold off
axis([0 16 0 2e2])
legend('nsxhnew','mcphac','location','best')
pause(0.05)
end
pause(0.5)
for i = 201:301
plot(data(:,i))
hold on
plot(data1(:,i))
hold off
axis([0 16 0 2e1])
legend('nsxhnew','mcphac','location','best')
pause(0.05)
end
