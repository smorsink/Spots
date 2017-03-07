clear, close all

%cd /Users/kitung/Desktop/thesis_material/Spot/new_atmo_test
cd /home/kitung/Spot/new_atmo_test

%{
data = load('first_test_bb.txt');
figure,
plot(data(:,1),data(:,5))
%}

%{
data = load('first_test_nsatmos.txt');
figure,
plot(data(:,1),data(:,5))
hold on,

data = load('first_test_nsxhe.txt');
plot(data(:,1),data(:,5))
hold off
%}

data = load('nsatmos_mono.txt');
figure,
plot(data(:,1),data(:,2))
hold on
data = load('nsxhe_mono.txt');
plot(data(:,1),data(:,2))
data = load('nsxh_mono.txt');
plot(data(:,1),data(:,2))
hold off
legend('nsatmos','nsxhe','nsxh')



for i = 1:30
%for i = linspace(1,29,15)
    data = load('nsatmos_inte.txt');
    figure,
    plot(data(:,1),data(:,i+1))
    hold on

    data = load('nsxh_inte.txt');
    plot(data(:,1),data(:,i+1))
    
    data = load('nsxhe_inte.txt');
    plot(data(:,1),data(:,i+1))
    
    legend('nsatmos','nsxh','nsx helium')
    name = [num2str(0.1*i-0.05),' keV to ', num2str(0.1*i+0.05), ' keV'];
    title(name)
    hold off
end


%{
data = load('nsatmos_inte.txt');
figure,
plot(data(:,1),data(:,2))
hold on

data = load('nsxh_inte.txt');
plot(data(:,1),data(:,2))
    
data = load('nsxhe_inte.txt');
plot(data(:,1),data(:,2))
    
legend('nsatmos','nsxh','nsx helium')
hold off
%}
