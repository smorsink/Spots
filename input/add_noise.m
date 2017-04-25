clear, close all

cd /Users/kitung/Desktop/thesis_material/Spot/input
data = load('nsxh_apr_19_inte.txt');


%adding noise goes here

noise = imnoise(data(:,2:end)*1e-12,'poisson');
ndat = data;

data(:,2:end) = noise*1e12;


newdata = zeros(32,31);
newdata(:,1) = data(:,1);
for i = 1:15
    newdata(:,2*i) = round(data(:,i+1));
    newdata(:,2*i+1) = sqrt(round(data(:,i+1)));
end

fid = fopen('nsxh_apr_19_withnoise.txt','w');
%fid = fopen('nsxh_apr_19_uncert.txt','w');
fprintf(fid, '%.6f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f\n',newdata');
fclose(fid)

figure,
plot(ndat(:,2))
hold on
plot(data(:,2))

data = load('SMMtest.txt');
data1 = load('SMMpoisson.txt');
noise = imnoise(data(1:32,3)*1e-12,'poisson')*1e12;

figure,
plot(data(1:32,1),data(1:32,3))
hold on
plot(data(1:32,1),noise)
figure,
plot(data(1:32,1),data(1:32,3))
hold on
plot(data1(1:32,1),data1(1:32,3))

