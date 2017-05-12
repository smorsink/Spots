clear, close all

cd /home/kitung/Spot/input
data = load('mcphacc_may12.txt');
data0 = data;
bg = load('mcphacc_may12_hbg.txt');
data(:,2:16) = data(:,2:16)+bg(:,2:16);

%adding noise goes here

noise = imnoise(data(:,2:end)*1e-12,'poisson');
ndat = data;

data(:,2:end) = noise*1e12;
data(:,2:end) = floor(data(:,2:end)+1);


newdata = zeros(32,31);
newdata(:,1) = data(:,1);
for i = 1:15
    newdata(:,2*i) = round(data(:,i+1));
    newdata(:,2*i+1) = sqrt(round(data(:,i+1)));
end

fid = fopen('mcphacc_may12_withnoise_bgh.txt','w');
fprintf(fid, '%.6f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f %d %.4f\n',newdata');
fclose(fid)

fid = fopen('mcphacc_may12_bgnormalized.txt','w');
fprintf(fid, '%1.4g %1.4g %1.4g %1.4g %1.4g %1.4g %1.4g %1.4g %1.4g %1.4g %1.4g %1.4g %1.4g %1.4g %1.4g',bg(1,2:16)/1e7*32');
fclose(fid)

figure,
plot(data0(:,2))
hold on
plot(ndat(:,2))
plot(data(:,2))
