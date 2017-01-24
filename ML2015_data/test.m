clear, close all

cd /Users/kitung/Desktop/thesis_material/Spot/ML2015_data

fake = load('newa-10-err.txt');
data = load('newa-10.txt');

plot(fake(:,1),fake(:,2));
hold on
plot(data(:,1),data(:,2));