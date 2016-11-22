clear, close all

cd /Users/kitung/Desktop/thesis_material/Spot/SpotShape

osn = load('1_sml_no.txt');
tsn = load('2_sml_no.txt');
oln = load('1_lrg_no.txt');
tln = load('2_lrg_no.txt');
figure,
plot(osn(:,1),osn(:,2))
figure,
plot(tsn(:,1),tsn(:,2))
figure,
plot(oln(:,1),oln(:,2))
figure,
plot(tln(:,1),tln(:,2))



sp = load('sml_pole.txt');
lp = load('lrg_pole.txt');
figure,
plot(sp(:,1),sp(:,2))
figure,
plot(lp(:,1),lp(:,2))

