clear, close all

defltable=fopen('Hatm8000dT0.05.bin');
a = fread(defltable,'double');   
fclose(defltable);
data = load('nsx_H_v170524_no_header.txt');
logg = 14.2;

for i = 1:15
logt = 5+i*0.1;
mgc = floor((logg-13.7)/0.1);
mtc = floor((logt-5.1)/0.05);
mcphac_first = (mtc*11+mgc)*5000;
ngc = floor((logg-13.7)/0.1);
ntc = floor((logt-5.1)/0.1);
nsxh_first = (ntc*11+ngc)*137;

cd /Users/kitung/Desktop/thesis_material/Spot/atmosphere

ev = 1.60217653e-12;
h = 6.6260693e-27;
k = 1.38e-16;

norm_ener = a([1:1595000]*5-2);
inte = a([1:1595000]*5);

loglog(10.^norm_ener(mcphac_first+linspace(50,5000,100))/ev/1000*k*10^logt,inte(mcphac_first+linspace(50,5000,100))*10^(3*logt))
xlabel('energy (keV)')
ylabel('intensity (cgs units)')
hold on
loglog(10.^(linspace(-1.32,1.4,137)),data(nsxh_first+[1:137],1))
hold off
axis([1e-2 1e2 1e-4 1e4])
pause(0.5)
end