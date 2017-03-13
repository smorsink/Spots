clear%, close all

cd /Users/kitung/Desktop/thesis_material/Spot/team_hw

temp = 5.69; %5.69 6.05 6.3
angle = 9090; % 9090 3080
ene = 2; % 1 to 15
typ = 1; % 1 is monochromatic, 2 is integrated

ener = ene+1;
if (typ == 1)
    type = 'mono';
    titl_name = ['monochromatic at ', num2str(ene*0.2), ' keV'];
    y_name = 'counts/kev';
elseif (typ == 2)
    type = 'inte';
    titl_name = ['integrated from ', num2str(ene*0.2-0.1), ' to ', num2str(ene*0.2+0.1), ' keV'];
    y_name = 'counts';
end

mcphac = load(['mcphac_',num2str(temp),'_',num2str(angle),'_',type,'.txt']);
nsxh = load(['nsxh_',num2str(temp),'_',num2str(angle),'_',type,'.txt']);
nsatmos = load(['nsatmos_',num2str(temp),'_',num2str(angle),'_',type,'.txt']);
nsxhe = load(['nsxhe_',num2str(temp),'_',num2str(angle),'_',type,'.txt']);

plot(nsxh(:,1),nsxh(:,ener),'linewidth',3);
hold on
plot(mcphac(:,1),mcphac(:,ener),'linewidth',3);
plot(nsatmos(:,1),nsatmos(:,ener),'linewidth',3);
plot(nsxhe(:,1),nsxhe(:,ener),'linewidth',3);
xlabel('phase')
ylabel(y_name)
title(titl_name)
ax = gca;
ax.FontSize = 24;
%axis([0 0.1 120 140])
%axis([0 0.1 450 600])

%nsxh as f2, the 'model'
%three columns are mcphac, nsatmos, and nsxhe
%compare only first phase bin
abs_diff = [mcphac(1,ener)-nsxh(1,ener), nsatmos(1,ener)-nsxh(1,ener), nsxhe(1,ener)-nsxh(1,ener)];
fra_diff = [(mcphac(1,ener)-nsxh(1,ener))/mcphac(1,ener), (nsatmos(1,ener)-nsxh(1,ener))/nsatmos(1,ener), (nsxhe(1,ener)-nsxh(1,ener))/nsxhe(1,ener)];

dex = [];
for i = 1:32
    if(nsxh(i,ener)>0)
        dex = [dex i];
    end
end

dex_size = size(dex);

chi_diff = [sum((mcphac(dex,ener)-nsxh(dex,ener)).^2./mcphac(dex,ener)), sum((nsatmos(dex,ener)-nsxh(dex,ener)).^2./nsatmos(dex,ener)), sum((nsxhe(dex,ener)-nsxh(dex,ener)).^2./nsxhe(dex,ener))];

for i = dex
    if((mcphac(i,ener)-nsxh(i,ener))^2/mcphac(i,ener) > 5*mean((mcphac(dex,ener)-nsxh(dex,ener)).^2./mcphac(dex,ener)))
        disp(['mcphac anomaly at ',num2str(i),' bin'])
    end
    
    if((nsatmos(i,ener)-nsxh(i,ener))^2/nsatmos(i,ener) > 5*mean((nsatmos(dex,ener)-nsxh(dex,ener)).^2./nsatmos(dex,ener)))
        disp(['nsatmos anomaly at ',num2str(i),' bin'])
    end
    
    if((nsxhe(i,ener)-nsxh(i,ener))^2/nsxhe(i,ener) > 5*mean((nsxhe(dex,ener)-nsxh(dex,ener)).^2./nsxhe(dex,ener)))
        disp(['nsxhe anomaly at ',num2str(i),' bin'])
    end
end

diff_matrix = [abs_diff; fra_diff; chi_diff];
disp('    mcphac   nsatmos    nsxhe')
disp(diff_matrix)
disp('row 1: absolute difference')
disp('row 2: fractional difference')
disp('row 3: chi-squared')
disp([num2str(dex_size(2)),' phase bins are nonzero.'])

legend('nsxh','mcphac','nsatmos','nsxhe','location','best')
hold off


