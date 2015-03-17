clear all 
close all

mydirnames = dir('aE0.9*')


  for dircntr=1:length(mydirnames)
		cd(mydirnames(dircntr).name)

thisqval = 0.1;

data = load('fractioninitvsfinal.pco');
dataf = load('fractions.pco');
cellcontents = load('cellcontents.pco');

ks = 1;
phi = 0.1;
vol = 1e-19;
%NZ = 664.929;


nrX = 101;
nrY = 101;


NAVOG =  6.02e23;
%concZ = NZ/NAVOG/vol;

%X = data(1:end,1);
%Y = data(1:end,2);
X = data(:,1);
Y = data(:,2);
Z = data(:,3);

X = reshape(X,[nrY,nrX]);
Y = reshape(Y,[nrY,nrX]);
Z = reshape(Z,[nrY,nrX]);

Zm = min(Z,0.01);

% h = figure(1);
% A = axes;
% plot3(as,qs,Z,'b-')
h = figure(1);
A = axes;
surfc(X,Y,Zm)
xlabel('fraction C initial')
ylabel('fraction C daughter cells')
zlabel('n')
view(0,90)
%set(A,'XScale','log');
%ylim([min(qs),max(qs)])
colorbar 

figure(2);
X2 = zeros([51,51]);
Y2 = zeros([51,51]);
Z2 = zeros([51,51]);
for i=1:50
    for j=1:50
        X2(i,j) = sum(sum(X((2*i-1):(2*i),(2*j-1):(2*j))))/4;
        Y2(i,j) = sum(sum(Y((2*i-1):(2*i),(2*j-1):(2*j))))/4;
        Z2(i,j) = sum(sum(Z((2*i-1):(2*i),(2*j-1):(2*j))));
    end
    X2(i,51) = sum(X((2*i-1):(2*i),101))/2;
    Y2(i,51) = sum(Y((2*i-1):(2*i),101))/2;
    Z2(i,51) = sum(Z((2*i-1):(2*i),101));
end
for j=1:50
    X2(51,j) = sum(X(101,(2*j-1):(2*j)))/2;
    Y2(51,j) = sum(Y(101,(2*j-1):(2*j)))/2;
    Z2(51,j) = sum(Z(101,(2*j-1):(2*j)));
end
X2(51,51)=X(101,101);
Y2(51,51)=Y(101,101);
Z2(51,51)=Z(101,101);

filename = 'fig1.pdf';
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-dpdf', filename)



Z2 = min(Z2,0.01);

A = axes;
%subplot(1,2,1)
%mesh(X2,Y2,Z2)
%xlabel('fraction C initial')
%ylabel('fraction C daughter cells')
%zlabel('n')
%view(0,90)
%colorbar 

%subplot(1,2,2)
imagesc([0:2:100],[0:2:100],Z2);
set(gca,'YDir','normal');
xlabel('fraction C initial')
ylabel('fraction C daughter cells')
zlabel('n')
colormap(flipud(hot))
colorbar 

hold on
plot3([0 100],[0 100], [1 1], 'k:');

filename = 'fig2.pdf';
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-dpdf', filename)


figure(3)
fx = dataf(:,1);
fy = dataf(:,4);
plot(fx,fy,'b-');
hold on
plot([0:1:100],sum(Z,2),'r:')
xlabel('fraction C in cell after division')
ylabel('fraction of cells')

filename = 'fig3.pdf';
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-dpdf', filename)

figure(4)
bar([0:1:100],sum(Z,1))
xlim([0 100])
xlabel('initial fraction C in cell')
ylabel('fraction of cells reaching division')

filename = 'fig4.pdf';
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-dpdf', filename)


figure(5)
fx = dataf(:,1);
fy = dataf(:,4);
bar(fx,fy,'b');
xlabel('fraction C in cell')
ylabel('fraction of cells in system')

filename = 'fig5.pdf';
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-dpdf', filename)

figure(6)
m=25;
%myhot = hot(66);
%myhot(64:65,:)=[];
myhot = hot(130);
myhot(128:129,:)=[];
myhot = hot(256);
myhot = flipud(myhot);
cellcontents = reshape(cellcontents,[101,101]);
%cellcontents(m,m) = -0.007*max(max(cellcontents));
imagesc([0:1:m-1],[0:1:m-1],log10(cellcontents(1:m,1:m)));
%plot3(X2,Y2,cellcontents(1:m+1,1:m+1));
set(gca,'YDir','normal');
xlabel('nr B')
ylabel('nr E')
zlabel('n')
%colormap(flipud(hot))
colormap(myhot)
colorbar 

filename = 'fig6.pdf';
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-dpdf', filename)

nrs = zeros([m,m]);
for i = 1:m
    for j=1:m
        nrs(i,j) = i+j-2;
    end
end


figure(61)
mymaxcval = 0.02
cellcontentsmass = cellcontents(1:m,1:m).*nrs(1:m,1:m);
masstot = sum(sum(cellcontentsmass));
cellcontentsmassM = min(cellcontentsmass/masstot,mymaxcval);
imagesc([0:1:m-1],[0:1:m-1],cellcontentsmassM);
set(gca,'YDir','normal');
set(gca,'FontSize',24)
set(gca, 'box', 'off')
set(gca,'TickDir','out')
set(gca, 'linewidth', 1.5)
caxis([0.0, mymaxcval])
axis image
%ht=text(18,20,['u=',num2str(1-thisqval)])
  %ht=text(18,20,'u=0.10')
%	   set(ht,'FontSize',28)
	   xlabel('nr B')
	   ylabel('nr E')
	   zlabel('n')
	   colormap(myhot)
	   %cbax = colorbar;
%set(cbax, 'fontsize', 20)
%set(cbax,'location', 'NorthOutside')

filename = 'fig6c.pdf';
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-dpdf', filename)

figure(69)
imagesc([0:1:m-1],[0:1:m-1],cellcontentsmassM);
set(gca,'YDir','normal');
set(gca,'FontSize',24)
set(gca, 'box', 'off')
set(gca,'TickDir','out')
set(gca, 'linewidth', 1.5)
%	   set(ht,'FontSize',28)
	   xlabel('nr B')
	   ylabel('nr E')
	   zlabel('n')
	   colormap(myhot)
	   cbax = colorbar;
set(cbax, 'fontsize', 20)
set(cbax,'location', 'NorthOutside')

filename = 'fig6c_withbar.pdf';
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-dpdf', filename)

figure(62)
imagesc([0:1:m-1],[0:1:m-1],log10(cellcontentsmass/masstot));
set(gca,'YDir','normal');
xlabel('nr B')
ylabel('nr E')
zlabel('n')
colormap(myhot)
colorbar 

filename = 'fig6d.pdf';
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-dpdf', filename)

figure(14)
imagesc([0:1:m-1],[0:1:m-1],cellcontents(1:m,1:m)>0);
set(gca,'YDir','normal');
xlabel('nr B')
ylabel('nr E')
colormap(myhot)
colorbar 

filename = 'fig6b.pdf';
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-dpdf', filename)


% figure(7) % fitness landscape
% % the payoff for C is: F_C(i,j)= [R(i-1)+Sj] / (i+j-1)
% % the payoff for D is: F_D(i,j)= [Ti+P(j-1)] / (i+j-1)
% % we can use the simplified payoff matrix:
% % R=b-c  S=-c  T=b  P=0  with b>c>0
% % the fitness of each type is 
% % f_C= exp (beta F_C)
% % f_D= exp (beta F_D)
% 
% b = 5;
% c = 1;
% beta = 1;
% R = b-c;  S = -c;  T = b;  P = 0;
% 
% fitnessC = zeros([m,m]);
% fitnessD = zeros([m,m]);
% fitness = zeros([m,m]);
% for n=1:m
%     for i = 0:n
%         nC = i;
%         nD = n-i;
%         fitnessC(nC+1,nD+1) = exp(beta*(R*(nC-1)+S*nD)/max(nC+nD-1,1));
%         fitnessD(nC+1,nD+1) = exp(beta*(T*nC+P*(nD-1))/max(nC+nD-1,1));
%         fitness(nC+1,nD+1) = nC*fitnessC(nC+1,nD+1) + nD*fitnessD(nC+1,nD+1);
%     end
% end
% %fitness(m,m) = -0.007*max(max(fitness));
% imagesc([0:1:m-1],[0:1:m-1],fitness);
% %plot3(X2,Y2,fitness);
% set(gca,'YDir','normal');
% xlabel('nr B')
% ylabel('nr E')
% colormap(flipud(hot))
% colorbar 
% title('fitness cell')

% figure(8)
% %fitnessC(m,m) = -0.007*max(max(fitnessC));
% imagesc([0:1:m-1],[0:1:m-1],fitnessC);
% set(gca,'YDir','normal');
% xlabel('nr B')
% ylabel('nr E')
% colormap(flipud(hot))
% colorbar 
% title('fitness C')
% 
% %filename = 'fig8.pdf';
% %set(gcf, 'PaperPositionMode', 'auto');
% %print(gcf, '-dpdf', filename)
% 
% figure(9)
% %fitnessD(m,m) = -0.007*max(max(fitnessD));
% imagesc([1:1:m-1],[1:1:m-1],fitnessD(2:m,2:m));
% %plot3(X2,Y2,fitnessD);
% set(gca,'YDir','normal');
% xlabel('nr B')
% ylabel('nr E')
% colormap(flipud(hot))
% colorbar 
% title('fitness D')
% 
% %filename = 'fig9.pdf';
% %set(gcf, 'PaperPositionMode', 'auto');
% %print(gcf, '-dpdf', filename)
% 
% 
% figure(10)
% growth = ks*fitness(1:m,1:m).*cellcontents(1:m,1:m).*concZ;
% imagesc([0:1:m-1],[0:1:m-1],growth);
% %plot3(X2,Y2,fitnessD);
% set(gca,'YDir','normal');
% xlabel('nr B')
% ylabel('nr E')
% colormap(flipud(hot))
% colorbar 
% title('growth')
% 
% sum(sum(growth))
% death = nrs.*cellcontents(1:m,1:m).*phi;
% sum(sum(death))
% 
% %filename = 'fig10.pdf';
% %set(gcf, 'PaperPositionMode', 'auto');
% %print(gcf, '-dpdf', filename)
% 
% 
% figure(11)
% nrsnonzero = nrs;
% nrsnonzero(find(nrsnonzero==0))=1;
% fitnesspermolincell = fitness(1:m,1:m)./nrsnonzero;
% imagesc([0:1:m-1],[0:1:m-1],fitnesspermolincell(1:m,1:m));
% set(gca,'YDir','normal');
% set(gca,'FontSize',32)
% set(gca, 'box', 'off')
% set(gca,'TickDir','out')
% set(gca, 'linewidth', 1.5)
% xlabel('nr B')
% ylabel('nr E')
% colormap(flipud(hot))
% cbax = colorbar;
% set(cbax, 'fontsize', 26)
% set(cbax,'location', 'NorthOutside')
% %title('average fitness per molecule')
% 
% filename = 'fig11.pdf';
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, '-dpdf', filename)
% 
% figure(12)
% q = ks*concZ*fitnesspermolincell(1:m,1:m)-phi;
% sum(sum(q.*nrs.*cellcontents(1:m,1:m)))
% imagesc([0:1:m-1],[0:1:m-1],q);
% hold on
% contour(X2(1:m,1:m),Y2(1:m,1:m),q,[0,0],'k');
% set(gca,'YDir','normal');
% xlabel('nr B')
% ylabel('nr E')
% colormap(flipud(hot))
% colorbar 
% title('growth positive or negative?')
% 
% filename = 'fig12.pdf';
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, '-dpdf', filename)


% figure(14)
% imagesc([0:1:m-1],[0:1:m-1],q>0);
% set(gca,'YDir','normal');
% set(gca,'FontSize',32)
% set(gca, 'box', 'off')
% set(gca,'TickDir','out')
% set(gca, 'linewidth', 1.5)
% xlabel('nr B')
% ylabel('nr E')
% colormap(flipud(hot))
% cbax = colorbar;
% set(cbax, 'fontsize', 26)
% set(cbax,'location', 'NorthOutside')
% %title('growth positive or negative?')
% 
% filename = 'fig14.pdf';
% set(gcf, 'PaperPositionMode', 'auto');
% print(gcf, '-dpdf', filename)


figure(13);
Z3 = Z2 + 0.0;
t = zeros([51,1]);
for i=1:51
    t(i) = sum(Z3(:,i));
    if t(i) > 0
       Z3(:,i) = Z3(:,i)/t(i);
    end
end

Z3 = min(Z3,0.25);

A = axes;
%surfc(X,Y,Z3)
imagesc([0:2:100],[0:2:100],Z3)
set(gca,'YDir','normal');
set(gca,'FontSize',32)
set(gca, 'box', 'off')
set(gca,'TickDir','out')
set(gca, 'linewidth', 1.5)
xlabel('fraction C initial')
ylabel('fraction C daughter cells')
colormap(flipud(hot))
cbax = colorbar;
set(cbax, 'fontsize', 26)
set(cbax,'location', 'NorthOutside')
hold on
plot3([0 100],[0 100], [1 1], 'k:');

filename = 'fig13.pdf';
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-dpdf', filename)

cd('..')
end
