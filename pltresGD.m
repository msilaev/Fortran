%%% This Maltab code plots Gilbert damping as function of temperature 
%%%  in FI/S'S structure

clear all

h1=figure
h2=figure

for indF=[8 14  16]

for indDiff=10:10


   if (indF==1)
        dRsuper=1.5;
%fnameComm='h3Rmax021D3d0DN1d5DiffNonConstCplVar';
fnameComm=['h3Rmax021D3d0DN1d5Diff' num2str(indDiff) 'CplVar'];
     end
      if (indF==2)
          dRsuper=1.4;
fnameComm='h3Rmax021D3d0DN1d4DiffNonConstCplVar';
fnameComm=['h3Rmax021D3d0DN1d4Diff' num2str(indDiff) 'CplVar'];
      end
       if (indF==3)
         dRsuper=1.3;
fnameComm='h3Rmax021D3d0DN1d3DiffNonConstCplVar';
fnameComm=['h3Rmax021D3d0DN1d3Diff' num2str(indDiff) 'CplVar'];
       end
        if (indF==4)
          dRsuper=1.2;
 fnameComm='h3Rmax021D3d0DN1d2DiffNonConstCplVar';
 fnameComm=['h3Rmax021D3d0DN1d2Diff' num2str(indDiff) 'CplVar'];
        end
         if (indF==5)
          dRsuper=1.1;
 fnameComm='h3Rmax021D3d0DN1d1DiffNonConstCplVar';
 fnameComm=['h3Rmax021D3d0DN1d1Diff' num2str(indDiff) 'CplVar'];
         end
          if (indF==6)
           dRsuper=1.0;
 fnameComm='h3Rmax021D3d0DN1d0DiffNonConstCplVar';
 fnameComm=['h3Rmax021D3d0DN1d0Diff' num2str(indDiff) 'CplVar'];
            end          
            if (indF==7)
           dRsuper=0.9;
 fnameComm='h3Rmax021D3d0DN0d9DiffNonConstCplVar';
 fnameComm=['h3Rmax021D3d0DN0d9Diff' num2str(indDiff) 'CplVar'];
            end
            if (indF==8)
           dRsuper=0.8;
 fnameComm='h3Rmax021D3d0DN0d8DiffNonConstCplVar';
 fnameComm=['h3Rmax021D3d0DN0d8Diff' num2str(indDiff) 'CplVar'];
            end
               if (indF==9)
              dRsuper=0.7;
 fnameComm='h3Rmax021D3d0DN0d7DiffNonConstCplVar';
 fnameComm=['h3Rmax021D3d0DN0d7Diff' num2str(indDiff) 'CplVar'];
               end
               if (indF==10)
            dRsuper=0.6;
 fnameComm='h3Rmax021D3d0DN0d6DiffNonConstCplVar';
 fnameComm=['h3Rmax021D3d0DN0d6Diff' num2str(indDiff) 'CplVar'];
               end
               if (indF==11)
           dRsuper=0.5;
 fnameComm='h3Rmax021D3d0DN0d5DiffNonConstCplVar';
 fnameComm=['h3Rmax021D3d0DN0d5Diff' num2str(indDiff) 'CplVar'];
               end
               if (indF==12)
           dRsuper=0.4;
 fnameComm='h3Rmax021D3d0DN0d4DiffNonConstCplVar';
 fnameComm=['h3Rmax021D3d0DN0d4Diff' num2str(indDiff) 'CplVar'];
               end
               if (indF==13)
           dRsuper=0.3;
 fnameComm='h3Rmax021D3d0DN0d3DiffNonConstCplVar';
%fnameComm='h3Rmax041D3d0DN0d3DiffNonConstCplVar';
fnameComm=['h3Rmax021D3d0DN0d3Diff' num2str(indDiff) 'CplVar'];
               end
               if (indF==14)
           dRsuper=0.2;
 fnameComm='h3Rmax021D3d0DN0d2DiffNonConstCplVar';
 fnameComm=['h3Rmax021D3d0DN0d2Diff' num2str(indDiff) 'CplVar'];
               end
               if (indF==15)
           dRsuper=0.1;
 fnameComm='h3Rmax021D3d0DN0d1DiffNonConstCplVar';
 fnameComm=['h3Rmax021D3d0DN0d1Diff' num2str(indDiff) 'CplVar'];
               end
               if (indF==16)
           dRsuper=0.0;
 fnameComm='h3Rmax021D3d0DN0d0DiffNonConstCplVar';
 fnameComm=['h3Rmax021D3d0DN0d0Diff' num2str(indDiff) 'CplVar'];
               end
            
fname0=['c:\Users\Mike\Work\Jyvaskyla2016\Artjom\SpinRelaxation\CheckSpinRelaxation\ResVarDGamma0001\'];
fname0Res=['c:\Users\Mike\Work\Jyvaskyla2016\SpinSuscpetibilityLinear\Inhomogeneous\numericsVarDGamma0001\'];


fileSave=[fname0Res fnameComm 'tSO1ResKin'];

 load(fileSave,'TT', 'DDiss');
DissKin=real(DDiss);

fileSave=[fname0Res fnameComm 'tSO1ResSpec'];

 load(fileSave,'TT', 'DDiss');
DissSpec=real(DDiss);

F=DissKin+DissSpec;

if (min(F)<0)
F= F - min(F);
end

figure(h2)
if(indF==8)    
C=[0.9290    0.6940    0.1250];
%plot(TT/TT(length(TT)),F /F(length(TT)),'LineWidth',2,'Color',C)
plot(TT/TT(length(TT)),F /F(length(TT)),'LineWidth',2)
end
if(indF==14)
%plot(TT/TT(length(TT)),F /F(length(TT)),'LineWidth',2,'Color','m')
plot(TT/TT(length(TT)),F /F(length(TT)),'LineWidth',2)

end
if(indF==16)
%plot(TT/TT(length(TT)),F /F(length(TT)),'LineWidth',2,'Color','r')
plot(TT/TT(length(TT)),F /F(length(TT)),'LineWidth',2)

end
hold on

% figure(h1)
% beta=0*TT + (16-indF)*0.1;
% plot3(beta,TT,abs(F/F(length(TT))),'color','b' )
% hold on
end

end 
figure(h1);     
set(gca,'PlotBoxAspectRatio',[1 1 1],'FontSize',15 )
%xlim([0,4])
%ylim([0,2.5])
ylabel(gca,'T/T_{c0}','FontSize',15)
xlabel(gca,'d_N/\xi_0','FontSize',15)

zlim([0 2.7]) 
view([-70 20])
grid on


figure(h2);     
set(gca,'PlotBoxAspectRatio',[1 1 1],'FontSize',26 )

%lgd=legend(' 0.8','0.2','0')
%hTitle = legendTitle ( lgd, '$d_N/\xi_0=$','interpreter','latex','FontSize',22);
%lgd.Location='best'

lgd=legend(' 0.8','0.2','0')
%hTitle = legendTitle ( lgd, '$d_N/\xi_0=$','interpreter','latex','FontSize',25);
lgd.Location='northwest'
lgd.Orientation='vertical'

ylim([0 5])
%xtick([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
ax = gca;
%ax.XTick = [0 0.2 0.4 0.6 0.8 1];
%ax.XTickLabel = [0  0.5  1];

grid on
grid minor

xlabel('$T/T_{c0}$','interpreter','latex','FontSize',29)
  ylabel('$\delta\alpha/\delta\alpha_N$','interpreter','latex','FontSize',29)
  
%  xlabel('$x/\xi_0$','interpreter','latex','FontSize',29)
% ylabel('$\varepsilon/\Delta_0$','interpreter','latex','FontSize',29)

fname=['GD1dVarDGamma0001']

 fname1=[fname '.png']
 print(gcf,fname1,'-dpng','-r300')
 


% beta=0*x + dT*0.01;
% figure(h1)
% plot3(beta,x,Op,'color','b')
% hold on


%ylim([0, max(DissKin-DissSpec)])