function [h,logPWTD,PWTD1]=Dopplerplot(X,espace_xx,espace_zz,test,FigFeatures) 
    if nargin<5
        FigFeatures.title=0;
        FigFeatures.print=0;
        FigFeatures.color1 ='b';
        FigFeatures.color2 = '--b';
        FigFeatures.mm=0;
    end  
   
    h1=figure(); 
    [~,~,Nt] =size(X);
    %hold on
    set(gcf,'Color',[1,1,1]);
    colormap hot                                                    ;  
    Amp = 35                                                        ;
    PW = 1/Nt*sum(abs(X).^2,3)                                ;
    PWTD = PW/max(max(PW))                                          ;
    PWTD1=max(PWTD,10^(-Amp/10));
    logPWTD = 10*log10(max(PWTD,10^(-Amp/10)))                      ;
    if test==1
        h = imagesc(logPWTD(:,:),[-Amp,0])      ;
        xlabel('N_X [nb]')                                  ; 
        ylabel('N_Z [nb]')                                  ; 
        set(gca,'XTick', 1:40:espace_xx(end));
        set(gca,'YTick', 1:100:espace_zz(end));       
        if FigFeatures.mm
            hold on
            plot(espace_xx, 228*ones(length(espace_xx),1),FigFeatures.color1,'LineWidth',1.5);
            plot(95*ones(1,length(espace_zz)),espace_zz,FigFeatures.color2,'LineWidth',1.5);
        end
    else
        h = imagesc(espace_xx,espace_zz,logPWTD(:,:),[-Amp,0])      ;                                               ; 
        xlabel('X [mm]','FontSize',14)                                  ; 
        ylabel('Z [mm]','FontSize',14)                                  ; 
        axis ij equal tight     
    end  
    set(gca, 'FontSize', 18, 'fontName','Arial','LineWidth',1.5)    ;    
    if FigFeatures.bar==1
        hh=colorbar                                                      ;
        caxis([-Amp 0])
        %xlabel(hh,'dB')                         ;
    end
    if FigFeatures.title
        newStr = replace(FigFeatures.nomtest,'_', '-');
        title(sprintf('%s', newStr)); %
    end   
    drawnow              ;
    hold off 
   
    if FigFeatures.print
    export_fig(h1, ... % figure handle
                sprintf('%s/%s.pdf', FigFeatures.result_folder, FigFeatures.nomtest),... % name of output file without extension
                '-painters', ...      % renderer
                '-transparent', ...   % renderer
                '-pdf', ...         % file format
                '-r5000');             % resolution in dpi    
    end
end