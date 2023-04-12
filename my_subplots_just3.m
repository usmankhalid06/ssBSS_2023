function my_subplots_just3(nA,S,v,w,TCcorr,SMcorr,rTC,rSM)
    N = size(rTC{1},1);
    axis off
    set(gca,'Units','normalized','Position',[0 0 1 1]);
    vec = 0.03:1/nA:1; 
%     text(0.0300,0.99, 'Sources','Color','m','FontSize',12)
%     text(0.1950,0.99, ['gICA, \Sigma SMs =' num2str(round(sum(SMcorr(1,:)),2)) ',   \Sigma TCs =' num2str(sum(TCcorr(1,:)))],'Color','m','FontSize',12)
%     text(0.2190,0.99, ['CODL,    \Sigma SMs =' num2str(round(sum(SMcorr(1,:)),2)) ',   \Sigma TCs =' num2str(round(sum(TCcorr(1,:)),2))],'Color','m','FontSize',12)
    text(0.0395,0.99, ['\Sigma cSM_s =' num2str(round(sum(SMcorr(1,:)),2),'%0.2f') ',   \Sigma cTC_s =' num2str(round(sum(TCcorr(1,:)),2),'%0.2f')],'Color','m','FontSize',12)
    text(0.3605,0.99, ['\Sigma cSM_s =' num2str(round(sum(SMcorr(2,:)),2),'%0.2f') ',   \Sigma cTC_s =' num2str(round(sum(TCcorr(2,:)),2),'%0.2f')],'Color','m','FontSize',12)
    text(0.6835,0.99, ['\Sigma cSM_s =' num2str(round(sum(SMcorr(3,:)),2),'%0.2f') ',   \Sigma cTC_s =' num2str(round(sum(TCcorr(3,:)),2),'%0.2f')],'Color','m','FontSize',12)
    text(0.0695,0.02, '(a)','FontSize',9)
    text(0.3905,0.02, '(b)','FontSize',9)
    text(0.7135,0.02, '(c)','FontSize',9)
%     text(vec(7),0.98, 'Proposed2','Color','m','FontSize',12)
    
    for i =1:nA 
        for j=1:S
            ihs = 0.02;  %initial_horizontal_shift
            ivs = 1.30;  %initial_vertical_shift
            shz = 0.085; %subplot_horizontal_size
            svs = S;  %subplot_vertical_size (more the better)
            hs  = 1.40;  %horizontal shift of subplots
            nR  = S+0.7; %No. of rows
            shifter = 0.23;


            %%
            zscore_rxSM = abs(zscore(rSM{i}(j,:)));
            hax=axes();
            imagesc(flipdim(reshape(zscore_rxSM,v,w),1)); %colormap('gray')
            newPos=[hs*(mod(j-1,1)+ihs+shifter*(i-1)),   (1-1/nR)-(1/nR)*(fix((j-1)/1)+ivs-1),   shz,   1/svs];
            set(gca,'outer',newPos),
            set(gca,'XTickLabel','')
            set(gca,'YTickLabel','')
            xlabel(['\gamma',' = ',num2str(round(SMcorr(i,j),2),'%0.2f')],'color','r')

            hax=axes();
            plot(zscore(rTC{i}(:,j))); axis([0 N -3 3]);
            newPos=[hs*(mod(j-1,1)+ihs+shifter*(i-1)+0.04),   (1-1/nR)-(1/nR)*(fix((j-1)/1)+ivs-1),   3*shz,   1/svs];
            set(gca,'outer',newPos),
            set(gca,'XTickLabel','')
            set(gca,'YTickLabel','')
            xlabel(['\gamma',' = ',num2str(round(TCcorr(i,j),2),'%0.2f')],'color','r')
  
        end
    end
