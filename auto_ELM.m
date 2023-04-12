function [H,F] = auto_ELM(Y, hN, aF, C, type)

% Original Author of this code is Yimin Yang
% His code can be downloaded at this link https://www.yiminyang.com/code.html
% I have modified his code to spatiotemporal case

iN=size(Y',1);                                          
% hN = number of hidden neurons
% iN = number of input neurons
% aF = activation fuction

for j=1:2
    
    %%%%%%%%%%% Random generate input weights InputWeight (w_i) and biases BiasofHiddenNeurons (b_i) of hidden neurons
    if j==1
        rng(1);
        bf=orth(rand(hN,1));   %random bias of the nuerons bf
        af=rand(hN,iN)*2-1; % random weights of the neurons af
        
        if hN > iN
            af = orth(af);
        else
            af = orth(af')';
        end
        
        YYM_H=af*Y';
        YYM_H= bsxfun(@minus, YYM_H', bf.')';
        %%%%%%%%%%% Or directly inherit the previous parameters
        
        H =mapminmax(YYM_H,-1,1);
%         H= (H'*diag(1./sqrt(sum(H'.*H'))))';
        F =H*inv(H'*H+eye(size(H,2))/C)*Y;
        
    else     
        YYM_H= bsxfun(@minus, (an*Y')', bn.')';
        H =mapminmax(YYM_H,-1,1);
        F= bsxfun(@minus, (an_s*Y)', bn_s.')';
        
    end
    
    
    if type==1
        switch lower(aF)
            case {'sig','sigmoid'}
                H   = 1 ./ (1 + exp(-H));
                F   = 1 ./ (1 + exp(-F));
            case {'sin','sine'}
                H = sin(H);
                F = sin(F);
            case  {'radbas'}
                H = radbas(H);
                F = radbas(F);
            case {'hardlim'}
                H = double(hardlim(H));
                F = double(hardlim(F));
        end
    else
        H = H;
        F = F;
    end
    
    

    if j ==1
        switch lower(aF)
            case {'sig','sigmoid'}
                [tmpG,PS(1)]=mapminmax(Y',0.01,0.99);
                tmpG=(-log((1./tmpG)-1));
                tmpJ=(-log((1./Y)-1));
            case {'sin','sine','radbas','hardlim'}
                [tmpG,PS(1)]=mapminmax(Y',-1,1);
                tmpG=asin(tmpG);
                tmpJ=asin(Y);
        end
    end


    G=real(tmpG);
    an=(eye(size(H,1))/C+H * H') \ H *G';
    BB = (sum(H'*an-G'))/size(G,2);
%     BB= sqrt(mse(YJX-Y4'));
    bn=BB(1);
    
    J=real(tmpJ);
    an_s=(eye(size(F,1))/C+F * F') \ F *J';
    BB2 = (sum(F'*an_s-J'))/size(J,2);
%     BB= sqrt(mse(YJX-Y4'));
    bn_s=BB2(1);
    
end

H = H; %(H'*diag(1./sqrt(sum(H'.*H'))))';
F = F;

