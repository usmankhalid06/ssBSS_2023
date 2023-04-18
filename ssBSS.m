 function [T,S,Err]= ssBSS(Y,Dp,params)
 
K = params.K;
P = params.P;
lam1 = params.lam1; 
lam2 = params.lam2; 
lam3 = params.lam3; 
zeta = params.zeta;
Kp = params.Kp;
nIter = params.nIter;
algo = params.algo;
upd = params.upd;
alpha = params.alpha;

if algo=='PC'
    [F, V, D]= svds(Y,K);
    Xt = V*D'*Y';
    Xs = F'*Y;
    Xt = (Xt'*diag(1./sqrt(sum(Xt'.*Xt'))))'; 
    rng(1,'twister'); T = randn(P,size(Y,1));   T = (T'*diag(1./sqrt(sum(T'.*T'))))';  S = pinv(T')*Y;
else
    [Xt,Xs] = auto_ELM(Y,K,'sine',2^20,0); %1/alpha
    rng(1,'twister'); T = randn(P,size(Y,1));   T = (T'*diag(1./sqrt(sum(T'.*T'))))';  S = pinv(T')*Y;
end


Dp2 = Dp(:,1:Kp);
U = zeros(K,P);
W = zeros(K,P);

for j= 1:nIter
    Told = T; 
    
    T = (Y*S'*inv((S*S') + alpha*speye(size(S,1))))'; T= T'; T = T*diag(1./sqrt(sum(T.*T))); T= T';  
    
    if upd == 'blk'
        A = zeros(size(Dp2,2),P);
        tmp3 = Xt*T'*inv((T*T') + alpha*speye(size(T,1)));
        for i =1:P
            U(:,i) = sign(tmp3(:,i)).*max(0, bsxfun(@minus,abs(tmp3(:,i)),lam1/2));
        end
        T = inv(U'*U+ alpha*speye(size(U,2)))*U'*Xt;
        T= T'; 
        for k = 1:P
            [~,bb]= sort(abs(Dp2'*T(:,k)),'descend');
            ind = bb(1:zeta);
            A(ind,k)= (Dp2(:,ind)'*Dp2(:,ind))\Dp2(:,ind)'*T(:,k);
            A(:,k) = A(:,k)./norm(Dp2*A(:,k));
            T(:,k) = Dp2*A(:,k);
        end
        T = T'; 
    else
        A = zeros(size(Dp2,2),P);
        for i =1:P
            U(:,i) = 0;
            E = Xt-U*T;
            tmp3(:,i)=T(i,:)*E';
            spa = lam1./abs(tmp3(:,i));
            U(:,i) = sign(tmp3(:,i)).*max(0, bsxfun(@minus,abs(tmp3(:,i)),spa/2));
            T(i,:) = U(:,i)'*E;
            [~,bb]= sort(abs(Dp2'*T(i,:)'),'descend');
            ind = bb(1:zeta);
            A(ind,i)= (Dp2(:,ind)'*Dp2(:,ind))\Dp2(:,ind)'*T(i,:)';
            A(:,i) = A(:,i)./norm(Dp2*A(:,i));
            T(i,:) = A(:,i)'*Dp2';
        end
    end

    
    S = inv((T*T') + alpha*speye(size(T,1)))*T*Y;
    
    tmp1 = Xs*S'*inv((S*S') + alpha*speye(size(S,1))); 
    for i =1:P
       W(:,i) = sign(tmp1(:,i)).*max(0, bsxfun(@minus,abs(tmp1(:,i)),lam2/2));
    end
    tmp2 = inv(W'*W+ alpha*speye(size(W,2)))*W'*Xs;
    for i =1:P 
       S(i,:) = sign(tmp2(i,:)).*max(0, bsxfun(@minus,abs(tmp2(i,:)),lam3/2));
       if (length(find(S(i,:)))<1)
           [~,ind]= max(sum(Y-T'*S.^2));
           T(i,:)= (Y(:,ind)./norm(Y(:,ind),2))';
           S(i,:) = 0;
           S(i,:) = sign(T(i,:)*Y).*max(0, bsxfun(@minus,abs(T(i,:)*Y),lam3/2));
       end
    end

    
    Err(j) = (sqrt(trace((T-Told)'*(T-Told)))/sqrt(trace(Told'*Told)));     
    
end


        