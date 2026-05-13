function output = ML_Estimate(u,sec,l,Es)
    staux=[sec zeros(1,l)];
    stEst=zeros(l+1,length(staux));
    stEst(1,:)=staux;
    for i=1:l
     aux=[zeros(1,i) staux];
     aux(length(staux)+1:size(aux,2))=[];
     stEst(i+1,:)=aux;
    end    
    estimation = zeros(1,l);
    uEst=[u zeros(1,size(stEst,2)-length(u))];
    
    for k = 1:l
        estimation(k) = sum(uEst.*conj(stEst(k,:)))/Es;
    end
    output=estimation;
end