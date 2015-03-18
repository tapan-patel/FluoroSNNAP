function [I,HA,HB] = MutualInformation(a,b)
% Give two signals (a and b), determine the mutual information - very crude
% implementation. a and b are row vectors.
% tau is time delay for the second signal
I = 0;
if(size(a,1)>1)
    a = a';
end
if(size(b,1)>1)
    b = b';
end

% NaN pad to make sizes equal
if(length(a)>length(b))
    b = [b NaN(1,length(a)-length(b))];
end
if(length(b)>length(a))
    a = [a NaN(1,length(b)-length(a))];
end
% bins = max([ceil(length(a)^(1/2)) 10]);
bins = 100:1e3:1e5;
pa = hist(a,bins); pa = pa/sum(pa);
pb = hist(b,bins); pb = pb/sum(pb);
pab = hist3([a',b'],[{bins} {bins}]); pab = pab/sum(pab(:));

HA = pa.*log(pa);
HA = -sum(HA(~isnan(HA)));

HB = pb.*log(pb);
HB = -sum(HB(~isnan(HB)));

for i=1:length(pa)
    for j=1:length(pb)
        tmp=pab(i,j)*log(pab(i,j)/(pa(i)*pb(j)));
        if(~isnan(tmp))
            I = I + tmp;
        end
    end
end