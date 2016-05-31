

function [Ke,f,ubar,time1]=FEMUnitMatrices()

GlobalConst
InitFEMMesh;
tic %Assembly timer

%Assembly K matrix and f vector internal cells
Ke = cell(mCon.m*mCon.nG^2,1);
f=zeros(2*mCon.n,1);

for ic=1:mCon.m                                         % Iterations over the internal cells
    for ip=1:cells(ic).ni                               % Iterations over the cell Gauss points
        B=zeros(3,2*length(cells(ic).nen));             % cells(ic).int(ip).nen = neighboring nodes
        F=zeros(2*length(cells(ic).nen),1);
        en=zeros(1,2*length(cells(ic).nen));
        coord = 2*(cells(ic).int(ip).x-cells(ic).x)./cells(ic).dx;
        [phi,dphidx,dphidy]=FEMShape(coord,length(cells(ic).nen));
        B(1,1:2:end-1)=dphidx;
        B(2,2:2:end)=dphidy;
        B(3,1:2:end-1)=dphidy;
        B(3,2:2:end)=dphidx;
        F(1:2:end-1)=phi*cells(ic).int(ip).cv(1);       % cells(ic).int(ip).cv = body force vector
        F(2:2:end)=phi*cells(ic).int(ip).cv(2);
        en(1:2:end-1)=2*[cells(ic).nen]-1;      % x index of neighboring cells
        en(2:2:end)=2*[cells(ic).nen];          % y index
        Ke{(ic-1)*mCon.nG^2+ip}(en,en)= sparse(B'*pCon.D*B*cells(ic).int(ip).w*cells(ic).J);
        f(en)=f(en)+F*cells(ic).int(ip).w*cells(ic).J;
    end
end


%Assembly G matrix and q vector boundary cells
for ic=1:mCon.mb+mCon.mp
    if bcells(ic).BC==2 || bcells(ic).BC==4
        for ip=1:bcells(ic).ni
            nc = bcells(ic).int(ip).nec;
            T=zeros(2*length(cells(nc).nen),1);
            en=zeros(1,2*length(cells(nc).nen));
            coord = 2*(bcells(ic).int(ip).x-cells(nc).x)./cells(nc).dx;
            phi=FEMShape(coord,length(cells(nc).nen));
            for neni=1:length(cells(nc).nen)
                T(2*neni-1:2*neni)=phi(neni).*bcells(ic).int(ip).bv;
                en(2*neni-1:2*neni)=[2*cells(nc).nen(neni)-1; 2*cells(nc).nen(neni)];
            end
            f(en)=f(en)+T*bcells(ic).int(ip).w*bcells(ic).J;
        end  
    end
end
time1=toc; %Assembly timer
%disp([num2str(time1),' seconds to assemble the matrices'])
tic %Solve timer

% Imposed displacement
ubar = nan(2*mCon.n,1);
for i = 1 : length(bnodes)
    
    indx = 2*bnodes(1,i)-1;
    indy = 2*bnodes(1,i);
    if ~isnan(bnodes(2,i))
        ubar(indx) = bnodes(2,i);
    end
    if ~isnan(bnodes(3,i))
        ubar(indy) = bnodes(3,i);
    end
end

end