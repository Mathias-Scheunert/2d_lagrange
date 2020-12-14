function s = senscalcfriedel_msh(mesh,nw,a,m,n,b)
%SENSCALCFRIEDEL
% 2D sensitivity of a cell within mesh
% Usage: S=senscalc(mesh,nw,a,m[,n[,b]])
% nw quadrature order
% a,m,n,b: vectors of electrode positions

if length(a)~=length(m)
    error('length(a) ~= length(m)');
end
[x,w] = Quad.getQuadratureRule(nw,2); %gausspkt und gewichte in ref dreieck
                                     %(order, dim)
xx = cell(length(mesh.maps),1);
for ii = 1:length(mesh.maps)
    for jj= 1: size(x,1) % xx: Matrix [n(quad order) x 2]
        xx{ii}(jj,:)= mesh.maps{ii}.B*x(jj,:)' + mesh.maps{ii}.b;
    end
end

s = zeros(length(xx), length(a));
for hh = 1:length(a) % loop über alle quellen
    for ii = 1:length(xx) % A : M
        S=integral(a(hh),m(hh),xx{ii}(:,1),xx{ii}(:,2));
         if nargin>=5     % A : M & N
             S=S-integral(a(hh),n(hh),xx{ii}(:,1),xx{ii}(:,2));
             if nargin>=6 % (A : M & N) & (B : M & N)
                 S=S-integral(b(hh),m(hh),xx{ii}(:,1),xx{ii}(:,2))+integral(b(hh),n(hh),xx{ii}(:,1),xx{ii}(:,2));
             end
         end
        % Summation (Quadratur) und Flächenwichtung
        s(ii,hh)= sum(w.*S);%./abs(0.5*mesh.maps{ii}.detB); %0.5*detB = area
    end
end

if nargout==0
    mm=max(abs(S(:)));
    mm=mm/2;
    surf(xx,zz,S);
    view(0,-90)
    shading flat
    colormap('bluewhitered');
    caxis([-mm mm]);
    colorbar
end

function I=integral(xa,xm,X,Z)
% a ... x Koo der Quelle
% m ... x Koo des Empf.
% X ... Vektor, x Koo der Quad.Pkte.
% Z ... Vektor, z Koo der Quad.Pkte.

A = (X-xa).^2+Z.^2;
B = (X-xm).^2+Z.^2;

n=find(B>A);
    DD=A(n);
    A(n)=B(n);
    B(n)=DD;

[K,E]=ellipke(1-B./A);

T2=2./(sqrt(A).*(A-B).^2);
T1=T2./B.*((A+B).*E-2*B.*K);
T2=T2   .*((A+B).*K-2*A.*E);

n=find(abs(A-B) < 1e-10);
    T1(n)= 3*pi ./ (8* sqrt(A(n)).^5);
    T2(n)=   pi ./ (8* sqrt(A(n)).^3);

I = 1./(4*pi^2).*(T1.*((X-xa).*(X-xm)+Z.^2)+T2);
