function s=sensitivitaet_msh(mesh,el,ngl)
% SENSITIVITAET
%  function s=sensitivitaet(x,z,rho,el,ngl)
% el: Elektroden in Reihenfolge [A M N B]
% ngl: Anzahl der St�tzstellen f�r Gau�-Legendre
if(nargin<3)
   ngl=6;
end
nel=length(el);
if(nel<2)
   warning('Minimum zwei Elektroden!');
end

xa=el(1);
xm=el(2);
if(nel>2) %at least 3 ele
   xn=el(3);
end
if(nel==4) %4 ele
   xb=el(4);
end


% A M N B
switch nel
 case 4
    sen=Sens.senscalcfriedel_msh(mesh,ngl,xa,xm,xn,xb);
 case 3
    sen=Sens.senscalcfriedel_msh(mesh,ngl,xa,xm,xn);
 case 2
    sen=Sens.senscalcfriedel_msh(mesh,ngl,xa,xm);
end
s=sen;
