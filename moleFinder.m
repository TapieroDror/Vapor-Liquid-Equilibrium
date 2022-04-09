function f = moleFinder(n,x,y,z,moles)
%n=[n1l n2l ... ncompl n1v n2v ... ncompv]
comp=length(z);
f=zeros(1,2*comp+comp);
nl=zeros(1,comp);
nv=zeros(1,comp);
for i=1:length(n)
    if i<=comp
        nl(i)=n(i);
    else
        nv(i-comp)=n(i);
    end
end
for i=1:comp-1
    f(i)=x(i)-nl(i)/sum(nl);
end
for i=1:comp-1
    f(i+(comp-1))=y(i)-nv(i)/sum(nv);
end
for i=1:comp
    f(i+(comp-1)*2)=z(i)-(nl(i)+nv(i))/(sum(nl)+sum(nv));
end
for i=1:comp
    f(i+3*comp-2)=moles(i)-nl(i)-nv(i);
end
    