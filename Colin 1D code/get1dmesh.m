function [coord,conec,nelem,nnode]= ...
    get1dmesh(domain,divide,L)

% compute the number of total node and element

    nelem=sum(divide);    
    nnode=nelem+1;
    coord1=linspace(0.0d0,domain(1),divide(1)+1)-L/2; 
    
    coord2=linspace(domain(1),sum(domain(1:2)),divide(2)+1)-L/2;
    
    coord3=linspace(sum(domain(1:2)),sum(domain(1:3)),divide(3)+1)-L/2;
    
    coord4=linspace(sum(domain(1:3)),sum(domain(1:4)),divide(4)+1)-L/2;
    
    coord5=linspace(sum(domain(1:4)),sum(domain(1:5)),divide(5)+1)-L/2;
    %{
    coord6=linspace(sum(domain(1:5)),sum(domain(1:6)),divide(6)+1)-L/2;
    
    coord7=linspace(sum(domain(1:6)),sum(domain(1:7)),divide(7)+1)-L/2;
    
    coord8=linspace(sum(domain(1:7)),sum(domain(1:8)),divide(8)+1)-L/2;
    
    coord9=linspace(sum(domain(1:8)),sum(domain(1:9)),divide(9)+1)-L/2;
    
    coord=[coord1 coord2(2:divide(2)+1) coord3(2:divide(3)+1) coord4(2:divide(4)+1) coord5(2:divide(5)+1)  coord6(2:divide(6)+1) coord7(2:divide(7)+1) coord8(2:divide(8)+1) coord9(2:divide(9)+1)];
%}
    coord=[coord1 coord2(2:divide(2)+1) coord3(2:divide(3)+1) coord4(2:divide(4)+1) coord5(2:divide(5)+1)];

% element connectivity
conec=[ 1:(nnode-1); 2:nnode ]';