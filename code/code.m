%Computer Project - Spring 2014
%Group 20

clear
clc

%Constant and Known
h=30;
k=3;
TempB=50;
TempInfinite=20;
length=0.14;
hight=0.12;

for meshCase=1:4  %for switch 1=2x2 2=2x1 3=1x2 4=1x1
    switch meshCase
        case 1
            numberOfLength=8; %number of mesh index = total mesh number +1
            numberOfHight=7;
        case 2
            numberOfLength=8;
            numberOfHight=13;
        case 3
            numberOfLength=15;
            numberOfHight=7;
        case 4
            numberOfLength=15;
            numberOfHight=13;  
    end

    %get some useful coordinate
    topOfFin=1+(numberOfHight-1)/3;
    bottomOfFin=topOfFin+(numberOfHight-1)/3;
    edgeOfBase=1+(numberOfLength-1)*(4/14);

    dl=length/numberOfLength;   %delta length
    dh=hight/numberOfHight;     %delta hight

    %create a temp array and set all elements as temp of air
    T=TempInfinite*ones(numberOfHight,numberOfLength);

    %initial array Temp by boundary
    T(1:numberOfHight,1)=TempB;

    totalRecursiveTimes=1000;
    for times=1:totalRecursiveTimes

        %index i control Y axis,and index j control X axis
        for i=1:numberOfHight
            for j=2:numberOfLength  %j==1 is boundary(==50) , need not to calculate

                %Encapsulate Boolean
                isMiddle= ((j>1&&j<edgeOfBase)&&(i>1&&i<numberOfHight)) || ((j>=edgeOfBase&&j<numberOfLength)&&(i>topOfFin&&i<bottomOfFin));
                isHorizontalEdge= (i==topOfFin||i==bottomOfFin) && (j>edgeOfBase&&j<numberOfLength);
                isVerticalEdge= (j==edgeOfBase&&((i>1&&i<topOfFin)||(i>bottomOfFin&&i<numberOfHight)))||(j==numberOfLength&&(i>topOfFin&&i<bottomOfFin));
                isCorner= ((i==topOfFin||i==bottomOfFin)&&(j==edgeOfBase));
                isVertexOfFin= ((i==topOfFin||i==bottomOfFin)&&j==numberOfLength);
                isTopOfFin= (i==topOfFin);

                if isMiddle %middle
                T(i,j)=( k*(dl/dh)*(T(i-1,j)+T(i+1,j))+ k*(dh/dl)*(T(i,j-1)+T(i,j+1)))/(2*(k*dh/dl)+2*(k*dl/dh));

                elseif isHorizontalEdge %up and down edge of fin without corner				
                    if isTopOfFin
                    T(i,j)=((0.5*k*dh/dl)*(T(i,j+1)+T(i,j-1))+(k*dl/dh)*T(i+1,j)+h*dl*TempInfinite)/(k*dh/dl+k*dl/dh+h*dl);
                    else
                    T(i,j)=((0.5*k*dh/dl)*(T(i,j+1)+T(i,j-1))+(k*dl/dh)*T(i-1,j)+h*dl*TempInfinite)/(k*dh/dl+k*dl/dh+h*dl);
                    end

                elseif isVerticalEdge %left conduction right convection
                    T(i,j)=((k*dh/dl)*T(i,j-1)+(0.5*k*dl/dh)*(T(i+1,j)+T(i-1,j))+h*dh*TempInfinite)/(k*dh/dl+k*dl/dh+h*dh);

                elseif isCorner  %corner of base and fin
                    if isTopOfFin
                        T(i,j)=( (k*dh/dl)*(T(i,j-1)+0.5*T(i,j+1)) + (k*dl/dh)*(T(i+1,j)+0.5*T(i-1,j) ) + TempInfinite*h*(dh+dl)/2 ) / ( 1.5*k*(dl/dh+dh/dl) + (dl+dh)*h*0.5);
                    else
                        T(i,j)=( (k*dh/dl)*(T(i,j-1)+0.5*T(i,j+1)) + (k*dl/dh)*(T(i-1,j)+0.5*T(i+1,j) ) + TempInfinite*h*(dh+dl)/2 ) / ( 1.5*k*(dl/dh+dh/dl) + (dl+dh)*h*0.5);
                    end

                elseif isVertexOfFin %vertex of fin
                    if isTopOfFin
                        T(i,j)=(0.5*k*(dh/dl)*T(i,j-1)+0.5*k*(dl/dh)*T(i+1,j)+0.5*h*TempInfinite*(dh+dl))/(0.5*k*(dh/dl+dl/dh)+0.5*h*(dh+dl));
                    else
                        T(i,j)=(0.5*k*(dh/dl)*T(i,j-1)+0.5*k*(dl/dh)*T(i-1,j)+0.5*h*TempInfinite*(dh+dl))/(0.5*k*(dh/dl+dl/dh)+0.5*h*(dh+dl));
                    end
                end
            end
        end	

        %because Symmetry ,some mesh point will be the same
        for i=1:numberOfLength
            T(1,i)=T(bottomOfFin,i);
            T(numberOfHight,i)=T(topOfFin,i);
        end
    end

    %for debug
    %T

    %figure
    switch meshCase
        case 1
            figureTitle='2x2';
        case 2
            figureTitle='2x1';
        case 3
            figureTitle='1x2';
        case 4
            figureTitle='1x1';
    end
    subplot(2,2,meshCase);
    contourf(T,30)
    axis([1 numberOfLength 1 numberOfHight]);
    shading flat;
    hold on;
    [a,b]=gradient(T);
    quiver(-a,-b,'k')
    colorbar
    title(figureTitle)
end

%calculate efficiency by eff.m
eff