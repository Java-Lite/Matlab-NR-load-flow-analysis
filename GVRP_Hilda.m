function [TDC]=GVRP_Hilda(FoodPosition,matrix_jarak,Demand,Capacity,v,KPL,Cf,Ce,Cv,Lt,Ult,p,M)
[R urut]= sort(FoodPosition,2,'descend');
Route= [1 urut+1 1];
%Formulating distance matrix===============================================
%matrix_jarak untuk jumlah retail dan DC (matrikx nxn)
%permintaan dibuat n baris 1 kolom
%Lt= loading time untuk retail
%Ult= unloading time untuk retail
d = matrix_jarak;
%==========================================================================
[aa,bb]=size(Route);
for hh=1:aa
    ee=1;
    cc=2;
    while   hh<=aa && cc<bb
            ff=2;
            Routes{hh,ee}(1)=1;
            Load(hh,ee)=0;
            while   Load(hh,ee)<=Capacity && cc<bb
                if  (Load(hh,ee)+Demand(Route(hh,cc)))>Capacity 
                    break
                end     
                Routes{hh,ee}(ff)=Route(hh,cc);
                Load(hh,ee)=Load(hh,ee)+Demand(Route(hh,cc));
                cc=cc+1;
                ff=ff+1;
            end
            Routes{hh,ee}(ff)=1;
            ee=ee+1;
    end
    hh=hh+1;
end
%==========================================================================
NumberOfVehicles=length(Routes);
RouteSets=Routes;
;
for i=1:NumberOfVehicles
    SubRute=flip(RouteSets{i});
    angkut(i,1)=0;
    for j=1:length(SubRute)-1
            LPH (i,j)=v(SubRute(j),SubRute(j+1))./KPL;
            F (i,j)=  LPH (i,j).*( d(SubRute(j),SubRute(j+1))./v(SubRute(j),SubRute(j+1)));
            Vehicle_usage (i,j)=d(SubRute(j),SubRute(j+1))./v(SubRute(j),SubRute(j+1));
        end
    Fuel_Consumption (i,1)= F(i,1);
        for j=2:length(SubRute)-1
            angkut(i,j)=angkut(i,j-1)+Demand(SubRute(1,j),1);
            Fuel_Consumption (i,j)= Fuel_Consumption(i,j-1)+ (F(i,j).*(1+(p.*(angkut(i,j)./M))));
            Load_unload_time(i,j-1)=Lt(SubRute(1,j),1)+Ult(SubRute(1,j),1);
        end
        Fuel_cost(1,i)= max(Fuel_Consumption(i,:)).*Cf;
        Emission_cost (1,i)= max(Fuel_Consumption(i,:)).*Ce;
        Vehicle_Cost(1,i)=(sum(Load_unload_time(i,:))+sum(Vehicle_usage(i,:))).*Cv;
end
TFC=sum(Fuel_cost); %total biaya bahan bakar
TEC=sum(Emission_cost);% total biaya emisi
TVC=sum(Vehicle_Cost); % vehicle cost
TDC=TFC+TEC+TVC;


