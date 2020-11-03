% %       Bus Bus  Voltage  Angle    ---Load---- -----Generator-----Static Mvar
% %       No  code   Mag.   Degree    MW    Mvar   MW    Mvar   Qmin  Qmax Qc/-Ql
busdata=  [1   1    1.040     0        0     0     0.0     0     0      0    0;
           2   2    1.025     0        0     0     163     0     0      0    0;
           3   2    1.025     0        0     0      85     0     0      0    0;
           4   0    1.000     0       0.0   0.0    0.0     0     0      0    0;
           5   0    1.000     0       125    50    0.0     0     0      0    0;
           6   0    1.000     0        90    30    0.0     0     0      0    0;
           7   0    1.000     0       0.0   0.0    0.0     0     0      0    0;
           8   0    1.000     0       100    35    0.0     0     0      0    0;
           9   0    1.000     0       0.0   0.0    0.0     0     0      0    0];
% 
% %      Line Code                                  
% %      Bus  bus     R        X       1/2 B    
% %       fb  tb     p.u.     p.u.     p.u.       T
linedata=[1   4    0.00000   0.05760  0.00000     1;
          2   7    0.00000   0.06250  0.00000     1;
          3   9    0.00000   0.05860  0.00000     1;
          4   5    0.01000   0.08500  0.17600/2   1;
          4   6    0.01700   0.09200  0.15800/2   1;
          5   7    0.03200   0.16100  0.30600/2   1;
          6   9    0.03900   0.17000  0.35800/2   1;
          7   8    0.00850   0.07200  0.14900/2   1;
          8   9    0.01190   0.10080  0.20900/2   1];
%% Data arranged for Linedata in the different vector 
fb=linedata(:,1);tb=linedata(:,2);
r=linedata(:,3);x=linedata(:,4);
b=linedata(:,5);a=linedata(:,6);
z=r+1i*x; 						% Impedance of branch 
y=1./z;
b=1i*b;  				% admittance of branch 
nl=length(fb);					% No of branch 
No_of_Bus=max(max(fb),max(tb));		% No of Bus 

%% Formation of YBus matrix 

Y=zeros(No_of_Bus,No_of_Bus);				% Initialize of YBus/k=41 30 coloum 30 baris 
for k=1:nl
    Y(fb(k),tb(k))=Y(fb(k),tb(k))-y(k)/a(k);  %mendapatkan hasil tanpa penjumlahan diagonal
    Y(tb(k),fb(k))=Y(fb(k),tb(k));
end
for m=1:No_of_Bus                            %m=41 30 coloum 1 baris
    for n=1:nl
        if fb(n)==m
            Y(m,m)=Y(m,m)+y(n)/a(n)^2;%+b(n); %b gak guna
        elseif tb(n)==m
            Y(m,m)=Y(m,m)+y(n);%+b(n);%b gak guna
        end
    end
end
y=(Y);
AP=real(Y);
RP=imag(Y);
G=abs(Y);
B=angle(Y)/pi*180;% Separation of YBus
%% Data arranged for Linedata in the different vector
BMva=100;
busNo=busdata(:,1);type=busdata(:,2);V=busdata(:,3);del=busdata(:,4);
Pg=busdata(:,5)/BMva;Qg=busdata(:,6)/BMva;Pl=busdata(:,7)/BMva;
Ql=busdata(:,8)/BMva;Qmin=busdata(:,9)/BMva;Qmax=busdata(:,10)/BMva;
PV_Bus=find(type==2|type==1);PQ_Bus=find(type==3);  	% type1(Slack),type2(PV_Bus Bus),type3(PQ_Bus Bus )
No_of_PQ_Bus=length(PQ_Bus);No_of_PV_Bus=length(PV_Bus);
Active_Power=Pg-Pl;
Reactive_Power=Qg-Ql;% Net Power flow through different node

P = zeros(nl,1);
Q = zeros(nl,1);
      % Calculate P and Q
      for i = 1:nl
          for k = 1:nl
              P(i)=P(i)+V(i)*V(k)*G(i,k)*cosd(B(i,k)-del(i)+del(k));
              Q(i)=Q(i)-V(i)*V(k)*G(i,k)*sind(B(i,k)-del(i)+del(k));
          end
      end
Psch=Active_Power-P;
Qsch=Reactive_Power-Q;

%Jacobian 1
J1 = zeros(nl-1,nl-1);
    for i = 1:(nl-1)
        m = i+1;
        for k = 1:(nl-1)
            n = k+1;
            if n == m
                 for n=1:No_of_Bus
                    J1(i,k) = J1(i,k) + V(m)*V(n)*G(m,n)*sind(B(m,n)-del(m)+del(n));
                 end
                    %J1(i,k) = J1(i,k)- V(n)*V(m)*G(n,m)*sind(B(n,m)-del(n)+del(m));
            else
                    J1(i,k) = J1(i,k)-V(n)*V(m)*G(n,m)*sind(B(n,m)-del(n)+del(m));
            end
        end
    end
 J1; 
 yuyu=diag(J1);
for h=1:nl-1
    p=h+1;
    for g=1
        yuyu(h,g)=yuyu(h,g)-V(p)*V(p)*G(p,p)*sind(B(p,p)-del(p)+del(p));
    end
end
 v = yuyu;
 n = size(J1,1);
 J1(1:(n+1):end) = v