
% %          Bus Bus  Voltage Angle   ---Load---- -----Generator-----Static Mvar
% %          No  code Mag.   Degree   MW    Mvar      MW    Mvar    Qmin  Qmax Qc/-Ql
busdata=  [1   1    1.040     0       0      0        0       0       0     0    0;
           2   2    1.025     0       0      0      1.630     0     -0.99  0.99  0;
           3   2    1.025     0       0      0      0.850     0     -0.99  0.99  0;
           4   0    1.000     0       0      0        0       0       0     0    0;
           5   0    1.000     0      1.25   0.5       0       0       0     0    0;
           6   0    1.000     0      0.90   0.3       0       0       0     0    0;
           7   0    1.000     0       0      0        0       0       0     0    0;
           8   0    1.000     0      1.0    0.35      0       0       0     0    0;
           9   0    1.000     0       0      0        0       0       0     0    0];
% 
% %         Line Code                                  
% %      Bus bus      R         X      1/2 B    
% %       nl  nr     p.u.      p.u.     p.u.     T
linedata=[1   4    0.00000   0.05760   0.00000    1;
          2   7    0.00000   0.06250   0.00000    1;
          3   9    0.00000   0.05860   0.00000    1;
          4   5    0.01000   0.08500   0.17600    1;
          4   6    0.01700   0.09200   0.15800    1;
          5   7    0.03200   0.16100   0.30600    1;
          6   9    0.03900   0.17000   0.35800    1;
          7   8    0.00850   0.07200   0.14900    1;
          8   9    0.01190   0.10080   0.20900    1];

%% Data arranged for Linedata in the different vector 
fb=linedata(:,1);tb=linedata(:,2);
r=linedata(:,3);x=linedata(:,4);
b=linedata(:,5);a=linedata(:,6);
z=r+1i*x; 						% Impedance of branch 
y=1./z;b=1i*b;  				% admittance of branch 
nl=length(fb);					% No of branch 
No_of_Bus=max(max(fb),max(tb));		% No of Bus 

%% Formation of YBus matrix 

Y=zeros(No_of_Bus,No_of_Bus);				% Initialize of YBus/k=41 30 coloum 30 baris 
for k=1:nl
    Y(fb(k),tb(k))=Y(fb(k),tb(k))-y(k)/a(k);
    Y(tb(k),fb(k))=Y(fb(k),tb(k));
end
for m=1:No_of_Bus                            %m=41 30 coloum 1 baris
    for n=1:nl
        if fb(n)==m
            Y(m,m)=Y(m,m)+y(n)/a(n)^2+b(n);
        elseif tb(n)==m
            Y(m,m)=Y(m,m)+y(n)+b(n);
        end
    end
end
G=real(Y);B=imag(Y);			% Separation of YBus 
%% Data arranged for Linedata in the different vector
BMva=100;
busNo=busdata(:,1);type=busdata(:,2);V=busdata(:,3);del=busdata(:,4);
Pg=busdata(:,5)/BMva;Qg=busdata(:,6)/BMva;Pl=busdata(:,7)/BMva;
Ql=busdata(:,8)/BMva;Qmin=busdata(:,9)/BMva;Qmax=busdata(:,10)/BMva;
PV_Bus=find(type==2|type==1);PQ_Bus=find(type==3);  	% type1(Slack),type2(PV_Bus Bus),type3(PQ_Bus Bus )
No_of_PQ_Bus=length(PQ_Bus);No_of_PV_Bus=length(PV_Bus);
Active_Power_specified=Pg-Pl;Reactive_Power_specified=Qg-Ql; % Net Power flow through different node 
Iter=1;Tol=1; % Iterantion And tolerance 
%% Newton Raphson Load Flow 
while Tol>1e-5
    P=zeros(No_of_Bus,1);
    Q=zeros(No_of_Bus,1);
    for i=1:No_of_Bus
        for j=1:No_of_Bus
 P(i)=P(i)+V(i)*V(j)*(G(i,j)*cos(del(i)-del(j))+B(i,j)*sin(del(i)-del(j)));
 Q(i)=Q(i)+V(i)*V(j)*(G(i,j)*sin(del(i)-del(j))-B(i,j)*cos(del(i)-del(j)));
        end
    end
	% Verification of limit violation for reactive power 
    if Iter>2 && Iter<=7
        for n=2:No_of_Bus
            if type(n)==2;
                QG=Q(n)+Ql(n);
                if QG > Qmax(n)
                    V(n)=V(n)-0.01;
                elseif QG < Qmin(n)
                    V(n)=V(n)+0.01;
                end
            end
        end
    end
    dPa=Active_Power_specified-P;
    dQa=Reactive_Power_specified-Q;
    dP=dPa(2:No_of_Bus);
    k=1;
    dQ=zeros(No_of_PQ_Bus,1);
    for i=1:No_of_Bus
        if type(i)==3
            dQ(k,1)=dQa(i);
            k=k+1;
        end
    end
    M=[dP;dQ];% delta Matrix 
	%% Formation Fo Jacobian Matrix[J1 J2;J3 J4]
	%% Formation Of J1 
    J1=zeros(No_of_Bus-1,No_of_Bus-1);
    for i=1:No_of_Bus-1
        m=i+1;
        for j=1:No_of_Bus-1;
            n=j+1;
            if m==n
                for n=1:No_of_Bus 
                J1(i,j)=J1(i,j)+V(m)*V(n)*(-G(m,n)*sin(del(m)-del(n))+B(m,n)*cos(del(m)-del(n)));
                end
                J1(i,j)=J1(i,j)-V(m)^2*B(m,m);
            else
                J1(i,j)=V(m)*V(n)*(G(m,n)*sin(del(m)-del(n))-B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
	%% Formation Of J2
    J2=zeros(No_of_Bus-1,No_of_PQ_Bus);
    for i=1:No_of_Bus-1
        m=i+1;
        for j=1:No_of_PQ_Bus
            n=PQ_Bus(j);
            if m==n
                for n=1:No_of_Bus
                    J2(i,j)=J2(i,j)+V(n)*(G(m,n)*cos(del(m)-del(n))+B(m,n)*sin(del(m)-del(n)));
                end
                J2(i,j)=J2(i,j)+V(m)*G(m,m);
            else
                J2(i,j)=V(m)*(G(m,n)*cos(del(m)-del(n))+B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
	%% Formation Of J3
    J3=zeros(No_of_PQ_Bus,No_of_Bus-1);
    for i=1:No_of_PQ_Bus
        m=PQ_Bus(i);
        for j=1:No_of_Bus-1
            n=j+1;
            if m==n
                for n=1:No_of_Bus
                    J3(i,j)=J3(i,j)+V(m)*V(n)*(G(m,n)*cos(del(m)-del(n))+B(m,n)*sin(del(m)-del(n)));
                end
                J3(i,j)=J3(i,j)-V(m)^2*G(m,m);
            else
                J3(i,j)=V(m)*V(n)*(-G(m,n)*cos(del(m)-del(n))-B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
	%% Formation Of J4
    J4=zeros(No_of_PQ_Bus,No_of_PQ_Bus);
    for i=1:No_of_PQ_Bus
        m=PQ_Bus(i);
        for j=1:No_of_PQ_Bus
            n=PQ_Bus(j);
            if m==n
                for n=1:No_of_Bus
                J4(i,j)=J4(i,j)+V(n)*(G(m,n)*sin(del(m)-del(n))-B(m,n)*cos(del(m)-del(n)));
                end
                J4(i,j)=J4(i,j)-V(m)*B(m,m);
            else
                J4(i,j)=V(m)*(G(m,n)*sin(del(m)-del(n))-B(m,n)*cos(del(m)-del(n)));
            end
        end              
    end
   J=[J1 J2;J3 J4]; % Jacobian Matrix 
   X=inv(J)*M;
    dTh=X(1:No_of_Bus-1); % Change in angle 
    dV=X(No_of_Bus:end);	% change in volatge mag 
    del(2:No_of_Bus)=del(2:No_of_Bus)+dTh; % Voltage angle update 
	% voltage mag update 
    k=1;
    for n=2:No_of_Bus
        if type(n)==3
            V(n)=V(n)+dV(k);
            k=k+1;
        end
    end
    Iter=Iter+1;
    Tol=max(abs(M));
end
Q=zeros(No_of_Bus,1);
 for i=1:No_of_Bus
        for j=1:No_of_Bus
            P(i)=P(i)+V(i)*V(j)*(G(i,j)*cos(del(i)-del(j))+B(i,j)*sin(del(i)-del(j)));
            Q(i)=Q(i)+V(i)*V(j)*(G(i,j)*sin(del(i)-del(j))-B(i,j)*cos(del(i)-del(j)));
        end
 end
 for i=1:No_of_Bus
     del(i)=180*del(i)/pi; % Converion radian to degree 
 end
 %% Load Flow Solution 
disp('----------------------------------------');
disp('  Newton Raphson Loadflow Solution    ');
disp('----------------------------------------');
disp(' |Bus |   |Voltage|    |Angle |');
disp(' | No.|   |pu     |    |Degree|');
disp('----------------------------------------');
for m=1:No_of_Bus 
    fprintf(' %3g   ' ,m);
    fprintf(' %8.3f    ' ,V(m));
    fprintf(' %8.3f  ' ,del(m));
     fprintf(' %8.3f  ',Pg(m)*BMva);
    if type(m)==2
    fprintf(' %8.3f  ',(Q(m)+Ql(m))*BMva);
     
    end
    fprintf('\n');
end
disp('----------------------------------------');
fprintf( 'Number Of Iteration %3g \n',Iter)
  