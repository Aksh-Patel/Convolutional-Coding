%%----------------------------BSC----------------------------------------------------------
clearvars; 
close all;
rate=1/2;
snr_db=0:0.5:8;
    snr_linear=10.^(snr_db./10);
sigma2=1./(2*rate.*snr_linear);
    sigma = sqrt(sigma2);
        sigmaset=sigma;
Nsim =2000;
sizeofsigma=size(sigma);
nerr=zeros(1,sizeofsigma(1,2));
perr=qfunc(sqrt(2*rate*snr_linear));
i11=0;
sum1=zeros(1,sizeofsigma(1,2));
% Monte-Carlo Experimental Trials
for sig=sigmaset
    i11=i11+1;
    for ksim = 1:Nsim
        %disp("%----------------------------------------------------------------------------------------");
        k=500;
        info_bits = randi([0 1],[1 k+2]);
        info_bits(1,k+1)=0;
        info_bits(1,k+2)=0;
        %disp(info_bits);
        e1=zeros(1,2*(k+2));
        %-----------------------------------------------Encoding------------------------------------------------
        for i=1:k+2
            if(i==1)
                e1(1,1)=info_bits(1,1);
                e1(1,2)=info_bits(1,1);
            end
            if(i==2)
                e1(1,3)=mod(info_bits(1,1)+info_bits(1,2),2);
                e1(1,4)=info_bits(1,2);
            end
            if(i>=3)
                e1(1,i*2-1)=mod(info_bits(1,i)+info_bits(1,i-1)+info_bits(1,i-2),2);
                e1(1,i*2)=mod(info_bits(1,i)+info_bits(1,i-2),2);
            end
        end
        %disp(e1);
        %--------------------Introducing bsc noise!!!!!!!!!!!!-------------------------------------------------------------------------------
        
        for i=1:2*(k+2)
            error=rand<perr(1,i11);
            if(error)
                if e1(1,i)==0
                    e1(1,i)=1;
                else
                    e1(1,i)=0;
                end
                sum1(1,i11)=sum1(1,i11)+1;
            end
        
        end
        %disp(e1);
        %----------------------------------------Path Matrix Calculation-------------------------------
        pm=zeros(k+2,4);
        %----------------------Calculation of first row of path matrix--
        pm(1,1)=abs(0-e1(1,1))+abs(0-e1(1,2));
        pm(1,2)=abs(1-e1(1,1))+abs(1-e1(1,2));
        pm(1,3)=10000;
        pm(1,4)=10000;
       %----------------------Calculation of rest rows of pth matrix -------
        i1=0;
        for i=2:k+2
            pm(i,1) = min(pm(i-1,1) + e1(1,i1+3)-0 + e1(1,i1+4)-0,pm(i-1,3) + 1-e1(1,i1+3) + 1-e1(1,i1+4));
            pm(i,2) = min(pm(i-1,1) + 1-e1(1,i1+3) + 1-e1(1,i1+4),pm(i-1,3) + e1(1,i1+3)-0 + e1(1,i1+4)-0);
            pm(i,3) = min(pm(i-1,2) + 1-e1(1,i1+3) + e1(1,i1+4)-0,pm(i-1,4) + e1(1,i1+3)-0 + 1-e1(1,i1+4));
            pm(i,4) = min(pm(i-1,2) + e1(1,i1+3)-0 + 1-e1(1,i1+4),pm(i-1,4) + 1-e1(1,i1+3) + e1(1,i1+4)-0);
            i1=i1+2;
        end
        %------------f is the final matrix of calculation of path matrix--------------%
        f=[0 100000 100000 100000;pm];
        %disp(f);
        %----------------------------------------Decoding----------------------------------%
        decoded_bits=zeros(1,k+2);
        m=1;
        i1=2*k+4-1;
        for i=k+2:-1:2
            if m==1
                decoded_bits(1,i)=0;
                q = f(i,1) + e1(1,i1)-0 + e1(1,i1+1)-0;
                r = f(i,3) + 1-e1(1,i1) + 1-e1(1,i1+1);
                if q<=r
                    m=1;
                     i1=i1-2;
                    continue;
                else
                    m=3;
                     i1=i1-2;
                    continue;
                end
            end
            if m==2
                decoded_bits(1,i)=1;
                q = f(i,1) + 1-e1(1,i1) + 1-e1(1,i1+1);
                r = f(i,3) + e1(1,i1)-0 + e1(1,i1+1)-0;
                if q<=r
                    m=1;
                     i1=i1-2;
                    continue;
                else
                    m=3;
                     i1=i1-2;
                    continue;
                end
            end
            if m==3
                decoded_bits(1,i)=0;
                q = f(i,2) + 1-e1(1,i1) + e1(1,i1+1)-0;
                r = f(i,4) + e1(1,i1)-0 + 1-e1(1,i1+1);
                if q<=r
                    m=2;
                     i1=i1-2;
                    continue;
                else
                    m=4;
                     i1=i1-2;
                    continue;
                end
            end
            if m==4
                decoded_bits(1,i)=1;
                q = f(i,2) + e1(1,i1)-0 + 1-e1(1,i1+1);
                r = f(i,4) + 1-e1(1,i1) + e1(1,i1+1)-0;
                if q<=r
                    m=2;
                     i1=i1-2;
                    continue;
                else
                    m=4;
                     i1=i1-2;
                    continue;
                end
            end
        end
        
        if m==1
            decoded_bits(1,1)=0;
        else
            decoded_bits(1,1)=1;
        end
        
        %disp(decoded_bits);
        %--------------------------------Error Calculation------------------------------------ 
        nerr(1,i11)=nerr(1,i11)+sum(xor(info_bits,decoded_bits));
    end
end
tmp=(Nsim*k);
pberrbsc=nerr./tmp;
%disp(sum1);
disp(pberrbsc);
%disp(perr);
figure(1);
semilogy(snr_db,pberrbsc,'o-','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on;


%%

%-----------------------------------------------------GAUSSIAN CHANNEL----------------------------------------------------------------------------
clearvars; 
s=1;
rate=1/2;
snr_db=0:0.5:8;
%snr_db=5;
    snr_linear=10.^(snr_db./10);
sigma2=1./(2*rate.*snr_linear);
    sigma = sqrt(sigma2);
        sigmaset=sigma;
Nsim=2000;
sizeofsigma=size(sigma);
nerr=zeros(1,sizeofsigma(1,2));
i11=0;
% Monte-Carlo Experimental Trials
for sig=sigmaset
    i11=i11+1;
    for ksim = 1:Nsim
        %disp("%----------------------------------------------------------------------------------------");
        k=500;
        info_bits = randi([0 1],[1 k+2]);
        info_bits(1,k+1)=0;
        info_bits(1,k+2)=0;
        %disp(info_bits);
        e1=zeros(1,2*(k+2));
        %-----------------------------------------------Encoding------------------------------------------------
        for i=1:k+2
            if(i==1)
                e1(1,1)=info_bits(1,1);
                e1(1,2)=info_bits(1,1);
                if e1(1,1)==1
                    e1(1,1)=s;
                else
                    e1(1,1)=-s;
                end
                
                if e1(1,2)==1
                    e1(1,2)=s;
                else
                    e1(1,2)=-s;
                end
            end
            if(i==2)
                e1(1,3)=mod(info_bits(1,1)+info_bits(1,2),2);
                e1(1,4)=info_bits(1,2);
                if e1(1,3)==1
                    e1(1,3)=s;
                else
                    e1(1,3)=-s;
                end
                
                if e1(1,4)==1
                    e1(1,4)=s;
                else
                    e1(1,4)=-s;
                end
            end
            if(i>=3)
                e1(1,i*2-1)=mod(info_bits(1,i)+info_bits(1,i-1)+info_bits(1,i-2),2);
                e1(1,i*2)=mod(info_bits(1,i)+info_bits(1,i-2),2);
                if e1(1,i*2-1)==1
                    e1(1,i*2-1)=s;
                else
                    e1(1,i*2-1)=-s;
                end
                
                if e1(1,i*2)==1
                    e1(1,i*2)=s;
                else
                    e1(1,i*2)=-s;
                end
            end
        end
        %disp(e1);
        %--------------------Introducing gaussian noise!!!!!!!!!!!!-------------------------------------------------------------------------------
        for i=1:2*k+4
            n = sigma(1,i11)*randn;
            e1(1,i) = e1(1,i) + n;
        end
        %disp(e1);
        %----------------------------------------Path Matrix Calculation-------------------------------
        pm=zeros(k+2,4);
        
        pm(1,1)=sqrt(abs(-1-e1(1,1))^2+abs(-1-e1(1,2))^2);
        pm(1,2)=sqrt(abs(1-e1(1,1))^2+abs(1-e1(1,2))^2);
        pm(1,3)=10;
        pm(1,4)=10;
       
        i1=0;
        for i=2:k+2
            pm(i,1) = min(pm(i-1,1) + sqrt(abs(e1(1,i1+3)-(-1))^2 + abs(e1(1,i1+4)-(-1))^2) , pm(i-1,3) + sqrt(abs(1-e1(1,i1+3))^2 + abs(1-e1(1,i1+4))^2));
            pm(i,2) = min(pm(i-1,1) + sqrt(abs(1-e1(1,i1+3))^2 + abs(1-e1(1,i1+4))^2) , pm(i-1,3) + sqrt(abs(e1(1,i1+3)-(-1))^2 + abs(e1(1,i1+4)-(-1))^2));
            pm(i,3) = min(pm(i-1,2) + sqrt(abs(1-e1(1,i1+3))^2 + abs(e1(1,i1+4)-(-1))^2) , pm(i-1,4) + sqrt(abs(e1(1,i1+3)-(-1))^2 + abs(1-e1(1,i1+4))^2));
            pm(i,4) = min(pm(i-1,2) + sqrt(abs(e1(1,i1+3)-(-1))^2 + abs(1-e1(1,i1+4))^2) , pm(i-1,4) + sqrt(abs(1-e1(1,i1+3))^2 + abs(e1(1,i1+4)-(-1))^2));
            i1=i1+2;
        end
        
        %------------f is the final matrix of calculation of path matrix--------------%
        
        f=[0 10 10 10;pm];
        %disp(f);
        
        %----------------------------------------Decoding----------------------------------%
        decoded_bits=zeros(1,k+2);
        m=1;
        i1=2*k+4-1;
        for i=k+2:-1:2
            if m==1
                decoded_bits(1,i)=0;
                q = f(i,1) + abs(e1(1,i1)-(-1)) + abs(e1(1,i1+1)-(-1));
                r = f(i,3) + abs(1-e1(1,i1)) + abs(1-e1(1,i1+1));
                if q<=r
                    m=1;
                     i1=i1-2;
                    continue;
                else
                    m=3;
                     i1=i1-2;
                    continue;
                end
            end
            if m==2
                decoded_bits(1,i)=1;
                q = f(i,1) + abs(1-e1(1,i1)) + abs(1-e1(1,i1+1));
                r = f(i,3) + abs(e1(1,i1)-(-1)) + abs(e1(1,i1+1)-(-1));
                if q<=r
                    m=1;
                     i1=i1-2;
                    continue;
                else
                    m=3;
                     i1=i1-2;
                    continue;
                end
            end
            if m==3
                decoded_bits(1,i)=0;
                q = f(i,2) + abs(1-e1(1,i1)) + abs(e1(1,i1+1)-(-1));
                r = f(i,4) + abs(e1(1,i1)-(-1)) + abs(1-e1(1,i1+1));
                if q<=r
                    m=2;
                     i1=i1-2;
                    continue;
                else
                    m=4;
                     i1=i1-2;
                    continue;
                end
            end
            if m==4
                decoded_bits(1,i)=1;
                q = f(i,2) + abs(e1(1,i1)-(-1)) + abs(1-e1(1,i1+1));
                r = f(i,4) + abs(1-e1(1,i1)) + abs(e1(1,i1+1)-(-1));
                if q<=r
                    m=2;
                     i1=i1-2;
                    continue;
                else
                    m=4;
                     i1=i1-2;
                    continue;
                end
            end
        end
        
        if m==1
            decoded_bits(1,1)=0;
        else
            decoded_bits(1,1)=1;
        end
        
        %disp(decoded_bits);
        %--------------------------------Error Calculation------------------------------------ 
        nerr(1,i11)=nerr(1,i11)+sum(xor(info_bits,decoded_bits));
    end
end
tmp=(Nsim*k);
pberrAWGN=nerr./tmp;
disp(pberrAWGN);
%disp(perr);
figure(1);
semilogy(snr_db,pberrAWGN,'^-','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
hold on;

%%
%-----------------------------------------------BEC CHANNEL---------------------------------
clearvars; 
rate=1/2;
snr_db=0:0.5:8;
%snr_db=5;
    snr_linear=10.^(snr_db./10);
sigma2=1./(2*rate.*snr_linear);
    sigma = sqrt(sigma2);
        sigmaset=sigma;
Nsim =2000;
sizeofsigma=size(sigma);
nerr=zeros(1,sizeofsigma(1,2));
perr=qfunc(sqrt(2*rate*snr_linear));
i11=0;
sum1=zeros(1,sizeofsigma(1,2));
% Monte-Carlo Experimental Trials
for sig=sigmaset
    i11=i11+1;
    for ksim = 1:Nsim
        %disp("%----------------------------------------------------------------------------------------");
        k=500;
        info_bits = randi([0 1],[1 k+2]);
        info_bits(1,k+1)=0;
        info_bits(1,k+2)=0;
        %disp(info_bits);
        e1=zeros(1,2*(k+2));
        %-----------------------------------------------Encoding------------------------------------------------
        for i=1:k+2
            if(i==1)
                e1(1,1)=info_bits(1,1);
                e1(1,2)=info_bits(1,1);
            end
            if(i==2)
                e1(1,3)=mod(info_bits(1,1)+info_bits(1,2),2);
                e1(1,4)=info_bits(1,2);
            end
            if(i>=3)
                e1(1,i*2-1)=mod(info_bits(1,i)+info_bits(1,i-1)+info_bits(1,i-2),2);
                e1(1,i*2)=mod(info_bits(1,i)+info_bits(1,i-2),2);
            end
        end
        %disp(e1);
        %--------------------Introducing bec noise!!!!!!!!!!!!-------------------------------------------------------------------------------
        
        for i=1:2*(k+2)
            error=rand<perr(1,i11);
            if(error)
                e1(1,i)=-100;
                sum1(1,i11)=sum1(1,i11)+1;
            end
        
        end
        %disp(e1);
        %----------------------------------------Path Matrix Calculation-------------------------------
        pm=zeros(k+2,4);
        
        pm(1,1)=abs(0-e1(1,1))+abs(0-e1(1,2));
        pm(1,2)=abs(1-e1(1,1))+abs(1-e1(1,2));
        pm(1,3)=10000;
        pm(1,4)=10000;
       
        i1=0;
        for i=2:k+2
             pm(i,1) = min(pm(i-1,1) + abs(e1(1,i1+3)-0) + abs(e1(1,i1+4)-0) , pm(i-1,3) + abs(1-e1(1,i1+3)) + abs(1-e1(1,i1+4)));
             pm(i,2) = min(pm(i-1,1) + abs(1-e1(1,i1+3)) + abs(1-e1(1,i1+4)) , pm(i-1,3) + abs(e1(1,i1+3)-0) + abs(e1(1,i1+4)-0));
             pm(i,3) = min(pm(i-1,2) + abs(1-e1(1,i1+3)) + abs(e1(1,i1+4)-0) , pm(i-1,4) + abs(e1(1,i1+3)-0) + abs(1-e1(1,i1+4)));
             pm(i,4) = min(pm(i-1,2) + abs(e1(1,i1+3)-0) + abs(1-e1(1,i1+4)) , pm(i-1,4) + abs(1-e1(1,i1+3)) + abs(e1(1,i1+4)-0));
             i1=i1+2;        
        end
        
        %------------f is the final matrix of calculation of path matrix--------------%
        
        f=[0 100000 100000 100000;pm];
        %disp(f);
        
        %----------------------------------------Decoding----------------------------------%
        decoded_bits=zeros(1,k+2);
        m=1;
        i1=2*k+4-1;
        for i=k+2:-1:2
            if m==1
                decoded_bits(1,i)=0;
                q = f(i,1) + abs(e1(1,i1)-0) + abs(e1(1,i1+1)-0);
                r = f(i,3) + abs(1-e1(1,i1)) + abs(1-e1(1,i1+1));
                if q<=r
                    m=1;
                     i1=i1-2;
                    continue;
                else
                    m=3;
                     i1=i1-2;
                    continue;
                end
            end
            if m==2
                decoded_bits(1,i)=1;
                q = f(i,1) + abs(1-e1(1,i1)) + abs(1-e1(1,i1+1));
                r = f(i,3) + abs(e1(1,i1)-0) + abs(e1(1,i1+1)-0);
                if q<=r
                    m=1;
                     i1=i1-2;
                    continue;
                else
                    m=3;
                     i1=i1-2;
                    continue;
                end
            end
            if m==3
                decoded_bits(1,i)=0;
                q = f(i,2) + abs(1-e1(1,i1)) + abs(e1(1,i1+1)-0);
                r = f(i,4) + abs(e1(1,i1)-0) + abs(1-e1(1,i1+1));
                if q<=r
                    m=2;
                     i1=i1-2;
                    continue;
                else
                    m=4;
                     i1=i1-2;
                    continue;
                end
            end
            if m==4
                decoded_bits(1,i)=1;
                q = f(i,2) + abs(e1(1,i1)-0) + abs(1-e1(1,i1+1));
                r = f(i,4) + abs(1-e1(1,i1)) + abs(e1(1,i1+1)-0);
                if q<=r
                    m=2;
                     i1=i1-2;
                    continue;
                else
                    m=4;
                     i1=i1-2;
                    continue;
                end
            end
        end
        
        if m==1
            decoded_bits(1,1)=0;
        else
            decoded_bits(1,1)=1;
        end
        
        %disp(decoded_bits);
        %--------------------------------Error Calculation------------------------------------ 
        nerr(1,i11)=nerr(1,i11)+sum(xor(info_bits,decoded_bits));
    end
end
tmp=(Nsim*k);
pberrBEC=nerr./tmp;
%disp(sum1);
disp(pberrBEC);
%disp(perr);
figure(1);
semilogy(snr_db,pberrBEC,'d-','linewidth',2,'color',[0 0.4 0.9],'markerfacecolor',[0 0.4 0.9],'markeredgecolor',[0 0.4 0.9]);
hold on;
grid on;
%legend('BSC','GAUSSIAN NOISE','BEC');
%%xlabel('SNR per Bit in dB'); 
%ylabel('Probability of Bit Error');
%axis([0 8 1e-7 1]); 
%set(gca,'xtick',0:0.5:8);
%saveas(gcf,'projectBasic.jpg','jpg');

%%
%----------------------------------------------BSC-------------------------------------------------
clearvars; 

rate=1/2;
snr_db=0:0.5:8;
%snr_db=5;
    snr_linear=10.^(snr_db./10);
sigma2=1./(2*rate.*snr_linear);
    sigma = sqrt(sigma2);
        sigmaset=sigma;
Nsim = 2000;
sizeofsigma=size(sigma);
nerr=zeros(1,sizeofsigma(1,2));
perr=qfunc(sqrt(2*rate*snr_linear));
i11=0;
%sum1=zeros(1,sizeofsigma(1,2));
% Monte-Carlo Experimental Trials
for sig=sigmaset
    i11=i11+1;
    for ksim = 1:Nsim
        %disp("%----------------------------------------------------------------------------------------");
        k=500;
        info_bits = randi([0 1],[1 k+3]);
        info_bits(1,k+1)=0;
        info_bits(1,k+2)=0;
        info_bits(1,k+3)=0;
        %disp(info_bits);
        e1=zeros(1,2*(k+3));
        %-----------------------------------------------Encoding------------------------------------------------
        for i=1:k+3
            if(i==1)
                e1(1,1)=info_bits(1,1);
                e1(1,2)=info_bits(1,1);
            elseif(i==2)
                e1(1,3)=mod(info_bits(1,1)+info_bits(1,2),2);
                e1(1,4)=mod(info_bits(1,1)+info_bits(1,2),2);
            elseif(i==3)
                e1(1,5)=mod(info_bits(1,2)+info_bits(1,3),2);
                e1(1,6)=mod(info_bits(1,1)+info_bits(1,2)+info_bits(1,3),2);
            elseif(i>=4)
                e1(1,i*2-1)=mod(info_bits(1,i)+info_bits(1,i-1)+info_bits(1,i-3),2);
                e1(1,i*2)=mod(info_bits(1,i)+info_bits(1,i-1)+info_bits(1,i-2),2);
            end
        end
        %disp(e1);
        %--------------------Introducing bsc noise!!!!!!!!!!!!-------------------------------------------------------------------------------
        for i=1:2*(k+3)
            error=rand<perr(1,i11);
            if(error)
                if e1(1,i)==0
                    e1(1,i)=1;
                else
                    e1(1,i)=0;
                end
                %sum1(1,i11)=sum1(1,i11)+1;
            end
        end
        %disp(e1);
        %----------------------------------------Path Matrix Calculation-------------------------------
        pm=zeros(k+4,8);
        %----------------------Calculation of first row of path matrix--
        pm(1,1)=0;
        pm(1,2)=100;
        pm(1,3)=100;
        pm(1,4)=100;
        pm(1,5)=100;
        pm(1,6)=100;
        pm(1,7)=100;
        pm(1,8)=100;
       %----------------------Calculation of rest rows of pth matrix -------
        i1=1;
        for i=2:k+4
            pm(i,1) = min(pm(i-1,1) + e1(1,i1)-0 + e1(1,i1+1)-0,pm(i-1,3) + 1-e1(1,i1) + e1(1,i1+1)-0);
            pm(i,2) = min(pm(i-1,5) + 1-e1(1,i1) + 1-e1(1,i1+1),pm(i-1,7) + e1(1,i1)-0 + 1-e1(1,i1+1));
            pm(i,3) = min(pm(i-1,2) + e1(1,i1)-0 + 1-e1(1,i1+1),pm(i-1,4) + 1-e1(1,i1) + 1-e1(1,i1+1));
            pm(i,4) = min(pm(i-1,6) + 1-e1(1,i1) + e1(1,i1+1)-0,pm(i-1,8) + e1(1,i1)-0 + e1(1,i1+1)-0);
            pm(i,5) = min(pm(i-1,1) + 1-e1(1,i1) + 1-e1(1,i1+1),pm(i-1,3) + e1(1,i1)-0 + 1-e1(1,i1+1));
            pm(i,6) = min(pm(i-1,5) + e1(1,i1)-0 + e1(1,i1+1)-0,pm(i-1,7) + 1-e1(1,i1) + e1(1,i1+1)-0);
            pm(i,7) = min(pm(i-1,2) + 1-e1(1,i1) + e1(1,i1+1)-0,pm(i-1,4) + e1(1,i1)-0 + e1(1,i1+1)-0);
            pm(i,8) = min(pm(i-1,6) + e1(1,i1)-0 + 1-e1(1,i1+1),pm(i-1,8) + 1-e1(1,i1) + 1-e1(1,i1+1));
            i1=i1+2;
        end
        
        
        %disp(pm);
        %-----------------------------------Decoding -----------------------
        decoded_bits=zeros(1,k+3);
        m=1;
        i1=2*(k+3)-1;
        for i=k+4:-1:3
            if m==1
                decoded_bits(1,i-1)=0;
                q = pm(i-1,1) + e1(1,i1)-0 + e1(1,i1+1)-0;
                r = pm(i-1,3) + 1-e1(1,i1) + e1(1,i1+1)-0;
                if q<=r
                    m=1;
                else 
                    m=3;
                end
                i1=i1-2;
                continue;
            end
            if m==2
                decoded_bits(1,i-1)=0;
                q = pm(i-1,5) + 1-e1(1,i1) + 1-e1(1,i1+1);
                r = pm(i-1,7) + e1(1,i1)-0 + 1-e1(1,i1+1);
                if q<=r
                    m=5;
                else
                    m=7;
                end
                i1=i1-2;
                continue;
            end
            if m==3
                decoded_bits(1,i-1)=0;
                q = pm(i-1,2) + e1(1,i1)-0 + 1-e1(1,i1+1);
                r = pm(i-1,4) + 1-e1(1,i1) + 1-e1(1,i1+1);
                if q<=r
                    m=2;
                else
                    m=4;
                end
                i1=i1-2;
                continue;
            end
            if m==4
                decoded_bits(1,i-1)=0;
                q = pm(i-1,6) + 1-e1(1,i1) + e1(1,i1+1)-0;
                r = pm(i-1,8) + e1(1,i1)-0 + e1(1,i1+1)-0;
                if q<=r
                    m=6;
                else
                    m=8;
                end
                i1=i1-2;
                continue;
            end
            if m==5
                decoded_bits(1,i-1)=1;
                q = pm(i-1,1) + 1-e1(1,i1) + 1-e1(1,i1+1);
                r = pm(i-1,3) + e1(1,i1)-0 + 1-e1(1,i1+1);
                if q<=r
                    m=1;
                else
                    m=3;
                end
                i1=i1-2;
                continue;
            end
            if m==6
                decoded_bits(1,i-1)=1;
                q = pm(i-1,5) + e1(1,i1)-0 + e1(1,i1+1)-0;
                r = pm(i-1,7) + 1-e1(1,i1) + e1(1,i1+1)-0;
                if q<=r
                    m=5;
                else
                    m=7;
                end
                i1=i1-2;
                continue;
            end
            if m==7
                decoded_bits(1,i-1)=1;
                q = pm(i-1,2) + 1-e1(1,i1) + e1(1,i1+1)-0;
                r = pm(i-1,4) + e1(1,i1)-0 + e1(1,i1+1)-0;
                if q<=r
                    m=2;
                else
                    m=4;
                end
                i1=i1-2;
                continue;
            end
            if m==8
                decoded_bits(1,i-1)=1;
                q = pm(i-1,6) + e1(1,i1)-0 + 1-e1(1,i1+1);
                r = pm(i-1,8) + 1-e1(1,i1) + 1-e1(1,i1+1);
                if q<=r
                    m=6;
                else
                    m=8;
                end
                i1=i1-2;
                continue;
            end
        end
        if m==1
            decoded_bits(1,1)=0;
        else
            decoded_bits(1,1)=1;
        end
        %disp(decoded_bits);
        nerr(1,i11)=nerr(1,i11)+sum(xor(info_bits,decoded_bits));
        
        
    end
end
tmp=(Nsim*k);
pberrbsc=nerr./tmp;
disp(pberrbsc);
figure(1);
semilogy(snr_db,pberrbsc,'o:','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on;


%%
%----------------------------------------Gaussian----------------------------------------------------
clearvars; 
s=1;
r_decision_boundary=0;
rate=1/2;
snr_db=0:0.5:8;
%snr_db=5;
    snr_linear=10.^(snr_db./10);
sigma2=1./(2*rate.*snr_linear);
    sigma = sqrt(sigma2);
        sigmaset=sigma;
Nsim = 2000;
sizeofsigma=size(sigma);
nerr=zeros(1,sizeofsigma(1,2));
%perr=qfunc(sqrt(2*rate*snr_linear));
i11=0;
%sum1=zeros(1,sizeofsigma(1,2));
% Monte-Carlo Experimental Trials
for sig=sigmaset
    i11=i11+1;
    for ksim = 1:Nsim
        %disp("%----------------------------------------------------------------------------------------");
        k=500;
        info_bits = randi([0 1],[1 k+3]);
        info_bits(1,k+1)=0;
        info_bits(1,k+2)=0;
        info_bits(1,k+3)=0;
        
        %disp(info_bits);
        
        e1=zeros(1,2*(k+3));
        %-----------------------------------------------Encoding------------------------------------------------
        for i=1:k+3
            if(i==1)
                e1(1,1)=info_bits(1,1);
                e1(1,2)=info_bits(1,1);
                if e1(1,1)==1
                    e1(1,1)=s;
                else
                    e1(1,1)=-s;
                end
                
                if e1(1,2)==1
                    e1(1,2)=s;
                else
                    e1(1,2)=-s;
                end
            elseif(i==2)
                e1(1,3)=mod(info_bits(1,1)+info_bits(1,2),2);
                e1(1,4)=mod(info_bits(1,1)+info_bits(1,2),2);
                if e1(1,3)==1
                    e1(1,3)=s;
                else
                    e1(1,3)=-s;
                end
                
                if e1(1,4)==1
                    e1(1,4)=s;
                else
                    e1(1,4)=-s;
                end
            elseif(i==3)
                e1(1,5)=mod(info_bits(1,2)+info_bits(1,3),2);
                e1(1,6)=mod(info_bits(1,1)+info_bits(1,2)+info_bits(1,3),2);
                if e1(1,5)==1
                    e1(1,5)=s;
                else
                    e1(1,5)=-s;
                end
                
                if e1(1,6)==1
                    e1(1,6)=s;
                else
                    e1(1,6)=-s;
                end
            elseif(i>=4)
                e1(1,i*2-1)=mod(info_bits(1,i)+info_bits(1,i-1)+info_bits(1,i-3),2);
                e1(1,i*2)=mod(info_bits(1,i)+info_bits(1,i-1)+info_bits(1,i-2),2);
                if e1(1,i*2-1)==1
                    e1(1,i*2-1)=s;
                else
                    e1(1,i*2-1)=-s;
                end
                
                if e1(1,i*2)==1
                    e1(1,i*2)=s;
                else
                    e1(1,i*2)=-s;
                end
            end
        end
        %disp(e1);
        %--------------------Introducing Gaussian noise!!!!!!!!!!!!-------------------------------------------------------------------------------
        for i=1:2*(k+3)
            n = sigma(1,i11)*randn;
            e1(1,i) = e1(1,i) + n;
        end
        %disp(e1);
        %----------------------------------------Path Matrix Calculation-------------------------------
        pm=zeros(k+4,8);
        %----------------------Calculation of first row of path matrix--
        pm(1,1)=0;
        pm(1,2)=100;
        pm(1,3)=100;
        pm(1,4)=100;
        pm(1,5)=100;
        pm(1,6)=100;
        pm(1,7)=100;
        pm(1,8)=100;
       %----------------------Calculation of rest rows of pth matrix -------
        i1=1;
        for i=2:k+4
            pm(i,1) = min(pm(i-1,1) + (e1(1,i1)-(-1))^2 + (e1(1,i1+1)-(-1))^2  ,  pm(i-1,3) + (1-e1(1,i1))^2 + (e1(1,i1+1)-(-1))^2);
            pm(i,2) = min(pm(i-1,5) + (1-e1(1,i1))^2 + (1-e1(1,i1+1))^2  ,  pm(i-1,7) + (e1(1,i1)-(-1))^2 + (1-e1(1,i1+1))^2);
            pm(i,3) = min(pm(i-1,2) + (e1(1,i1)-(-1))^2 + (1-e1(1,i1+1))^2  ,  pm(i-1,4) + (1-e1(1,i1))^2 + (1-e1(1,i1+1))^2);
            pm(i,4) = min(pm(i-1,6) + (1-e1(1,i1))^2 + (e1(1,i1+1)-(-1))^2  ,  pm(i-1,8) + (e1(1,i1)-(-1))^2 + (e1(1,i1+1)-(-1))^2);
            pm(i,5) = min(pm(i-1,1) + (1-e1(1,i1))^2 + (1-e1(1,i1+1))^2  ,  pm(i-1,3) + (e1(1,i1)-(-1))^2 + (1-e1(1,i1+1))^2);
            pm(i,6) = min(pm(i-1,5) + (e1(1,i1)-(-1))^2 + (e1(1,i1+1)-(-1))^2  ,  pm(i-1,7) + (1-e1(1,i1))^2 + (e1(1,i1+1)-(-1))^2);
            pm(i,7) = min(pm(i-1,2) + (1-e1(1,i1))^2 + (e1(1,i1+1)-(-1))^2  ,  pm(i-1,4) + (e1(1,i1)-(-1))^2 + (e1(1,i1+1)-(-1))^2);
            pm(i,8) = min(pm(i-1,6) + (e1(1,i1)-(-1))^2 + (1-e1(1,i1+1))^2  ,  pm(i-1,8) + (1-e1(1,i1))^2 + (1-e1(1,i1+1))^2);
            i1=i1+2;
        end
        
        %disp(pm);
        %-----------------------------------Decoding -----------------------
        decoded_bits=zeros(1,k+3);
        m=1;
        i1=2*k+6-1;
        for i=k+4:-1:3
            if m==1
                decoded_bits(1,i-1)=0;
                q = pm(i-1,1) + (e1(1,i1)-(-1))^2 + (e1(1,i1+1)-(-1))^2;
                r = pm(i-1,3) + (1-e1(1,i1))^2 + (e1(1,i1+1)-(-1))^2;
                if q<=r
                    m=1;
                else 
                    m=3;
                end
                i1=i1-2;
                continue;
            end
            if m==2
                decoded_bits(1,i-1)=0;
                q = pm(i-1,5) + (1-e1(1,i1))^2 + (1-e1(1,i1+1))^2;
                r = pm(i-1,7) + (e1(1,i1)-(-1))^2 + (1-e1(1,i1+1))^2;
                if q<=r
                    m=5;
                else
                    m=7;
                end
                i1=i1-2;
                continue;
            end
            if m==3
                decoded_bits(1,i-1)=0;
                q = pm(i-1,2) + (e1(1,i1)-(-1))^2 + (1-e1(1,i1+1))^2;
                r = pm(i-1,4) + (1-e1(1,i1))^2 + (1-e1(1,i1+1))^2;
                if q<=r
                    m=2;
                else
                    m=4;
                end
                i1=i1-2;
                continue;
            end
            if m==4
                decoded_bits(1,i-1)=0;
                q = pm(i-1,6) + (1-e1(1,i1))^2 + (e1(1,i1+1)-(-1))^2;
                r = pm(i-1,8) + (e1(1,i1)-(-1))^2 + (e1(1,i1+1)-(-1))^2;
                if q<=r
                    m=6;
                else
                    m=8;
                end
                i1=i1-2;
                continue;
            end
            if m==5
                decoded_bits(1,i-1)=1;
                q = pm(i-1,1) + (1-e1(1,i1))^2 + (1-e1(1,i1+1))^2;
                r = pm(i-1,3) + (e1(1,i1)-(-1))^2 + (1-e1(1,i1+1))^2;
                if q<=r
                    m=1;
                else
                    m=3;
                end
                i1=i1-2;
                continue;
            end
            if m==6
                decoded_bits(1,i-1)=1;
                q = pm(i-1,5) + (e1(1,i1)-(-1))^2 + (e1(1,i1+1)-(-1))^2;
                r = pm(i-1,7) + (1-e1(1,i1))^2 + (e1(1,i1+1)-(-1))^2;
                if q<=r
                    m=5;
                else
                    m=7;
                end
                i1=i1-2;
                continue;
            end
            if m==7
                decoded_bits(1,i-1)=1;
                q = pm(i-1,2) + (1-e1(1,i1))^2 + (e1(1,i1+1)-(-1))^2;
                r = pm(i-1,4) + (e1(1,i1)-(-1))^2 + (e1(1,i1+1)-(-1))^2;
                if q<=r
                    m=2;
                else
                    m=4;
                end
                i1=i1-2;
                continue;
            end
            if m==8
                decoded_bits(1,i-1)=1;
                q = pm(i-1,6) + (e1(1,i1)-(-1))^2 + (1-e1(1,i1+1))^2;
                r = pm(i-1,8) + (1-e1(1,i1))^2 + (1-e1(1,i1+1))^2;
                if q<=r
                    m=6;
                else
                    m=8;
                end
                i1=i1-2;
                continue;
            end
        end
        if m==1
            decoded_bits(1,1)=0;
        else
            decoded_bits(1,1)=1;
        end
        %disp(decoded_bits);
        nerr(1,i11)=nerr(1,i11)+sum(xor(info_bits,decoded_bits));
        
        
    end
end
tmp=(Nsim*k);
pberrAWGN=nerr./tmp;
disp(pberrAWGN);
figure(1);
semilogy(snr_db,pberrAWGN,'o:','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
hold on;


%%
%------------------------------------------BEC------------------------------------
clearvars;
rate=1/2;
snr_db=0:0.5:8;
%snr_db=5;
    snr_linear=10.^(snr_db./10);
sigma2=1./(2*rate.*snr_linear);
    sigma = sqrt(sigma2);
        sigmaset=sigma;
Nsim = 2000;
sizeofsigma=size(sigma);
nerr=zeros(1,sizeofsigma(1,2));
perr=qfunc(sqrt(2*rate*snr_linear));
i11=0;
%sum1=zeros(1,sizeofsigma(1,2));
% Monte-Carlo Experimental Trials
for sig=sigmaset
    i11=i11+1;
    for ksim = 1:Nsim
        %disp("%----------------------------------------------------------------------------------------");
        k=500;
        info_bits = randi([0 1],[1 k+3]);
        info_bits(1,k+1)=0;
        info_bits(1,k+2)=0;
        info_bits(1,k+3)=0;
        %disp(info_bits);
        e1=zeros(1,2*(k+3));
        %-----------------------------------------------Encoding------------------------------------------------
        for i=1:k+3
            if(i==1)
                e1(1,1)=info_bits(1,1);
                e1(1,2)=info_bits(1,1);
            elseif(i==2)
                e1(1,3)=mod(info_bits(1,1)+info_bits(1,2),2);
                e1(1,4)=mod(info_bits(1,1)+info_bits(1,2),2);
            elseif(i==3)
                e1(1,5)=mod(info_bits(1,2)+info_bits(1,3),2);
                e1(1,6)=mod(info_bits(1,1)+info_bits(1,2)+info_bits(1,3),2);
            elseif(i>=4)
                e1(1,i*2-1)=mod(info_bits(1,i)+info_bits(1,i-1)+info_bits(1,i-3),2);
                e1(1,i*2)=mod(info_bits(1,i)+info_bits(1,i-1)+info_bits(1,i-2),2);
            end
        end
        %disp(e1);
        %--------------------Introducing bec noise!!!!!!!!!!!!-------------------------------------------------------------------------------
        for i=1:2*(k+3)
            error=rand<perr(1,i11);
            if(error)
                e1(1,i)=10;
                %sum1(1,i11)=sum1(1,i11)+1;
            end
        end
        %disp(e1);
        %----------------------------------------Path Matrix Calculation-------------------------------
        pm=zeros(k+4,8);
        %----------------------Calculation of first row of path matrix--
        pm(1,1)=0;
        pm(1,2)=100;
        pm(1,3)=100;
        pm(1,4)=100;
        pm(1,5)=100;
        pm(1,6)=100;
        pm(1,7)=100;
        pm(1,8)=100;
       %----------------------Calculation of rest rows of pth matrix -------
        i1=1;
        for i=2:k+4
            pm(i,1) = min(pm(i-1,1) + abs(e1(1,i1)-0) + abs(e1(1,i1+1)-0),pm(i-1,3) + abs(1-e1(1,i1)) + abs(e1(1,i1+1)-0));
            pm(i,2) = min(pm(i-1,5) + abs(1-e1(1,i1)) + abs(1-e1(1,i1+1)),pm(i-1,7) + abs(e1(1,i1)-0) + abs(1-e1(1,i1+1)));
            pm(i,3) = min(pm(i-1,2) + abs(e1(1,i1)-0) + abs(1-e1(1,i1+1)),pm(i-1,4) + abs(1-e1(1,i1)) + abs(1-e1(1,i1+1)));
            pm(i,4) = min(pm(i-1,6) + abs(1-e1(1,i1)) + abs(e1(1,i1+1)-0),pm(i-1,8) + abs(e1(1,i1)-0) + abs(e1(1,i1+1)-0));
            pm(i,5) = min(pm(i-1,1) + abs(1-e1(1,i1)) + abs(1-e1(1,i1+1)),pm(i-1,3) + abs(e1(1,i1)-0) + abs(1-e1(1,i1+1)));
            pm(i,6) = min(pm(i-1,5) + abs(e1(1,i1)-0) + abs(e1(1,i1+1)-0),pm(i-1,7) + abs(1-e1(1,i1)) + abs(e1(1,i1+1)-0));
            pm(i,7) = min(pm(i-1,2) + abs(1-e1(1,i1)) + abs(e1(1,i1+1)-0),pm(i-1,4) + abs(e1(1,i1)-0) + abs(e1(1,i1+1)-0));
            pm(i,8) = min(pm(i-1,6) + abs(e1(1,i1)-0) + abs(1-e1(1,i1+1)),pm(i-1,8) + abs(1-e1(1,i1)) + abs(1-e1(1,i1+1)));
            i1=i1+2;
        end
        
        
        %disp(pm);
        %-----------------------------------Decoding -----------------------
        decoded_bits=zeros(1,k+3);
        m=1;
        i1=2*k+6-1;
        for i=k+4:-1:3
            if m==1
                decoded_bits(1,i-1)=0;
                q = pm(i-1,1) + abs(e1(1,i1)-0) + abs(e1(1,i1+1)-0);
                r = pm(i-1,3) + abs(1-e1(1,i1)) + abs(e1(1,i1+1)-0);
                if q<=r
                    m=1;
                else 
                    m=3;
                end
                i1=i1-2;
                continue;
            end
            if m==2
                decoded_bits(1,i-1)=0;
                q = pm(i-1,5) + abs(1-e1(1,i1)) + abs(1-e1(1,i1+1));
                r = pm(i-1,7) + abs(e1(1,i1)-0) + abs(1-e1(1,i1+1));
                if q<=r
                    m=5;
                else
                    m=7;
                end
                i1=i1-2;
                continue;
            end
            if m==3
                decoded_bits(1,i-1)=0;
                q = pm(i-1,2) + abs(e1(1,i1)-0) + abs(1-e1(1,i1+1));
                r = pm(i-1,4) + abs(1-e1(1,i1)) + abs(1-e1(1,i1+1));
                if q<=r
                    m=2;
                else
                    m=4;
                end
                i1=i1-2;
                continue;
            end
            if m==4
                decoded_bits(1,i-1)=0;
                q = pm(i-1,6) + abs(1-e1(1,i1)) + abs(e1(1,i1+1)-0);
                r = pm(i-1,8) + abs(e1(1,i1)-0) + abs(e1(1,i1+1)-0);
                if q<=r
                    m=6;
                else
                    m=8;
                end
                i1=i1-2;
                continue;
            end
            if m==5
                decoded_bits(1,i-1)=1;
                q = pm(i-1,1) + abs(1-e1(1,i1)) + abs(1-e1(1,i1+1));
                r = pm(i-1,3) + abs(e1(1,i1)-0) + abs(1-e1(1,i1+1));
                if q<=r
                    m=1;
                else
                    m=3;
                end
                i1=i1-2;
                continue;
            end
            if m==6
                decoded_bits(1,i-1)=1;
                q = pm(i-1,5) + abs(e1(1,i1)-0) + abs(e1(1,i1+1)-0);
                r = pm(i-1,7) + abs(1-e1(1,i1)) + abs(e1(1,i1+1)-0);
                if q<=r
                    m=5;
                else
                    m=7;
                end
                i1=i1-2;
                continue;
            end
            if m==7
                decoded_bits(1,i-1)=1;
                q = pm(i-1,2) + abs(1-e1(1,i1)) + abs(e1(1,i1+1)-0);
                r = pm(i-1,4) + abs(e1(1,i1)-0) + abs(e1(1,i1+1)-0);
                if q<=r
                    m=2;
                else
                    m=4;
                end
                i1=i1-2;
                continue;
            end
            if m==8
                decoded_bits(1,i-1)=1;
                q = pm(i-1,6) + abs(e1(1,i1)-0) + abs(1-e1(1,i1+1));
                r = pm(i-1,8) + abs(1-e1(1,i1)) + abs(1-e1(1,i1+1));
                if q<=r
                    m=6;
                else
                    m=8;
                end
                i1=i1-2;
                continue;
            end
        end
        if m==1
            decoded_bits(1,1)=0;
        else
            decoded_bits(1,1)=1;
        end
        %disp(decoded_bits);
        nerr(1,i11)=nerr(1,i11)+sum(xor(info_bits,decoded_bits)); 
    end
end
tmp=(Nsim*k);
pberrbec=nerr./tmp;
disp(pberrbec);
figure(1);
semilogy(snr_db,pberrbec,'o:','linewidth',2,'color',[0 0.4 0.9],'markerfacecolor',[0 0.4 0.9],'markeredgecolor',[0 0.4 0.9]);
hold on;
grid on;
legend('BSC','GAUSSIAN NOISE','BEC','BSC k=4 r=1/2','GAUSSIAN NOISE k=4 r=1/2','BEC k=4 r=1/2');
xlabel('SNR per Bit in dB'); 
ylabel('Probability of Bit Error');
axis([0 8 1e-7 1]); 
set(gca,'xtick',0:0.5:8);
%saveas(gcf,'projectInter1.jpg','jpg');