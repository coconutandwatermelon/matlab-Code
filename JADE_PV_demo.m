function [recRMSE,bestP,optValue]=JADE_PV_demo(mType,Max_NFEs,NP)

    addpath('../CauchyFunction');

    %%��ʼ��
    [ UB,LB,Dim ] = Parameter(mType);
    G=0;%���õ���������ǰ����������
    uCR=0.5;%��ʼ���������
    uF=0.5;%��ʼ����������
    p0=0.05;
    top=round(p0*NP);%ÿ�������ŵ�top��
    A=[];%��ʼ�鵵��ȺΪ�ռ�
    t=1;%��¼�鵵��ȺA�еĸ������
    index=mType;%���Ժ�������
    c=0.1;
    
    UB=(repmat((UB),NP,1));
    LB=repmat((LB),NP,1);
    P=(UB-LB).*rand(NP,Dim)+LB;%���������ʼ��Ⱥ����
    
%     if index==1
%         P=load('..\convergence-curve\InitPoPdata\model1pop.mat','Pinit');
%     elseif index==2
%         P=load('..\convergence-curve\InitPoPdata\model2pop.mat','Pinit');
%     elseif index==3
%         P=load('..\convergence-curve\InitPoPdata\model3pop.mat','Pinit');
%     elseif index==4
%         P=load('..\convergence-curve\InitPoPdata\model4pop.mat','Pinit');
%     elseif index==5
%         P=load('..\convergence-curve\InitPoPdata\model5pop.mat','Pinit');
%     end
%     P=double(P.Pinit);
    
    for i=1:NP
        fitnessP(i)=PV_TestFunction(index,P(i,:));%������Ⱥ������Ӧֵ
    end
    NFEs=NP;
    
    [fitnessBestP,indexBestP]=min(fitnessP);
     bestP=P(indexBestP,:);
    recRMSE(1:NP)=fitnessP;
    
    %%�����ѭ��
    while NFEs<Max_NFEs
        
        Scr=[];%��ʼ�ɹ��μӱ���ĸ���Ľ������Ϊ�ռ�
        Sf=[];%��ʼ�ɹ��μӱ���ĸ������������Ϊ�ռ�
        n0=1;%��¼Scr��Sf�е�Ԫ�ظ���
        
        %���ݸ�����Ӧֵ�������򣬵õ���õ�ǰtop������
        [~,indexSortP]=sort(fitnessP);
        
        for i=1:NP
            
            %�Ե�G����ÿ��������㽻����ʺ���������
            CR(i)=normrnd(uCR,0.1);%������̬�ֲ�
            F(i)=cauchyrnd(uF,0.1);%���ӿ����ֲ�

            if(CR(i)>1)
               CR(i)=1;
            elseif CR(i)<0
               CR(i)=0; 
            end
            while (F(i)<=0)
                F(i)=cauchyrnd(uF,0.1);
            end
            if (F(i)>1)
                F(i)=1;
            end
            
        end

        %%�������+�������
        for i=1:NP   
            
            %%�������
            for j=1:top
               bestTopP(j,:)=P(indexSortP(j),:); 
            end
            %��top�����������ѡ��һ����ΪXpbest
            k0=randperm(top,1);
            Xpbest=bestTopP(k0,:);
            %�ӵ�ǰ��ȺP�����ѡ��P1
            k1=randi([1,NP]);
            P1=P(k1,:);
            while (k1==i)
                k1=randi([1,NP]);
                P1=P(k1,:); 
            end
            %��P��A�����ѡ��P2
            PandA=[P;A];
            [num,~]=size(PandA);
            k2=randi([1,num]);
            P2=PandA(k2,:);
            while (k2==i||k2==k1)
                k2=randi([1,num]);
                P2=PandA(k2,:); 
            end
            V(i,:)=P(i,:)+F(i).*(Xpbest-P(i,:))+F(i).*(P1-P2);   
       
             %%�������
            jrand=randi([1,Dim]); 
            for j=1:Dim
                k3=rand;
                if(k3<=CR(i)||j==jrand)
                    U(i,j)=V(i,j);
                else
                    U(i,j)=P(i,j);
                end
            end
            
        end
        
        
        %%�߽紦��
        for i=1:NP
           for j=1:Dim
              if (U(i,j)>UB(i,j)||U(i,j)<LB(i,j))
                 U(i,j)=(UB(i,j)-LB(i,j))*rand+LB(i,j); 
              end
           end
           fitnessU(i)=PV_TestFunction(index,U(i,:));
           NFEs=NFEs+1;

        %%ѡ�����

            if(fitnessU(i)<fitnessP(i))
                A(t,:)=P(i,:);
                P(i,:)=U(i,:);
                fitnessP(i)=fitnessU(i);
                Scr(n0)=CR(i);
                Sf(n0)=F(i);
                t=t+1;
                n0=n0+1;
                if(fitnessU(i)<fitnessBestP)
                   fitnessBestP=fitnessU(i);
                   bestP=U(i,:);
                end
            end
            
            recRMSE(NFEs)=fitnessP(i);
        end

        %�жϹ鵵��ȺA�Ĺ�ģ�Ƿ���NP֮�ڣ������ڣ�������Ƴ�����ʹ���ģ����NP
        [tA,~]=size(A);
        if tA>NP
            nRem=tA-NP;
            k4=randperm(tA,nRem);
            A(k4,:)=[]; 
            [tA,~]=size(A);
            t=tA+1;
        end

        %����Ӧ����������uCR��uF
        [~,ab]=size(Scr);
        if ab~=0
            newSf=(sum(Sf.^2))/(sum(Sf));
            uCR=(1-c)*uCR+c.*mean(Scr);
            uF=(1-c)*uF+c.*newSf;
        end
    

        G=G+1;
        
    end
    
    optValue=fitnessBestP;
    
end
