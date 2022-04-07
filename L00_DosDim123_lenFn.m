%function [nums]=L00_DosDim123_lenFn(L,ks, Lnums)
clc
close all
clear all

L=3;
L=2;
L=8;
ks=0.7;
Lnums=[1:55];

nums=[];
% clc
% clear all
% close all
disp('max_min DoS  filtering...')

%-------------------------------------------------------------------------
%generate the stick templates

%assign the radius value to L , and the length of 
%stick will be 2*L+1
%L=3; L=4; L=5;

Len_stick=2*L+1;
    %generate the kernel in the South_North direction
    SKxy=zeros(2*L+1, 2, 2*L+1);
    disp('the stick kernels in the y direction')
    for i=-L: 1: L
        %generate the offset array
        A=[];
        for j=-L:1:L
            A=[A  [j, round(j/L*i)]' ];            
        end   
        %store the stick offset or integer indices in a 3D array
        SKxy(i+L+1,:,:)=A;  
    end
    %generate the kernel in the West_East direction
    SKyx=zeros(2*L+1, 2, 2*L+1);
    disp('the stick kernels in the x direction')
    for i=-L: 1: L
        %generate the offset array
        A=[];
        for j=-L:1:L
            A=[A  [ round(j/L*i),j]' ];            
        end       
        %store the stick offset or integer indices in a 3D array
        SKyx(i+L+1,:,:)=A;  
    end

%-------------------------------------------------------------------------
%Set the filtering parameters

%A space parameter for the definition of derivatives
step=3; 

% %the coefficient of stad/mean  
% ks=0.3;
% ks=2.0;

% %Select the image number 
%  Lnums=[34:-1:32];


 
 %pre Paths
 preP_org='..\orgImage\';
 if ks<1
     preP_dos=['..\Len17DosK0' num2str(ks*10) '\'];
     pre1=['..\Len17DosK0' num2str(ks*10) '\Dim1\'];
     pre2=['..\Len17DosK0' num2str(ks*10) '\Dim2\'];
     pre3=['..\Len17DosK0' num2str(ks*10) '\Dim3\'];
 else     
     preP_dos=['..\Len17DosK' num2str(ks*10) '\'];
     pre1=['..\Len17DosK' num2str(ks*10) '\Dim1\'];
     pre2=['..\Len17DosK' num2str(ks*10) '\Dim2\'];
     pre3=['..\Len17DosK' num2str(ks*10) '\Dim3\'];
 end
     
 [SUCCESS,MESSAGE,MESSAGEID] = mkdir(pre1);
 [SUCCESS,MESSAGE,MESSAGEID] = mkdir(pre2);
 [SUCCESS,MESSAGE,MESSAGEID] = mkdir(pre3);
 
%%
%the outer Loop (Image number)

for filenum=1:1:length(Lnums)
         
    %%
      % the original image filename  
      nf=Lnums(filenum);
      if nf<10,
         fn1=[ preP_org 'LOLA11-0' num2str(nf) '.mhd'];
      else
        fn1=[ preP_org 'LOLA11-' num2str(nf) '.mhd'];
      end  
      disp('Reading: ')
      disp(fn1)
     %determine the existence of file
       FID = fopen(fn1);
      if FID==-1,
         disp('File doesnot exist!')
          continue;
      else
          fclose(FID);
      end
      %read the original image
      vol1= mha_read_volumeZ(fn1);     
      
      
      %%
      %=============================
      %Filtering in the Dim1 direction     
      
      %cast to double type
      I=single(vol1.pixelData);
      disp(size(I))             

     
     %%
      %max-min iteration
  
    for iter=1:2

      for dim=1 : 1 : 1
        %----------------------------------------------------------------------
        %The second level inner loop (Image slice)
        for Indz=1 : size(I,dim)
            disp(['dim=' num2str(dim) ' Slice=' num2str(Indz)])

            %Extract a 2D slice from the volume image
            switch dim
                case 1
                     I2D=I(Indz,:,:);   
                case 2
                     I2D=I(:,Indz,:);   
                case 3
                     I2D=I(:,:,Indz);  
            end                        
            I2D=squeeze(I2D); %Squeeze to 2D      
            Iorg=I2D; % Save as the original 2D image        
            [M, N]=size(I2D); %size of 2D slice        
            [X0,Y0]=meshgrid(1 :M,1:N);%grid coordinates

            %To save the results of different combinations        
            Imax_max=zeros(size(I2D));          

            %------------------------------------------------------------------
            %The third level inner loop (Stick orientation)

            %For sticks in the North_south directions
            for t1=1:(2*L+1)
                %The temporary variables 
                Imean=zeros(size(X0));  %for the mean along a stick        
                Istd=zeros(size(X0));   %for the standard deviation(STD)

                %Calculate the mean and STD
                for t2=1:(2*L+1)
                    %X offsets
                    X=X0+SKxy(t1,1,t2);
                    %Limit inside the image boundary
                    X(find(X<1))=1;
                    X(find(X>M))=M;

                    %Y offsets 
                    Y=Y0+SKxy(t1,2,t2);
                    %Limit inside the image boundary
                    Y(find(Y<1))=1;
                    Y(find(Y>N))=N;

                    %the current linear indices
                    inds=X+(Y-1)*M;    

                    %the sum of mean and STD
                    Imean=Imean+I2D(inds);      
                    Istd=Istd+I2D(inds).^2;   
                end
                %the mean value
                Imean=Imean/(2*L+1); 
                %the STD
                Istd=Istd/(2*L+1);      
                Istd=(Istd-Imean.^2).^0.5;           

                %calculate the derivative of stick (DoS) 

                %the linear indices of center points
                [Mc,Nc]=size(X0);
                [XC,YC]=meshgrid(1:Mc, 1:Nc);
                indsCenter=XC+(YC-1)*Mc;   

                %the linear indices of left points
                XL=XC-step;
                XL(find(XL<1))=1;
                indsLeft=XL+(YC-1)*Mc;             
                %the left derivatives
                IderivL=( Imean( indsCenter)-Imean( indsLeft));

                %the linear indices of right points
                XR=XC+step;
                XR(find(XR>Mc))=Mc;
                indsRight=XR+(YC-1)*Mc; 
                %the right derivative
                IderivR=(Imean( indsCenter)-Imean( indsRight));    

                %the nonlinear derivatives
                switch iter
                    case 1
                        %IderivES_min=min(IderivL,IderivR); 
                        IderivES_max=max(IderivL,IderivR);
                        %IderivES_max=(IderivL+IderivR)/2;   
                    case 2
                        %IderivES_min=(IderivL+IderivR)/2;   
                        %IderivES_max=(IderivL+IderivR)/2;   
                        IderivES_max=min(IderivL,IderivR);
                    case 3
                        %IderivES_min=(IderivL+IderivR)/2;   
                        %IderivES_max=(IderivL+IderivR)/2;        
                        IderivES_max=max(IderivL,IderivR);
                end               


                %Combine with the std               
                Ideriv_max=IderivES_max-ks*Istd( indsCenter);

                %Integration of responses in varying orientations                
                Imax_max=max(real(Imax_max), real(Ideriv_max));    
            end

            %For sticks in the East_west directions                   
            for t3=1:(2*L+1)
                %for the mean along a stick
                Imean=zeros(size(X0));  
                %for the standard deviation
                Istd=zeros(size(X0));  

                for t4=1:(2*L+1)
                    X=X0+SKyx(t3,1,t4);
                    X(find(X<1))=1;
                    X(find(X>M))=M;

                    Y=Y0+SKyx(t3,2,t4);
                    Y(find(Y<1))=1;
                    Y(find(Y>N))=N;
                    inds=X+(Y-1)*M;        
                    Imean=Imean+I2D(inds);      
                    Istd=Istd+I2D(inds).^2;   
                end

                Imean=Imean/(2*L+1);      

                Istd=Istd/(2*L+1);      
                Istd=(Istd-Imean.^2).^0.5;

                %the linear indices of center points
                [Mc,Nc]=size(X0);
                [XC,YC]=meshgrid(1:Mc, 1:Nc);
                indsCenter=XC+(YC-1)*Mc;                 

                %the derivative in the north direction
                YU=YC-step;
                YU(find(YU<1))=1;
                indsUp=XC+(YU-1)*Mc; 
                %up derivative
                IderivU=(Imean( indsCenter)-Imean( indsUp));

                %the derivative in the south direction
                YD=YC+step;
                YD(find(YD>Nc))=Nc;
                indsDown=XC+(YD-1)*Mc; 
                %up derivative
                IderivD=(Imean( indsCenter)-Imean(  indsDown));

                %the nonlinear derivative
                %IderivNS=min( IderivU, IderivD);             
                switch iter
                    case 1
                         %IderivNS_min=min( IderivU, IderivD); 
                         IderivNS_max=max( IderivU, IderivD); 
                         %IderivNS_max=(IderivU+IderivD)/2; 
                    case 2
                         %IderivNS_min=(IderivU+IderivD)/2; 
                         %IderivNS_max=(IderivU+IderivD)/2; 
                         IderivNS_max=min( IderivU, IderivD); 
                    case 3
                         %IderivNS_min=(IderivU+IderivD)/2; 
                         %IderivNS_max=(IderivU+IderivD)/2;    
                         IderivNS_max=max( IderivU, IderivD); 
                end                

                %Combine with the std            
                Ideriv_max=IderivNS_max-ks*Istd( indsCenter);

                %Integration of responses in varying orientations
                %Imin_max=max(real(Imin_max), real(Ideriv_min));    
                Imax_max=max(real(Imax_max), real(Ideriv_max));    
            end

            %------------------------------------------------------------------
            %Show the original and filtered images for debugging
            %figure(1)
            %subplot(1,2,1),  imshow(real(Iorg),[]), title('Original Image')  
            %Imax=Imax-min(min(Imax));
            %subplot(1,2,2),  imshow(real(Imax),[]), title('Filtered Image')  
            %pause

            %------------------------------------------------------------------
            %Save the 2D results as a 3D output
            switch dim
                case 3
                    %I3dMin(:,:,Indz)=single(Imin_max);    
                    %I3dMax(:,:,Indz)=single(Imax_max); 
                    I(:,:,Indz)=single(Imax_max); 
               case 2                
                    %I3dMin(:,Indz,:)=single(Imin_max);    
                    %I3dMax(:,Indz,:)=single(Imax_max); 
                    I(:,Indz,:)=single(Imax_max); 
                case 1                
                    %I3dMin(Indz,:,:)=single(Imin_max);    
                    %I3dMax(Indz,:,:)=single(Imax_max); 
                    I(Indz,:,:)=single(Imax_max); 
            end              

        end    
        I=single(I); %Cast to a single float type
      end


    end


    %%       
    vol2=vol1;
    vol2.pixelData=single(I);
    vol2.metaData.ElementType='MET_FLOAT';  
    vol2.metaData.CompressedData='True'; 

    if nf<10,
      fn3=[preP_dos '\Dim1\L0'  num2str(nf) 'maxMinDosDim1_L' num2str(Len_stick) 'ks' num2str(10*ks) '.mhd'];
      vol2.metaData.ElementDataFile = ['L0'  num2str(nf) 'maxMinDosDim1_L' num2str(Len_stick) 'ks' num2str(10*ks) '.zraw'];     
    else
      fn3=[preP_dos '\Dim1\L'  num2str(nf) 'maxMinDosDim1_L' num2str(Len_stick) 'ks' num2str(10*ks) '.mhd'];
      vol2.metaData.ElementDataFile = ['L'  num2str(nf) 'maxMinDosDim1_L' num2str(Len_stick) 'ks' num2str(10*ks) '.zraw'];      
    end      
    mha_write_volumeZ(fn3, vol2,1);
    disp('Write: ')     
    disp(fn3)
    
  %%
      %=============================
      %Filtering in the Dim2 direction     
      
     
      %cast to double type
      I=single(vol1.pixelData);
      disp(size(I))             

     
     %%
      %max-min iteration
  
    for iter=1:2

      for dim=2 : 1 : 2
        %----------------------------------------------------------------------
        %The second level inner loop (Image slice)
        for Indz=1 : size(I,dim)
            disp(['dim=' num2str(dim) ' Slice=' num2str(Indz)])

            %Extract a 2D slice from the volume image
            switch dim
                case 1
                     I2D=I(Indz,:,:);   
                case 2
                     I2D=I(:,Indz,:);   
                case 3
                     I2D=I(:,:,Indz);  
            end                        
            I2D=squeeze(I2D); %Squeeze to 2D      
            Iorg=I2D; % Save as the original 2D image        
            [M, N]=size(I2D); %size of 2D slice        
            [X0,Y0]=meshgrid(1 :M,1:N);%grid coordinates

            %To save the results of different combinations        
            Imax_max=zeros(size(I2D));          

            %------------------------------------------------------------------
            %The third level inner loop (Stick orientation)

            %For sticks in the North_south directions
            for t1=1:(2*L+1)
                %The temporary variables 
                Imean=zeros(size(X0));  %for the mean along a stick        
                Istd=zeros(size(X0));   %for the standard deviation(STD)

                %Calculate the mean and STD
                for t2=1:(2*L+1)
                    %X offsets
                    X=X0+SKxy(t1,1,t2);
                    %Limit inside the image boundary
                    X(find(X<1))=1;
                    X(find(X>M))=M;

                    %Y offsets 
                    Y=Y0+SKxy(t1,2,t2);
                    %Limit inside the image boundary
                    Y(find(Y<1))=1;
                    Y(find(Y>N))=N;

                    %the current linear indices
                    inds=X+(Y-1)*M;    

                    %the sum of mean and STD
                    Imean=Imean+I2D(inds);      
                    Istd=Istd+I2D(inds).^2;   
                end
                %the mean value
                Imean=Imean/(2*L+1); 
                %the STD
                Istd=Istd/(2*L+1);      
                Istd=(Istd-Imean.^2).^0.5;           

                %calculate the derivative of stick (DoS) 

                %the linear indices of center points
                [Mc,Nc]=size(X0);
                [XC,YC]=meshgrid(1:Mc, 1:Nc);
                indsCenter=XC+(YC-1)*Mc;   

                %the linear indices of left points
                XL=XC-step;
                XL(find(XL<1))=1;
                indsLeft=XL+(YC-1)*Mc;             
                %the left derivatives
                IderivL=( Imean( indsCenter)-Imean( indsLeft));

                %the linear indices of right points
                XR=XC+step;
                XR(find(XR>Mc))=Mc;
                indsRight=XR+(YC-1)*Mc; 
                %the right derivative
                IderivR=(Imean( indsCenter)-Imean( indsRight));    

                %the nonlinear derivatives
                switch iter
                    case 1
                        %IderivES_min=min(IderivL,IderivR); 
                        IderivES_max=max(IderivL,IderivR);
                        %IderivES_max=(IderivL+IderivR)/2;   
                    case 2
                        %IderivES_min=(IderivL+IderivR)/2;   
                        %IderivES_max=(IderivL+IderivR)/2;   
                        IderivES_max=min(IderivL,IderivR);
                    case 3
                        %IderivES_min=(IderivL+IderivR)/2;   
                        %IderivES_max=(IderivL+IderivR)/2;        
                        IderivES_max=max(IderivL,IderivR);
                end               


                %Combine with the std               
                Ideriv_max=IderivES_max-ks*Istd( indsCenter);

                %Integration of responses in varying orientations                
                Imax_max=max(real(Imax_max), real(Ideriv_max));    
            end

            %For sticks in the East_west directions                   
            for t3=1:(2*L+1)
                %for the mean along a stick
                Imean=zeros(size(X0));  
                %for the standard deviation
                Istd=zeros(size(X0));  

                for t4=1:(2*L+1)
                    X=X0+SKyx(t3,1,t4);
                    X(find(X<1))=1;
                    X(find(X>M))=M;

                    Y=Y0+SKyx(t3,2,t4);
                    Y(find(Y<1))=1;
                    Y(find(Y>N))=N;
                    inds=X+(Y-1)*M;        
                    Imean=Imean+I2D(inds);      
                    Istd=Istd+I2D(inds).^2;   
                end

                Imean=Imean/(2*L+1);      

                Istd=Istd/(2*L+1);      
                Istd=(Istd-Imean.^2).^0.5;

                %the linear indices of center points
                [Mc,Nc]=size(X0);
                [XC,YC]=meshgrid(1:Mc, 1:Nc);
                indsCenter=XC+(YC-1)*Mc;                 

                %the derivative in the north direction
                YU=YC-step;
                YU(find(YU<1))=1;
                indsUp=XC+(YU-1)*Mc; 
                %up derivative
                IderivU=(Imean( indsCenter)-Imean( indsUp));

                %the derivative in the south direction
                YD=YC+step;
                YD(find(YD>Nc))=Nc;
                indsDown=XC+(YD-1)*Mc; 
                %up derivative
                IderivD=(Imean( indsCenter)-Imean(  indsDown));

                %the nonlinear derivative
                %IderivNS=min( IderivU, IderivD);             
                switch iter
                    case 1
                         %IderivNS_min=min( IderivU, IderivD); 
                         IderivNS_max=max( IderivU, IderivD); 
                         %IderivNS_max=(IderivU+IderivD)/2; 
                    case 2
                         %IderivNS_min=(IderivU+IderivD)/2; 
                         %IderivNS_max=(IderivU+IderivD)/2; 
                         IderivNS_max=min( IderivU, IderivD); 
                    case 3
                         %IderivNS_min=(IderivU+IderivD)/2; 
                         %IderivNS_max=(IderivU+IderivD)/2;    
                         IderivNS_max=max( IderivU, IderivD); 
                end                

                %Combine with the std            
                Ideriv_max=IderivNS_max-ks*Istd( indsCenter);

                %Integration of responses in varying orientations
                %Imin_max=max(real(Imin_max), real(Ideriv_min));    
                Imax_max=max(real(Imax_max), real(Ideriv_max));    
            end

            %------------------------------------------------------------------
            %Show the original and filtered images for debugging
            %figure(1)
            %subplot(1,2,1),  imshow(real(Iorg),[]), title('Original Image')  
            %Imax=Imax-min(min(Imax));
            %subplot(1,2,2),  imshow(real(Imax),[]), title('Filtered Image')  
            %pause

            %------------------------------------------------------------------
            %Save the 2D results as a 3D output
            switch dim
                case 3
                    %I3dMin(:,:,Indz)=single(Imin_max);    
                    %I3dMax(:,:,Indz)=single(Imax_max); 
                    I(:,:,Indz)=single(Imax_max); 
               case 2                
                    %I3dMin(:,Indz,:)=single(Imin_max);    
                    %I3dMax(:,Indz,:)=single(Imax_max); 
                    I(:,Indz,:)=single(Imax_max); 
                case 1                
                    %I3dMin(Indz,:,:)=single(Imin_max);    
                    %I3dMax(Indz,:,:)=single(Imax_max); 
                    I(Indz,:,:)=single(Imax_max); 
            end              

        end    
        I=single(I); %Cast to a single float type
      end


    end


    %%       
    vol2=vol1;
    vol2.pixelData=single(I);
    vol2.metaData.ElementType='MET_FLOAT';  
    vol2.metaData.CompressedData='True'; 

     if nf<10,
      fn3=[preP_dos '\Dim2\L0'  num2str(nf) 'maxMinDosDim2_L' num2str(Len_stick) 'ks' num2str(10*ks) '.mhd'];
      vol2.metaData.ElementDataFile = ['L0'  num2str(nf) 'maxMinDosDim2_L' num2str(Len_stick) 'ks' num2str(10*ks) '.zraw'];     
    else
      fn3=[preP_dos '\Dim2\L'  num2str(nf) 'maxMinDosDim2_L' num2str(Len_stick) 'ks' num2str(10*ks) '.mhd'];
      vol2.metaData.ElementDataFile = ['L'  num2str(nf) 'maxMinDosDim2_L' num2str(Len_stick) 'ks' num2str(10*ks) '.zraw'];      
    end      
        
    mha_write_volumeZ(fn3, vol2,1);
    disp('Write: ')     
    disp(fn3)   
       
    
     %%
      %=============================
      %Filtering in the Dim3 direction     
      
     
      %cast to double type
      I=single(vol1.pixelData);
      disp(size(I))             

     
     %%
      %max-min iteration
  
    for iter=1:2

      for dim=3 : 1 : 3
        %----------------------------------------------------------------------
        %The second level inner loop (Image slice)
        for Indz=1 : size(I,dim)
            disp(['dim=' num2str(dim) ' Slice=' num2str(Indz)])

            %Extract a 2D slice from the volume image
            switch dim
                case 1
                     I2D=I(Indz,:,:);   
                case 2
                     I2D=I(:,Indz,:);   
                case 3
                     I2D=I(:,:,Indz);  
            end                        
            I2D=squeeze(I2D); %Squeeze to 2D      
            Iorg=I2D; % Save as the original 2D image        
            [M, N]=size(I2D); %size of 2D slice        
            [X0,Y0]=meshgrid(1 :M,1:N);%grid coordinates

            %To save the results of different combinations        
            Imax_max=zeros(size(I2D));          

            %------------------------------------------------------------------
            %The third level inner loop (Stick orientation)

            %For sticks in the North_south directions
            for t1=1:(2*L+1)
                %The temporary variables 
                Imean=zeros(size(X0));  %for the mean along a stick        
                Istd=zeros(size(X0));   %for the standard deviation(STD)

                %Calculate the mean and STD
                for t2=1:(2*L+1)
                    %X offsets
                    X=X0+SKxy(t1,1,t2);
                    %Limit inside the image boundary
                    X(find(X<1))=1;
                    X(find(X>M))=M;

                    %Y offsets 
                    Y=Y0+SKxy(t1,2,t2);
                    %Limit inside the image boundary
                    Y(find(Y<1))=1;
                    Y(find(Y>N))=N;

                    %the current linear indices
                    inds=X+(Y-1)*M;    

                    %the sum of mean and STD
                    Imean=Imean+I2D(inds);      
                    Istd=Istd+I2D(inds).^2;   
                end
                %the mean value
                Imean=Imean/(2*L+1); 
                %the STD
                Istd=Istd/(2*L+1);      
                Istd=(Istd-Imean.^2).^0.5;           

                %calculate the derivative of stick (DoS) 

                %the linear indices of center points
                [Mc,Nc]=size(X0);
                [XC,YC]=meshgrid(1:Mc, 1:Nc);
                indsCenter=XC+(YC-1)*Mc;   

                %the linear indices of left points
                XL=XC-step;
                XL(find(XL<1))=1;
                indsLeft=XL+(YC-1)*Mc;             
                %the left derivatives
                IderivL=( Imean( indsCenter)-Imean( indsLeft));

                %the linear indices of right points
                XR=XC+step;
                XR(find(XR>Mc))=Mc;
                indsRight=XR+(YC-1)*Mc; 
                %the right derivative
                IderivR=(Imean( indsCenter)-Imean( indsRight));    

                %the nonlinear derivatives
                switch iter
                    case 1
                        %IderivES_min=min(IderivL,IderivR); 
                        IderivES_max=max(IderivL,IderivR);
                        %IderivES_max=(IderivL+IderivR)/2;   
                    case 2
                        %IderivES_min=(IderivL+IderivR)/2;   
                        %IderivES_max=(IderivL+IderivR)/2;   
                        IderivES_max=min(IderivL,IderivR);
                    case 3
                        %IderivES_min=(IderivL+IderivR)/2;   
                        %IderivES_max=(IderivL+IderivR)/2;        
                        IderivES_max=max(IderivL,IderivR);
                end               


                %Combine with the std               
                Ideriv_max=IderivES_max-ks*Istd( indsCenter);

                %Integration of responses in varying orientations                
                Imax_max=max(real(Imax_max), real(Ideriv_max));    
            end

            %For sticks in the East_west directions                   
            for t3=1:(2*L+1)
                %for the mean along a stick
                Imean=zeros(size(X0));  
                %for the standard deviation
                Istd=zeros(size(X0));  

                for t4=1:(2*L+1)
                    X=X0+SKyx(t3,1,t4);
                    X(find(X<1))=1;
                    X(find(X>M))=M;

                    Y=Y0+SKyx(t3,2,t4);
                    Y(find(Y<1))=1;
                    Y(find(Y>N))=N;
                    inds=X+(Y-1)*M;        
                    Imean=Imean+I2D(inds);      
                    Istd=Istd+I2D(inds).^2;   
                end

                Imean=Imean/(2*L+1);      

                Istd=Istd/(2*L+1);      
                Istd=(Istd-Imean.^2).^0.5;

                %the linear indices of center points
                [Mc,Nc]=size(X0);
                [XC,YC]=meshgrid(1:Mc, 1:Nc);
                indsCenter=XC+(YC-1)*Mc;                 

                %the derivative in the north direction
                YU=YC-step;
                YU(find(YU<1))=1;
                indsUp=XC+(YU-1)*Mc; 
                %up derivative
                IderivU=(Imean( indsCenter)-Imean( indsUp));

                %the derivative in the south direction
                YD=YC+step;
                YD(find(YD>Nc))=Nc;
                indsDown=XC+(YD-1)*Mc; 
                %up derivative
                IderivD=(Imean( indsCenter)-Imean(  indsDown));

                %the nonlinear derivative
                %IderivNS=min( IderivU, IderivD);             
                switch iter
                    case 1
                         %IderivNS_min=min( IderivU, IderivD); 
                         IderivNS_max=max( IderivU, IderivD); 
                         %IderivNS_max=(IderivU+IderivD)/2; 
                    case 2
                         %IderivNS_min=(IderivU+IderivD)/2; 
                         %IderivNS_max=(IderivU+IderivD)/2; 
                         IderivNS_max=min( IderivU, IderivD); 
                    case 3
                         %IderivNS_min=(IderivU+IderivD)/2; 
                         %IderivNS_max=(IderivU+IderivD)/2;    
                         IderivNS_max=max( IderivU, IderivD); 
                end                

                %Combine with the std            
                Ideriv_max=IderivNS_max-ks*Istd( indsCenter);

                %Integration of responses in varying orientations
                %Imin_max=max(real(Imin_max), real(Ideriv_min));    
                Imax_max=max(real(Imax_max), real(Ideriv_max));    
            end

            %------------------------------------------------------------------
            %Show the original and filtered images for debugging
            %figure(1)
            %subplot(1,2,1),  imshow(real(Iorg),[]), title('Original Image')  
            %Imax=Imax-min(min(Imax));
            %subplot(1,2,2),  imshow(real(Imax),[]), title('Filtered Image')  
            %pause

            %------------------------------------------------------------------
            %Save the 2D results as a 3D output
            switch dim
                case 3
                    %I3dMin(:,:,Indz)=single(Imin_max);    
                    %I3dMax(:,:,Indz)=single(Imax_max); 
                    I(:,:,Indz)=single(Imax_max); 
               case 2                
                    %I3dMin(:,Indz,:)=single(Imin_max);    
                    %I3dMax(:,Indz,:)=single(Imax_max); 
                    I(:,Indz,:)=single(Imax_max); 
                case 1                
                    %I3dMin(Indz,:,:)=single(Imin_max);    
                    %I3dMax(Indz,:,:)=single(Imax_max); 
                    I(Indz,:,:)=single(Imax_max); 
            end              

        end    
        I=single(I); %Cast to a single float type
      end


    end


    %%       
    vol2=vol1;
    vol2.pixelData=single(I);
    vol2.metaData.ElementType='MET_FLOAT';  
    vol2.metaData.CompressedData='True'; 

    if nf<10,
      fn3=[preP_dos '\Dim3\L0'  num2str(nf) 'maxMinDosDim3_L' num2str(Len_stick) 'ks' num2str(10*ks) '.mhd'];
      vol2.metaData.ElementDataFile = ['L0'  num2str(nf) 'maxMinDosDim3_L' num2str(Len_stick) 'ks' num2str(10*ks) '.zraw'];     
    else
      fn3=[preP_dos '\Dim3\L'  num2str(nf) 'maxMinDosDim3_L' num2str(Len_stick) 'ks' num2str(10*ks) '.mhd'];
      vol2.metaData.ElementDataFile = ['L'  num2str(nf) 'maxMinDosDim3_L' num2str(Len_stick) 'ks' num2str(10*ks) '.zraw'];      
    end      
    mha_write_volumeZ(fn3, vol2,1);
    disp('Write: ')     
    disp(fn3)
    
end
 
 
 

return
