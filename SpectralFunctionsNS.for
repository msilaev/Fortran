!!!!  This Fortran code calculates \Delta (x) and spectral functions in S'S system with spatially-dependent diffusion constant
!!!!  Spectral functions are calculated (x,E) 
 
      program a1111

      real*8 Delta01
      real*8 Temperature,OmegaM,Oom,pi   
        real*8 En, Enr,Eh,h,a,b,cosa,chb, StepE,signH
        real*8 g3Om(10000),g1Om(10000),hEx(10000),SOC(10000),
     *signb, rhoVar(10000),hEx0(10000), Magn(10000),
     *Magn1(10000),FermiL(10000),DuffVar(10000),hExPrev,Dev2

        complex*16 Enc,j0,g1(10000),g3(10000),tetaF(10000), 
     * teta3(10000)
 
      complex*16 g13(10000),g30(10000),g33(10000),g11(10000),
     * g10(10000)
      complex*16 g3P(10000),g3M(10000), g1P(10000),g1M(10000)   
          
      complex*16 cteta(10000),steta(10000),
     *M0(10000),d2M0(10000)
      
      real*8 ro11,xa,xb,xc,xd
      real*8 Diff1, MF, dnE, IntTheta, IntThetaR,  IntThetaIm,
     *IntThetaRr,IntThetaImr,Mvort,Rsuper,RNl

      real*8 Gorkov1R,Sparam
      real*8 GGorkov1R, GammaDi
     
      real*8 Op1,OpPrev1,OpNew1
      real*8 OpArr1(2000),AssTetaR1,AssTetaIm1, OpMax
      real*8  DiffN

      real*8 Op1rr(10000),OpPrev1r(10000),OpNew1r(10000),
     *Gorkov1Rr(10000),GGorkov1Rr(10000),G1r(10000),Dev1
      
      real*8 DiffVar(10000) 
              
       integer MatsMax,Nm,Nm1,Temp,Iter,dG,Jos,dr,dRMax,
     *IndMaxR,IndMaxR1, dT, dE,NEmax,dE1,NOmMax,dRsuper,
     *dDiffN

       character*300 fileG,fileDOS,fileOP,fileOP1

       character*1 nameG1,nameT1,nameR1,nameDiffN1
       character*2 nameG2,nameT2,nameR2,nameG,nameDiffN
       character*3 nameG3,nameR,nameR3,nameT,nameDN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8 Rr(10000) 
        real*8 tetaR1(10000),tetaIm1(10000),
     * tetaOm(10000),tetaImPrev2(10000), tetaRPrev2(10000),
     * tetaImPrev1(10000), tetaRPrev1(10000)

        real*8 D2tetaR(10000), D1tetaR(10000)

        real*8 Rmax,Rmax1,Rmax0
        real*8 r,r1,StepR  
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

        common/signH/signH
        common/g1/g1
        common/g3/g3
        common/g13/g13
        common/g33/g33
        common/g11/g11
  
        common/g10/g10
        common/g30/g30

        common/cteta/cteta
        common/steta/steta
        common/M0/M0
        common/d2M0/d2M0
        
        common/g3P/g3P
        common/g3M/g3M
        common/g1P/g1P
        common/g1M/g1M 

        common/g3Om/g3Om
        common/g1Om/g1Om 
        common/tetaF/tetaF
        common/dT/dT
           
        common/Diff1/Diff1 
        common/DiffVar/DiffVar 
        common/IndIter/IndIter
 !       common/AssTetaR1/AssTetaR1                     
 !       common/AssTetaR2/AssTetaR2                     

        common/En/En
        common/Enr/Enr
        common/Eh/Eh
        common/hEx/hEx
        common/MF/MF
        common/SOC/SOC

        common/Mvort/Mvort  
        common/pi/pi    
        common/Op1rr/Op1rr
        common/StepR/StepR
        common/Rmax/Rmax                  
      
        common/tetaR1/tetaR1   
        common/tetaOm/tetaOm   
        common/tetaIm1/tetaIm1      
        common/D1teta1R/D1teta1R            
         
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
        common/OmegaM/OmegaM

        j0=(0.d0,1.d0)
        pi=3.141592654d0 
       
        StepR=0.01d0   
        StepE=0.05d0 
        StepE=0.02d0 
 !       IndMaxR1=Rmax1/StepR+1                   
 !       IndMaxR=Rmax/StepR+1
                
        Mvort=0.d0
        Diff1=1.d0    
        pi=3.141592654d0                                
      
        ro11= 0.2d0        
           
        MF=0.d0
            
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
!        do 222 dRmax=20,20

          do 222 dRmax=21,21
       
!        Rmax =1.d0+dRmax*0.1d0       
         Rmax =1.0d0+(dRmax-1.d0)*0.1d0       
   
         do 888 dDiffN=10,10
         if(dDiffN<10) then
            write(nameDiffN1,'(i1)') dDiffN
          nameDiffN='0'//nameDiffN1 
         endif
          
          if(dDiffN>9) then
            write(nameDiffN,'(i2)') dDiffN
         endif

  
          DiffN=0.05d0 + (dDiffN-1.d0)*0.1d0
        
          do 777 dRsuper=1,16
          if(dRsuper==1) then 
          nameDN='0d0'  
          Rsuper=Rmax-0.00000d0
          endif
          
          if(dRsuper==2) then 
          nameDN='0d1'  
          Rsuper=Rmax-0.10000d0
          endif
          
          if(dRsuper==3) then 
          nameDN='0d2'  
          Rsuper=Rmax-0.20000d0
          endif
 
          if(dRsuper==4) then 
          nameDN='0d3'  
          Rsuper=Rmax-0.300000d0
          endif
          
          if(dRsuper==5) then 
          nameDN='0d4'  
          Rsuper=Rmax-0.400000d0
          endif

          if(dRsuper==6) then 
          nameDN='0d5'  
          Rsuper=Rmax-0.500000d0
          endif
          
          if(dRsuper==7) then 
          nameDN='0d6'  
          Rsuper=Rmax-0.600000d0
          endif

          if(dRsuper==8) then 
          nameDN='0d7'  
          Rsuper=Rmax-0.700000d0
          endif
          
          if(dRsuper==9) then 
          nameDN='0d8'  
          Rsuper=Rmax-0.800000d0
          endif
          
          if(dRsuper==10) then 
          nameDN='0d9'  
          Rsuper=Rmax-0.900000d0
          endif
  
          if(dRsuper==11) then 
          nameDN='1d0'  
          Rsuper=Rmax-1.00000d0
          endif
          
          if(dRsuper==12) then 
          nameDN='1d1'  
          Rsuper=Rmax-1.100000d0
          endif

          if(dRsuper==13) then 
          nameDN='1d2'  
          Rsuper=Rmax-1.200000d0
          endif
          
          if(dRsuper==14) then 
          nameDN='1d3'  
          Rsuper=Rmax-1.300000d0
          endif
          
          if(dRsuper==15) then 
          nameDN='1d4'  
          Rsuper=Rmax-1.400000d0
          endif 
          
          if(dRsuper==16) then 
          nameDN='1d5'  
          Rsuper=Rmax-1.500000d0
          endif

                 
         RNl=0.8d0
                   
        Rmax1=Rmax
        Rmax0=Rmax
          
        IndMaxR1=Rmax1/StepR+1                    
        IndMaxR=Rmax/StepR+1

!        Rsuper=7.d0*Rmax/8.d0
         
!        Rsuper=Rmax+10.3d0
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(dRMax<10) then     
      write(nameR1,'(i1)') dRMax 
      nameR='00'//nameR1 
      endif 

      if(dRMax<100.AND.dRMax>9) then 
      write(nameR2,'(i2)') dRMax
      nameR='0'//nameR2 
      endif  

!      if(dRMax<1000.AND.dRMax>99) then 
!      write(nameR,'(i3)') dRMax
!      !nameR='0'//nameR2 
!      endif  

!       fileDOS='chGDCCh6Rmax'//nameR//'Gamma01DS.dat'             
!       open(104,file=fileDOS)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !      do 111 dT=1,220

         do 111 dT=1,25
         
        write(*,*) dT, dRsuper,dDiffN   

      Temp=(dT-1)*4+1
!          Temp=(57-dT)

!      Temp=(35+dT)

      if(Temp<10) then    
      write(nameT1,'(i1)') Temp      
      nameT='00'//nameT1                         
      endif
      
      if(Temp<100.AND.Temp>9) then          
      write(nameT2,'(i2)') Temp  
      nameT='0'//nameT2                   
      endif   

      if(Temp<1000.AND.Temp>99) then          
      write(nameT,'(i3)') Temp  
!      nameT='0'//nameT2                   
      endif      

 !    fileDOS='DOST'//nameT//'R'//nameR//'D1CC.dat'             
      
      GammaDi=0.1d0           

!      fileDOS=
!     *'DOST'//nameT//'h3Rmax'//nameR//
!     *'D3d0DN0d0DiffNonConst.dat'  

         
      fileDOS='... Path to the folder with results..'//
     *'ResVarDGamma01\DOST'//nameT//'h3Rmax'//nameR//
     *'D3d0DN'//nameDN//'Diff'//nameDiffN//'CplVar.dat'             


 
      open(101,file=fileDOS)    

 !     fileOP='chOPT'//nameT//'h6Rmax'//nameR//'.dat' 

!      fileOP='OPT'//nameT//'h3Rmax'//nameR//
!     *'D3d0DN0d0DiffNonConst.dat'   

      fileOP='... Path to the folder with results..'//
     *'ResVarDGamma01\OPT'//nameT//'h3Rmax'//nameR//
     *'D3d0DN'//nameDN//'Diff'//nameDiffN//'CplVar.dat'


                 
      Temperature=Temp*0.01d0
      Sparam=1.d0/ro11-dlog(Temperature)        
!!!!!!!!! Set h distr
        
      do 499 dr=1,IndMaxR
 !     r=(dr-1.d0)*StepR-Rmax/8.d0
 !     r1=(IndMaxR-1.d0)*StepR-Rmax/8.d0

      r=(dr-1.d0)*StepR-0.1d0
      r1=(dr-1.d0)*StepR-Rsuper
      
      hEx0(dr)=0.0d0*( 1.d0 -dtanh(r/0.025d0))/2.d0 +
     * 0.0d0*( 1.d0 + dtanh(r1/0.025d0))/2.d0 +0.d0

         if(r>Rsuper) then
        hEx0(dr)=0.d0
        endif
 
      r=(dr-1.d0)*StepR-Rsuper
      r1=(dr-1.d0)*StepR-RNl

      rhoVar(dr) =
     *0.95d0*ro11*( 1.d0 - dtanh(r/0.025d0))/2.d0 + 
     *0.05d0*ro11 
	 
       DiffN = 0.05d0
          
       DiffVar(dr) =
     * (1.d0-DiffN)*( 1.d0 - dtanh(r/0.025d0))/2.d0 + DiffN

     
!     rhoVar(dr) =
!     *0.95d0*ro11*(( 1.d0 - dtanh(r/0.025d0))/2.d0)* 
!     *(( 1.d0 + dtanh(r1/0.025d0))/2.d0) + 
!     *0.05d0*ro11 
     
!     rhoVar(dr) =ro11

        if (dRsuper==1) then 
      rhoVar(dr) =ro11 
       DiffVar(dr) =1.d0
       endif

      FermiL(dr)=0.d0
            
499   continue   
                                  
                       
!!!!!!!! Set initial OP distribution !!!!!!!!!!!!!!!!
 !     if(dT==1) then
 
          if(dT>0) then
 
      do 501 dr=1,IndMaxR
      r=(dr-1.d0)*StepR

      hEx(dr)= hEx0(dr)
       
      Op1=1.78d0

      Op1rr(dr)=Op1

        if(r>Rsuper) then
        Op1rr(dr)=0.d0
        endif
     
501    continue
 !      close(99)    
       endif

 !!!!!!!!!!!!!!!!   Order Parameter  !!!!!!!!!!!!!!!!!!!!!!!!!!!
!ccccc Iterations for self consistency  ccccccc       

       Iter=0    
       Dev2=1.d0      
!       if(Dev2>0.0001d0) then         
194    continue

      Dev1=1.d0
!      if(Dev1>0.0001d0) then         
           
193   Iter=Iter+1   
         
      do 505 dr=1,IndMaxR
      OpNew1r(dr) =0.d0    
      Magn(dr)=0.d0
505   continue

      NOmMax=30/(2.d0*pi*Temperature)
      do 162 Nm1=1,NOmMax

!      write(*,*) 'Nm', Nm

       Nm=NOmMax +1 -Nm1
                 
      OmegaM=(2.d0*Nm-1.d0)*pi*Temperature  
      Enr=OmegaM
      En=0.d0
      signH=1.d0
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

       AssTetaR1 = 0.d0
       AssTetaIm1 = AssTetaR1
                       
!!!!!!!!!!!!!! Set initial theta distribution!!!!!!!!!!!!!!!!!!!
      if (Nm1==1) then
      do 114 dr=1,IndMaxR                    
      r=(dr-1.d0)*StepR              
      tetaIm1(dr)=0.d0           
      tetaR1(dr)=0.d0              
114   continue
      endif                                    

      call DOS2(AssTetaR1)  

       do 127 dr=1,IndMaxR 
       g10(dr)=-g1(dr)  
127    continue
               
  !!!!!!!!!! Set initial theta distr !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!! Contribution to OP

      do 504 dr=1,IndMaxR
      r=(dr-1.d0)*StepR

      Gorkov1Rr(dr)=dreal(g10(dr))
      GGorkov1Rr(dr)= Gorkov1Rr(dr)-Op1rr(dr)/OmegaM                                 
                   
      OpNew1r(dr)=OpNew1r(dr)+
     *2.d0*pi*Temperature*(rhoVar(dr)*GGorkov1Rr(dr))
      
       Magn(dr) = Magn(dr) + 
     *2.d0*pi*Temperature*dimag(g3(dr))
 
504    continue   
 
162   continue   
!!!!!!!!!!!!!!!! END Summation Omega !!!!!!!!!!!!
       Dev1=0.d0
      
!!!!!! update OP 
        open(102,file=fileOP)    

        do 506 dr=1,IndMaxR
     
       OpNew1r(dr)=OpNew1r(dr) + rhoVar(dr)*Op1rr(dr)*Sparam
 
!!!!!!!!!!!!                 
!      OpPrev1r(dr)=OpNew1r(dr)
       OpPrev1r(dr)=Op1rr(dr)                                     
       Op1rr(dr)=OpNew1r(dr)
!!!!!!!!!!!
 
       Magn1(dr) = Magn(dr) + hEx(dr)
  !    hEx(dr)=hEx0(dr) + FermiL(dr)*Magn(dr)/(1.d0-FermiL(dr))

 !!!!!!!!!! comparison !!!!!!!!!!!!!!!!
       
       if(dabs(Op1rr(dr)-OpPrev1r(dr))>Dev1) then                        
       Dev1 = dabs(Op1rr(dr)-OpPrev1r(dr))
       endif   
        
        write(102,200) (dr-1.d0)*StepR, Op1rr(dr),Magn(dr),
     * hEx(dr)
!       write(102,200) (dr-1.d0)*StepR, Op1rr(dr),hEx(dr)
     
506    continue
         close(102)

                            
         if(Dev1>0.0001d0) then 
        goto 193     
        endif ! OP iteratrion

        Dev2=0.d0
        do 187 dr=1,IndMaxR 
        hExPrev=hEx(dr)
  !      hEx(dr)=(hEx0(dr) + FermiL(dr)*Magn(dr))/(1.d0-FermiL(dr))
        hEx(dr)= hEx0(dr) + FermiL(dr)*( hEx(dr)+ Magn(dr) )
           
         if(dabs(hEx(dr)-hExPrev)>Dev2) then                        
       Dev2 = dabs(hEx(dr)-hExPrev)
       endif  
       
187    continue
!         write(*,*) 'Dev2', Dev2, Temp       
             
         if(Dev2>0.0001d0) then  
        goto 194  
        endif  ! Exch F iteration

!!!!!!!!!!!!!!!!!!!!!!!!!!!111

!        write(*,*) 'Op',Temperature,Op1rr(1),Magn(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!!!!!!!!!!!!!!!! DOS  !!!!!!!!!!!!!!!!!!!!!!!!!!!      
      OmegaM=GammaDi
                    
      Enr=OmegaM

  !    IntTheta=0.d0
      IntThetaR=0.d0
      IntThetaIm=0.d0

 
          NEmax=(3.5d0/StepE)

 !       do 31 dE=1,1
       do 31 dE1=1,2*NEmax

         dE=NEmax+1-dE1
    
      En=(dE-1.d0)*StepE + 0.01d0     
           
!!!!!!!!!!!!!! Set initial theta distribution!!!!!!!!!!!!!!!!!!!

 
       AssTetaIm1 =0.d0
       AssTetaR1=0.d0
                       
       do 115 dr=1,IndMaxR                   
               
       if (dE1>1) then
       tetaIm1(dr)=   tetaImPrev1(dr) 
       tetaR1(dr)=   tetaRPrev1(dr)
        endif
       
       if (dE1<2) then
       tetaIm1(dr)=AssTetaIm1
       tetaR1(dr)=AssTetaR1
       endif
                                
115    continue           
     
      tetaIm1(IndMaxR)= tetaIm1(IndMaxR-1)    
      tetaR1(IndMaxR)= tetaR1(IndMaxR-1) 
      
  !    En=(dE-1.d0)*StepE+h*0.d0    
      signH=1.d0    
 
       call DOS2(AssTetaR1)    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
       do 117 dr=1,IndMaxR                    
       
       tetaImPrev1(dr)=tetaIm1(dr) 
       tetaRPrev1(dr)=tetaR1(dr)
                                       
117    continue           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do 6031 dr=1,IndMaxR            
       g13(dr)= g1(dr)
       g33(dr)= g3(dr)
       
       g10(dr)= g1(dr)
       g30(dr)= g3(dr)     
        
       g3P(dr)=g3(dr)
       g1P(dr)=g1(dr)
       
6031   continue   
           
 !      En=(dE-1.d0)*StepE- h*0.d0       
       signH=-1.d0 
 !      Enc= j0*En-Enr
           
 !     call DOS1(AssTetaR1)        
 
       do 113 dr=1,IndMaxR    

        if (dE1>1) then
       tetaIm1(dr)=   tetaImPrev2(dr) 
       tetaR1(dr)=   tetaRPrev2(dr)
        endif                
                   
       if (dE1<2) then
       tetaIm1(dr)=AssTetaIm1
       tetaR1(dr)=AssTetaR1
       endif
                                
113    continue        
      
        call DOS2(AssTetaR1)      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
       do 118 dr=1,IndMaxR                    
       
       tetaImPrev2(dr)=tetaIm1(dr) 
       tetaRPrev2(dr)=tetaR1(dr)
                                       
118    continue           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do 6032 dr=1,IndMaxR            
        g13(dr)=  (g13(dr) -g1(dr))/2.d0 
        g10(dr)=  (g10(dr) +g1(dr))/2.d0
   
        g30(dr)=  (g30(dr) +g3(dr))/2.d0 
        g33(dr)=  (g33(dr) -g3(dr))/2.d0   
 
        g3M(dr)=g3(dr)
        g1M(dr)=g1(dr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
        xa=dreal((g30(dr))**2/( (g30(dr))**2+(g10(dr))**2))
        xb=dimag((g30(dr))**2/( (g30(dr))**2+(g10(dr))**2))
            
          xc=( ((xa**2+xb**2)**(0.5) + xa)/2.d0 )**(0.5)
          xd=( ((xa**2+xb**2)**(0.5) - xa)/2.d0 )**(0.5)

          cteta(dr)=xc+j0*xd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          xa=dreal( (g30(dr))**2+(g10(dr))**2)
          xb=dimag( (g30(dr))**2+(g10(dr))**2)

          signb=(xb)/(dabs(xb)+0.000000001d0)
            
          xc=( ((xa**2+xb**2)**(0.5) + xa )/2.d0 )**(0.5)
          xd=( ((xa**2+xb**2)**(0.5) - xa )/2.d0 )**(0.5)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
 
6032    continue 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!ccccccc Write to file-1  
                   
        do 502 dr=1,IndMaxR

       write(101,200) En, (dr-1.d0)*StepR,
     * dreal(g30(dr)), dimag(g30(dr)) , 
     * dreal(g10(dr)), dimag(g10(dr)) !, 


502   continue 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

31     continue     
      close(101)   
!      close(106)  
200   format(15f20.10)
    
        
111   continue  
      close(103)
      close(104)
      
777    continue  ! edn scan by dN

888    continue  ! end scan by DiffN

222   continue        
                    
      end     
    
!  cccccccccccccccccccc SUBROUTINE ccccccccccccccccccccccccc
      subroutine DOS2(AssTetaR1)        
      double precision hEx(10000)
      complex*16 g1(10000),j0,g3(10000)
      double precision  signH,MF    

      real*8 DeviationFRIm,AssTetaR1,AssTetaIm1,OmegaM     
      integer dr,IndMaxR,k,k1,IndIter
     
      real*8 pi,Temperature,En,Enr,Op1rr(10000),StepR,Rmax
      real*8 tetaR1(10000),tetaIm1(10000),DiffVar(10000),Mvort
      complex*16 teta(10000),FRIm(10000)                        

      common/signH/signH
    
      common/OmegaM/OmegaM 
      common/IndIter/IndIter
      common/En/En
      common/Enr/Enr
      common/hEx/hEx
      common/MF/MF
      
      common/Mvort/Mvort   
      
      common/FRIm/FRIm           
      common/Op1rr/Op1rr      
      common/g3/g3      
      common/g1/g1
           
        common/DiffVar/DiffVar   
      common/StepR/StepR    
      common/tetaR1/tetaR1
      common/tetaIm1/tetaIm1
      common/teta/teta
      common/Rmax/Rmax
    
      j0=(0.d0,1.d0) 
      IndMaxR=Rmax/StepR+1                                          
    
501   call progonkaRIm()
                    
! Calculate the deviation 
      DeviationFRIm=0.d0      
               
      do 21 dr=2,IndMaxR-1 
      if(DeviationFRIm < cdabs(FRIm(dr))) then
      DeviationFRIm = cdabs(FRIm(dr))  
      endif
       
21    continue             
         
200   format(10f20.10)                   
  
!      write(*,*) 'devF', DeviationFRIm

      if((DeviationFRIm>0.00000000001d0)) then 
      goto 501  
      write(*,*) 'devFRIm', DeviationFRIm
      endif  
                  
      do 20 dr=1,IndMaxR
      g3(dr)=dcos(tetaR1(dr))*dcosh(tetaIm1(dr))-
     *j0*dsin(tetaR1(dr))*dsinh(tetaIm1(dr))

      g1(dr)=dsin(tetaR1(dr))*dcosh(tetaIm1(dr))+ 
     *j0*dcos(tetaR1(dr))*dsinh(tetaIm1(dr))
     
!     write(*,*) 'teta', tetaR1(dr)
20    continue
!     pause         
       end

!cccccccccccccccccccc SUBROUTINE ccccccccccccccccccccccccc
      subroutine progonkaRIm()    
      double precision Op1rr(10000)      
      double precision D2tetaR(10000),D1tetaR(10000),DtetaR(10000)
      double precision D2tetaIm(10000),D1tetaIm(10000),DtetaIm(10000) 
      double precision tetaR1(10000),tetaIm1(10000),hEx(10000)
     
      double precision En,Enr,MF,En1,signH,Mvort      
      
      integer IndMaxR,k1,k,IndIter

      double precision StepR,Rmax,OmegaM,Opp,r,Diff1
      real*8 AIm(10000),BIm(10000),CIm(10000)

      real*8 FIm1(10000),GIm(10000),FR1(10000),GR(10000)
      real*8 DiffVar(10000)
  
      complex*16 A(10000),B(10000),C(10000),FRIm(10000),G1(10000)
      complex*16 D2teta(10000),Dteta(10000),alpha(10000),beta(10000)
      complex*16 Ec,teta(10000),j0

      common/signH/signH                                     
      common/Op1rr/Op1rr
      common/OmegaM/OmegaM 
      common/Diff1/Diff1  
      common/hEx/hEx
      common/MF/MF

      common/DiffVar/DiffVar    
      common/Mvort/Mvort      
      common/IndIter/IndIter
      common/FRIm/FRIm
      
      common/StepR/StepR
      common/Rmax/Rmax
      common/tetaR1/tetaR1
      common/tetaIm1/tetaIm1
      common/teta/teta
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
                    
      common/En/En
      common/Enr/Enr

        IndMaxR=Rmax/StepR+1
         j0=(0.d0,1.d0)
         
        alpha(1)=1.d0
         beta(1)=0.d0
!!  for zero derivaive b.c. 
            alpha(2)=1.d0
            beta(2)=0.d0        
    
!!!        for zero b.c. at FM see also line 753
!                 alpha(2)=0.d0
!           beta(2)=0.d0        

        do 20 k=2,IndMaxR-1     
        r=(k-1.d0)*StepR    
         teta(k+1)=tetaR1(k+1) + j0*tetaIm1(k+1)
         teta(k-1)=tetaR1(k-1) + j0*tetaIm1(k-1)
         teta(k)=tetaR1(k) + j0*tetaIm1(k)
    
        Opp=Op1rr(k)    
        En1=En+signH*hEx(k)
        Ec=En1 +j0*Enr
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        D2teta(k)=(teta(k+1)+teta(k-1)-2.d0*teta(k))/StepR**2
        
         D2teta(k)=DiffVar(k+1)*(teta(k+1)-teta(k))/StepR**2 + 
     * DiffVar(k)*(teta(k-1)-teta(k))/StepR**2
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        G1(k) = j0*Ec*(- j0*dsin(tetaR1(k))*dsinh(tetaIm1(k)) + 
     * dcos(tetaR1(k))*dcosh(tetaIm1(k)) ) +
     * Opp*( dcosh(tetaIm1(k))*dsin(tetaR1(k)) + 
     * j0*dcos(tetaR1(k))*dsinh(tetaIm1(k)) ) 
     
        FRIm(k)=  ( -j0*Ec*( j0*dsinh(tetaIm1(k))*dcos(tetaR1(k)) + 
     * dcosh(tetaIm1(k))*dsin(tetaR1(k)) ) +
     *  Opp*( dcosh(tetaIm1(k))*dcos(tetaR1(k)) - 
     * j0*dsinh(tetaIm1(k))*dsin(tetaR1(k)) )  -
     * D2teta(k) )!*0.5d0

        A(k)=DiffVar(k)*StepR**(-2)
        B(k)=DiffVar(k+1)*StepR**(-2)
        C(k)=G1(k)-(DiffVar(k)+DiffVar(k+1))*StepR**(-2)
     
!    A(k)=StepR**(-2)
!   B(k)=StepR**(-2)
!   C(k)=G1(k)-2.d0*StepR**(-2)
    
      alpha(k+1)=-B(k)/(A(k)*alpha(k)+C(k))
      beta(k+1)=(FRIm(k)-A(k)*beta(k))/(A(k)*alpha(k)+C(k))
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        
20     continue
            
!!!       Zero derivative at N 
        Dteta(IndMaxR)=beta(IndMaxR)/(1.d0-alpha(IndMaxR))
        
!!!         Zero function at N
!        Dteta(IndMaxR)=0.d0
    
        do 21 k=2,IndMaxR-1 
        k1=IndMaxR+1-k      
       Dteta(k1)=alpha(k1+1)*Dteta(k1+1)+beta(k1+1)        
21     continue
! For zero derovaive       
         Dteta(1)=Dteta(2)     

!!!!    For zero b.c. at FM
!         Dteta(1)=0.d0     


 !       open(10,file='pressFRIM.dat')           
 
      do 22 k=1,IndMaxR   
      r=(k-1.d0)*StepR                
      teta(k)=teta(k)+Dteta(k)
        tetaR1(k)=dreal(teta(k))
        tetaIm1(k)=dimag(teta(k)) 
         
!        write(10,200) r, dreal(FRIm(k)), dimag(FRIm(k))                

22      continue      
200    format(10f20.10)    

!   close(10)    

 !       pause     
        end
    


    
