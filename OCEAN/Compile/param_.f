      module param
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=41,   MMm0=42,   N=32)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 N_sl
      parameter (N_sl=0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=1,  NP_ETA=4,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 stdout, Np, NpHz, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (NpHz=(N+1+N_sl))
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d,N3dHz, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=35)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      parameter (N3dHz=size_XI*size_ETA*NpHz)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, NTA, itemp, NTot
      integer*4   ntrc_temp, ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      integer*4   ntrc_subs, ntrc_substot
      parameter (itemp=1)
      parameter (ntrc_temp=1)
      parameter (ntrc_salt=1)
      parameter (ntrc_pas=0)
         parameter (ntrc_bio=9)
      parameter (ntrc_subs=0, ntrc_substot=0)
      parameter (ntrc_sed=0)
      parameter (NTA=itemp+ntrc_salt)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      parameter (NTot=NT)
      integer*4 NGLS
      parameter(NGLS=2)
      integer*4 itke
      parameter(itke=1)
      integer*4 igls
      parameter(igls=2)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_diapv
      integer*4   ntrc_diaeddy, ntrc_surf
     &          , itrc_bio
     &          , isalt
     &          , iDIC_, iTAL_, iOXY_, iCAL_, iPO4_
     &          , iPOC_, iSIL_, iPHY_, iZOO_, iDOC_
     &          , iDIA_, iMES_, iDSI_, iFER_
     &          , iBFE_, iGOC_, iSFE_, iDFE_, iGSI_
     &          , iNFE_, iNCH_, iDCH_, iNO3_, iNH4_
     &          , iLGW_, iDON_, iDOP_, iPON_, iPOP_
     &          , iNPH_, iPPH_, iNDI_, iPDI_, iPIC_
     &          , iNPI_, iPPI_, iPFE_, iPCH_, iGON_
     &          , iGOP_
     &          , Nhi,Nco3,Naksp,Netot,Nprorca
     &          , Nprorcad,Npronew,Npronewd
     &          , Nprobsi,Nprofed,Nprofen
     &          , Ngrapoc,Ngrapoc2
     &          , Nmico2,Nmeso2
     &          , Nnitrifo2,Nfixo2,Nremino2
     &          , Npronewo2,Nprorego2
     &          , Nfld,Nflu16,Nkgco2,Natcco2,Nsinking
     &          , Nsinkfer,Nsinksil,Nironsed
     &          , Nsinkcal,Nheup,Nnitrpot
     &          , Nirondep,Nsildep,Npo4dep
     &          , Nno3dep,Nnh4dep
     &          , NumFluxTerms,NumVSinkTerms,NumGasExcTerms
      parameter (isalt=itemp+1)
      parameter (itrc_bio=itemp+ntrc_salt+ntrc_pas+1)
      parameter (iDIC_=itrc_bio, iTAL_=iDIC_+1, iOXY_=iDIC_+2)
      parameter ( iPOC_=iDIC_+3,  iPHY_=iDIC_+4, iZOO_=iDIC_+5,
     &            iDOC_=iDIC_+6,  iNO3_=iDIC_+7, iFER_=iDIC_+8)
      parameter (iDON_=iDIC_+25, iDOP_=iDIC_+26, iPON_=iDIC_+27,
     &	         iPOP_=iDIC_+28, iNPH_=iDIC_+29, iPPH_=iDIC_+30,
     &	         iNDI_=iDIC_+31, iPDI_=iDIC_+32, iPIC_=iDIC_+33,
     &	         iNPI_=iDIC_+34, iPPI_=iDIC_+35, iPFE_=iDIC_+36,
     &	         iPCH_=iDIC_+37, iGON_=iDIC_+38, iGOP_=iDIC_+39)
      parameter (Nhi        = 1,
     &            Nco3      = 2,
     &            Naksp     = 3,
     &            Netot     = 4,
     &            Nprorca   = 5,
     &            Ngrapoc   = 6,
     &            Nmico2    = 7,
     &            Nremino2  = 8,
     &            Nfixo2    = 9,
     &            Nirondep  = 10,
     &            Nironsed  = 11,
     &            Npronewo2 = 12,
     &            Npronew   = 13,
     &            NumFluxTerms = Npronew)
       parameter (Nfld      = 1,
     &            Nflu16    = 2,
     &            Nkgco2    = 3,
     &            Natcco2   = 4,
     &            Nsinking  = 5,
     &            Nheup     = 6,
     &            Nno3dep   = 7,
     &            Nnitrpot  = 8,
     &            NumGasExcTerms = 0,
     &            NumVSinkTerms = Nnitrpot)
      parameter (ntrc_diabio=NumFluxTerms+
     &                       NumGasExcTerms+NumVSinkTerms)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=0)
      parameter (ntrc_surf=0)
      end module param
