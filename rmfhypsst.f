cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       calculate \lambda_c^+ hypernuclei
c       add \lambda_c^+ coulmb
c       2019/4/13      
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      calculate hypernuclei by RMF
c      modified from RMFQMF
c      add strange mesons "sigma*" and "phi"
c
c      shen 
c      2004/11/08
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      PROGRAM RMF1
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'BBCOUP.INC'
      INCLUDE 'CONDIT.INC'
      INCLUDE 'CONSTA.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'SPECIE.INC'
      INCLUDE 'YMSTAR.INC'
      INCLUDE 'YLAMDA.INC'
      INCLUDE 'LAMDA.INC'
      INCLUDE 'YSSMES.INC'
      DIMENSION F(5),DF(5)
C      
      DATA    PI,         ONE/        
     * 3.14159265358979,  1.0/
      DATA    HBARC,    HMC,      OMEG/
     *     197.32858, .02411,      41./
      DATA ITPAIR/
     *       -100/
      DATA ENINV/
     *      40.0/
C     
C
      OPEN(UNIT=10,FILE='output250.rmf',STATUS='UNKNOWN')
      OPEN(UNIT=11,FILE='output250.dat',STATUS='UNKNOWN')
      OPEN(UNIT=41,FILE='test250.rmf',STATUS='UNKNOWN')
      OPEN(UNIT=42,FILE='test2250.rmf',STATUS='UNKNOWN')
C
      OPEN(UNIT=05,FILE='input.rmf',STATUS='OLD')
      READ(05,*)
      READ(05,*) ENDCON,RSTEP,ADDNEW,ITERMX
      READ(05,*) 
      READ(05,*) IZPE,ISOVEC,IWRIT
      READ(05,*)
      READ(05,*) IPAIR,IFPAIR,DELPIN,DELNIN
      READ(05,*)
      READ(05,*) BARYMS,SMESMS,VMESMS,DMESMS,RMESMS
      READ(05,*)
      READ(05,*) BBCOUP(1,1),BBCOUP(1,2)
      READ(05,*) BBCOUP(2,1),BBCOUP(2,2)
      READ(05,*) BBCOUP(3,1),BBCOUP(3,2)
      READ(05,*) BBCOUP(4,1),BBCOUP(4,2)
      READ(05,*) BBCOUP(5,1),BBCOUP(5,2)
      READ(05,*) BBCOUP(6,1),BBCOUP(6,2)
      READ(05,*) POWER
      CLOSE(05)
C
      OPEN(UNIT=06,FILE='nucldt.rmf',STATUS='OLD')
      READ(06,*)
      READ(06,*) NNUCL
      IF(NNUCL.GT.NUCLMX) THEN
       WRITE(10,*) ' TOO MANY  NUCLEI !!'
       STOP
      ENDIF
      READ(06,*)
      DO 310 I=1,NNUCL
      READ(06,*) INPROT(I),INNEUT(I),RMAXIN(I),
     *           INSHEL(I,1),INSHEL(I,2)
     *          ,INLAMDA(I),IJSL(I),RGW(I)
  310 CONTINUE
C
      IF(IPAIR.EQ.0  .OR.  IPAIR.EQ.3) THEN
      READ(06,*)
       DO 320 I=1,NNUCL
        READ(06,*) 
        READ(06,*) (WEIGIN(I,J),J=1,INSHEL(I,1))
        READ(06,*) (WEIGIN(I,J),
     *              J=INSHEL(I,1)+1,INSHEL(I,1)+INSHEL(I,2))
  320  CONTINUE
      ENDIF
      CLOSE(06)
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C     ADD LAMDA

      OPEN(UNIT=04,FILE='qmf DIRAC250.DAT',STATUS='OLD')
      OPEN(UNIT=03,FILE='DIRAC250.CHK',STATUS='UNKNOWN')
C
      READ(04,*) NNY
      DO 100 I=1,NNY
      READ(04,*) YQM(I),YMS(I),YMSL(I)
      YQM(I)=YQM(I)/HBARC
      YMS(I)=YMS(I)/HBARC
      YMSL(I)=YMSL(I)/HBARC
  100 CONTINUE
C-------------------------
      YH=YQM(2)-YQM(1)
C
      I=3
      F(1)=YMS(I-2)
      F(2)=YMS(I-1)
      F(3)=YMS(I)
      F(4)=YMS(I+1)
      F(5)=YMS(I+2)
      CALL DERIV(F,DF,YH)
      YDMQ(1)=DF(1)
      YDMQ(2)=DF(2)
      YDMQ(3)=DF(3)
C
      DO 200 I=4,NNY-3
      F(1)=YMS(I-2)
      F(2)=YMS(I-1)
      F(3)=YMS(I)
      F(4)=YMS(I+1)
      F(5)=YMS(I+2)
      CALL DERIV(F,DF,YH)
      YDMQ(I)=DF(3)
  200 CONTINUE
C
      I=NNY-2
      F(1)=YMS(I-2)
      F(2)=YMS(I-1)
      F(3)=YMS(I)
      F(4)=YMS(I+1)
      F(5)=YMS(I+2)
      CALL DERIV(F,DF,YH)
      YDMQ(NNY-2)=DF(3)
      YDMQ(NNY-1)=DF(4)
      YDMQ(NNY)  =DF(5)
C-------------------------
C     NOTE YQM, YMS, YDMS ARE POSITIVE!
C          YSIG, YDMQ     ARE NEGATIVE!
C     YSIG,YMS  WITH UNIT FM^-1
C 
      GSQ=BBCOUP(1,1)
      DO 300 I=1,NNY
      YSIG(I)=-YQM(I)/GSQ
      YDMS(I)=-YDMQ(I)*GSQ
      WRITE(3,*) YSIG(I)*HBARC,YMS(I)*HBARC,YDMS(I)
  300 CONTINUE
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C     CALCULATE LAMDA
      YH=YQM(2)-YQM(1)
C
      I=3
      F(1)=YMSL(I-2)
      F(2)=YMSL(I-1)
      F(3)=YMSL(I)
      F(4)=YMSL(I+1)
      F(5)=YMSL(I+2)
      CALL DERIV(F,DF,YH)
      YDMQL(1)=DF(1)
      YDMQL(2)=DF(2)
      YDMQL(3)=DF(3)
C
      DO 400 I=4,NNY-3
      F(1)=YMSL(I-2)
      F(2)=YMSL(I-1)
      F(3)=YMSL(I)
      F(4)=YMSL(I+1)
      F(5)=YMSL(I+2)
      CALL DERIV(F,DF,YH)
      YDMQL(I)=DF(3)
  400 CONTINUE
C
      I=NNY-2
      F(1)=YMSL(I-2)
      F(2)=YMSL(I-1)
      F(3)=YMSL(I)
      F(4)=YMSL(I+1)
      F(5)=YMSL(I+2)
      CALL DERIV(F,DF,YH)
      YDMQL(NNY-2)=DF(3)
      YDMQL(NNY-1)=DF(4)
      YDMQL(NNY)  =DF(5)
C-------------------------
      DO 500 I=1,NNY
      YDMSL(I)=-YDMQL(I)*GSQ
      WRITE(3,*) YSIG(I)*HBARC,YMSL(I)*HBARC,YDMSL(I),YDMSL(I)/YDMS(I)
  500 CONTINUE
C

CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C     ADD strange mesons
C
      SSMESMS= 980.0D0/HBARC
	VSMESMS=1020.0D0/HBARC
      SSMASS2=SSMESMS**2
      VSMASS2=VSMESMS**2
C
      SQ2=DSQRT(2.0D0)
c      SSMESCPL=SQ2/3.0D0*bbcoup(1,1)
c      VSMESCPL=SQ2/3.0D0*bbcoup(2,1)
c      SSMESCPL=0.50*3.0*bbcoup(1,1)
      SSMESCPL=0.0D0*bbcoup(1,1)
      VSMESCPL=0.0D0*bbcoup(2,1)
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C     Start  calculation
C 
      DO 10 JS=1,JSMXMX
      XJ(JS)=FLOAT((JS-1)/2)+0.5D0
      XL(JS)=FLOAT(JS/2)
      XKAPPA(JS)=FLOAT((1-MOD(JS,2)*2)*((JS+1)/2))
      GPARIT(JS)=FLOAT(MOD(JS/2,2)*2-1)
      PHASFC(JS)=(XJ(JS)+XJ(JS)+1.D0)/(4.D0*PI)
   10 CONTINUE
C
C
C
      IF(IWRIT.GE.1) THEN
       WRITE(10,'(2(A,F10.4),/2X,4(A,F10.4))')
     *       '     PMASS=',BARYMS(1),'   NMASS=',BARYMS(2),
     *       '  SMESMS=',SMESMS,
     *       '  VMESMS=',VMESMS,'  DMESMS=',DMESMS,'  RMESMS=',RMESMS
       WRITE(10,'(A,F10.3,6(/,A,2F16.4))')
     *       '  POWER=',POWER,
     *       '  SMESCP=',BBCOUP(1,1),BBCOUP(1,2),
     *       '  VMESCP=',BBCOUP(2,1),BBCOUP(2,2),
     *       '  RMESCP=',BBCOUP(3,1),BBCOUP(3,2),
     *       '  SMSB2C=',BBCOUP(4,1),BBCOUP(4,2),
     *       '  SMSB3C=',BBCOUP(5,1),BBCOUP(5,2),
     *       '  VMCPC3=',BBCOUP(6,1),BBCOUP(6,2)
       WRITE(10,'(2(A,E11.3))') '  ENDCON=',ENDCON,
     *                           '  RSTEP=',RSTEP
       WRITE(10,'(A,I5)') '  IZPE=',IZPE
      ENDIF
      IF(ISOVEC.EQ.-1) THEN
       WRITE(10,'(/,A)') '  ISOVECOTR MESONS ARE TREATED UNNATURALLY !'
      ENDIF
C
      KERROR=0
      CALL RMFMN(KERROR)
C      
 1000 CONTINUE
      CLOSE(05)
      STOP
      END
C
C---   COULMB  ---------------------------------------------------------
C
      SUBROUTINE COULMB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C      
      INCLUDE 'CONSTA.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'WORKSP.INC'
      INCLUDE 'YLAMDA.INC'
C      
      DATA VPEQ3/6.0317561/
      DATA EZW,FSX/.0833333333,.833333333/
C
C
      VOR =3.0D0*VPEQ3*RSTEP/HBARC
      DO 20 I=1,NGRID
   20 AUX1(I)=VOR*(CRHO(I)+CRHOL(I))*(I-1)
      INTM=NGRID-1
      AUX2(1)=0.
      DO 21 I=2,INTM
   21 AUX2(I)=EZW*(AUX1(I+1)+AUX1(I-1))+FSX*AUX1(I)
     *       +AUX2(I-1)*(FLOAT(I-2))/(FLOAT(I-1))
      AUX1(NGRID)=(EZW*AUX1(INTM)+FSX*AUX1(NGRID)+AUX2(INTM)
     *       *(FLOAT(NGRID-2))/(FLOAT(INTM)))*(FLOAT(INTM))
     *          /(FLOAT(NGRID))*RSTEP*RSTEP
     *         +1.440000/HBARC*APROT*NGRID/(FLOAT(NGRIDP))
      FAC=RSTEP*RSTEP
      I=INTM
   22 AUX1(I)=(AUX2(I)*FAC+AUX1(I+1))*(I-1)/(FLOAT(I))
      I=I-1
      IF(I-1) 23,23,22
   23 CONTINUE
      AUX1(1)=0.
      COULPT(1)=(VDRH*AUX1(2)-ESXH*AUX1(3))
      DO 24 I=2,NGRID
   24 COULPT(I)=AUX1(I)/(I*RSTEP-RSTEP)
C
      RETURN
      END
C      
C---   DENSIT   --------------------------------------------------------
C
      SUBROUTINE DENSIT(ITER)      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'CONDIT.INC'
      INCLUDE 'CONSTA.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'WORKSP.INC'
      INCLUDE 'YLAMDA.INC'
      INCLUDE 'LAMDA.INC'
C
C
      IF(ITER.EQ.1) THEN
       ADDWGT=1.0
       FACRHO=0.0
      ELSE
       ADDWGT=ADDNEW
       FACRHO=(1.0-ADDNEW)/ADDNEW
      ENDIF
C
      DO 10 I=1,NGRID
      SRHO(I)=FACRHO*SRHO(I)
      VRHO(I)=FACRHO*VRHO(I)
      DRHO(I)=FACRHO*DRHO(I)
      RRHO(I)=FACRHO*RRHO(I)   
      CRHO(I)=VRHO(I)
   10 CONTINUE
C
C     NUCLEON SIDE
C
      DO 20 ISO=1,2
       FACISO=-3.0+2.0*ISO
       DO 25 JS=1,JSMAX(ISO)
       DO 25 NSHL=1,NSMAX(JS,ISO)
       NPLM=NPLACE(NSHL,JS,ISO)-1
       PHAS=PHASFC(JS)*WEIGHT(NPLM/NGRID2+1)
C
      IF(XKAPPA(JS).EQ.-1.) THEN
       DWFN=VDRH*WFN(NPLM+2)-ESXH*WFN(NPLM+3)
       SRHO(1)=PHAS*DWFN*DWFN+SRHO(1)
       VRHO(1)=PHAS*DWFN*DWFN+VRHO(1)
       DRHO(1)=PHAS*DWFN*DWFN*FACISO+DRHO(1)
       RRHO(1)=PHAS*DWFN*DWFN*FACISO+RRHO(1)
      ENDIF
      IF(XKAPPA(JS).EQ.1.) THEN
       DWFN=VDRH*WFN(NPLM+NGRID+2)-ESXH*WFN(NPLM+NGRID+3)
       SRHO(1)=-PHAS*DWFN*DWFN+SRHO(1)
       VRHO(1)=PHAS*DWFN*DWFN+VRHO(1)
       DRHO(1)=-PHAS*DWFN*DWFN*FACISO+DRHO(1)
       RRHO(1)=PHAS*DWFN*DWFN*FACISO+RRHO(1)
      ENDIF
C
      DO 30 I=2,NGRID
       PHASR2=PHAS/((I*RSTEP-RSTEP)*(I*RSTEP-RSTEP))
       SRHO(I)=PHASR2*(WFN(I+NPLM)*WFN(I+NPLM)
     *                -WFN(I+NPLM+NGRID)*WFN(I+NPLM+NGRID))
     *                 +SRHO(I)
       VRHO(I)=PHASR2*(WFN(I+NPLM)*WFN(I+NPLM)
     *                +WFN(I+NPLM+NGRID)*WFN(I+NPLM+NGRID))
     *                +VRHO(I)
       DRHO(I)=PHASR2*FACISO*(WFN(I+NPLM)*WFN(I+NPLM)
     *                -WFN(I+NPLM+NGRID)*WFN(I+NPLM+NGRID))
     *                +DRHO(I)
       RRHO(I)=PHASR2*FACISO*(WFN(I+NPLM)*WFN(I+NPLM)
     *                +WFN(I+NPLM+NGRID)*WFN(I+NPLM+NGRID))
     *                +RRHO(I)
   30 CONTINUE
C
   25 CONTINUE
      IF(ISO.EQ.1) THEN
       DO 40 I=1,NGRID
       CRHO(I)=VRHO(I)-CRHO(I)
   40  CONTINUE
      ENDIF
   20 CONTINUE
C
      IF(ITER.NE.1) THEN
       DO 60 I=1,NGRID
       SRHO(I)=ADDWGT*SRHO(I)
       VRHO(I)=ADDWGT*VRHO(I)
       DRHO(I)=ADDWGT*DRHO(I)
       RRHO(I)=ADDWGT*RRHO(I)
   60  CONTINUE
      ENDIF
C
      IF(ABS(VMESTC).GT.1.D-10) THEN
        DO 100 I=1,NGRID2
        AUX1(I)=0.0
  100   CONTINUE
        DO 110 ISO=1,2
          DO 120 JS=1,JSMAX(ISO)
          DO 120 NSHL=1,NSMAX(JS,ISO)
          NPLM=NPLACE(NSHL,JS,ISO)-1
          PHAS=PHASFC(JS)*WEIGHT(NPLM/NGRID2+1)
          DO 130 I=2,NGRID
          RA=(I*RSTEP-RSTEP)
          AUX1(I)=PHAS/(RA*RA*RA)*WFN(I+NPLM)*WFN(I+NPLM+NGRID)
     *           +AUX1(I)
  130     CONTINUE
C
  120     CONTINUE
  110   CONTINUE
        AUX1(1)=(4.0*AUX1(2)-AUX1(3))*0.33333330
        CALL DERIV0(AUX1,AUX1(NGRIDP),ONE)
        DO 140 I=1,NGRID
        RA=I*RSTEP-RSTEP
        IP=I+NGRID
        TERHO(I)=(3.00*AUX1(I)+RA*AUX1(IP)+FACRHO*TERHO(I))*ADDWGT
  140   CONTINUE
      ENDIF
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C     ADD LAMDA
C	
      IF (NLAMDA.EQ.0) GOTO 500
C
       JS=JSL
       NPLM=0
       PHAS=1.D0/(4.D0*PI)
C
      IF(XKAPPA(JS).EQ.-1.) THEN
       DWFN=VDRH*WFNL(NPLM+2)-ESXH*WFNL(NPLM+3)
       SRHOL(1)=PHAS*DWFN*DWFN
       VRHOL(1)=PHAS*DWFN*DWFN
      ENDIF
      IF(XKAPPA(JS).EQ.1.) THEN
       DWFN=VDRH*WFNL(NPLM+NGRID+2)-ESXH*WFNL(NPLM+NGRID+3)
       SRHOL(1)=-PHAS*DWFN*DWFN
       VRHOL(1)=PHAS*DWFN*DWFN
      ENDIF
C
      DO 230 I=2,NGRID
       PHASR2=PHAS/((I*RSTEP-RSTEP)*(I*RSTEP-RSTEP))
       SRHOL(I)=PHASR2*(WFNL(I+NPLM)*WFNL(I+NPLM)
     *                -WFNL(I+NPLM+NGRID)*WFNL(I+NPLM+NGRID))
       VRHOL(I)=PHASR2*(WFNL(I+NPLM)*WFNL(I+NPLM)
     *                +WFNL(I+NPLM+NGRID)*WFNL(I+NPLM+NGRID))
  230 CONTINUE
C
      IF(ABS(VMESTCL).GT.1.D-10) THEN
        DO 330 I=2,NGRID
        RA=(I*RSTEP-RSTEP)
        AUX1(I)=PHAS/(RA*RA*RA)*WFNL(I+NPLM)*WFNL(I+NPLM+NGRID)
  330   CONTINUE
C
        AUX1(1)=(4.0*AUX1(2)-AUX1(3))*0.33333330
        CALL DERIV0(AUX1,AUX1(NGRIDP),ONE)
        DO 340 I=1,NGRID
        RA=I*RSTEP-RSTEP
        IP=I+NGRID
        TERHOL(I)=3.00*AUX1(I)+RA*AUX1(IP)
C        TERHOL(I)=(3.00*AUX1(I)+RA*AUX1(IP)+FACRHO*TERHOL(I))*ADDWGT
  340   CONTINUE
      ENDIF
       
C
      DO 360 I=1,NGRID
      SRHOL(I) = NLAMDA*SRHOL(I)
      VRHOL(I) = NLAMDA*VRHOL(I)
      CRHOL(I)=VRHOL(I)
      TERHOL(I)= NLAMDA*TERHOL(I)
  360 CONTINUE
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
 500   CONTINUE
       RETURN
       END
C       
C---   DERIV0   --------------------------------------------------------
C
      SUBROUTINE DERIV0(WFIN,DWFOUT,PARITY)      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C      
      INCLUDE 'GRID.INC'
      DIMENSION WFIN(NGRDMX),DWFOUT(NGRDMX)
C
C
      IA=1
      IF(PARITY.GT.0.) THEN
      DWFOUT(IA)=0.
      ELSE
      DWFOUT(IA)=VDRH*WFIN(IA+1)-ESXH*WFIN(IA+2)
      ENDIF
      IA=IA+1
      DWFOUT(IA)=ZDRH*(WFIN(IA+1)-WFIN(IA-1))
     *           -EZWH*(WFIN(IA+2)-PARITY*WFIN(IA))
      IA=1+IA
      DO 10 I=IA,NGRID-2
      DWFOUT(I)=EZWH*(WFIN(I-2)-WFIN(I+2))
     *          +ZDRH*(WFIN(I+1)-WFIN(I-1))
   10 CONTINUE
      IA=NGRID-1
      DWFOUT(IA)=EVRH*WFIN(NGRID)+FSXH*WFIN(NGRID-1)-DR2H*WFIN(NGRID-2)
     *           +EHH*WFIN(NGRID-3)-EZWH*WFIN(NGRID-4)
      DWFOUT(NGRID)=FZWTH*WFIN(NGRID)-VRH*WFIN(NGRID-1)
     *              +DRH*WFIN(NGRID-2)-VDRH*WFIN(NGRID-3)
     *              +EVRH*WFIN(NGRID-4)
C
      RETURN
      END
C      
C---   DERIV1   --------------------------------------------------------
C
      SUBROUTINE DERIV1(WFIN,DWFOUT,PARITY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C      
      INCLUDE 'GRID.INC'
      DIMENSION WFIN(NGRDM2),DWFOUT(NGRDM2)
C
C
      IA=1
      IF(PARITY.GT.0.) THEN
      DWFOUT(IA)=0.
      ELSE
      DWFOUT(IA)=VDRH*WFIN(IA+1)-ESXH*WFIN(IA+2)
      ENDIF
      IA=IA+1
      DWFOUT(IA)=ZDRH*(WFIN(IA+1)-WFIN(IA-1))
     *           -EZWH*(WFIN(IA+2)-PARITY*WFIN(IA))
      IA=1+IA
      DO 10 I=IA,NGRID-2
      DWFOUT(I)=EZWH*(WFIN(I-2)-WFIN(I+2))
     *          +ZDRH*(WFIN(I+1)-WFIN(I-1))
   10 CONTINUE
      IA=NGRID-1
      DWFOUT(IA)=EVRH*WFIN(NGRID)+FSXH*WFIN(NGRID-1)-DR2H*WFIN(NGRID-2)
     *           +EHH*WFIN(NGRID-3)-EZWH*WFIN(NGRID-4)
      DWFOUT(NGRID)=FZWTH*WFIN(NGRID)-VRH*WFIN(NGRID-1)
     *              +DRH*WFIN(NGRID-2)-VDRH*WFIN(NGRID-3)
     *              +EVRH*WFIN(NGRID-4)
C
      IA=NGRIDP
      IF(PARITY.LT.0.) THEN
      DWFOUT(IA)=0.
      ELSE
      DWFOUT(IA)=VDRH*WFIN(IA+1)-ESXH*WFIN(IA+2)
      ENDIF
      IA=IA+1
      DWFOUT(IA)=ZDRH*(WFIN(IA+1)-WFIN(IA-1))
     *           -EZWH*(WFIN(IA+2)-PARITY*WFIN(IA))
      IA=1+IA
      DO 20 I=IA,NGRID2-2
      DWFOUT(I)=EZWH*(WFIN(I-2)-WFIN(I+2))
     *          +ZDRH*(WFIN(I+1)-WFIN(I-1))
   20 CONTINUE
      IA=NGRID2-1
      DWFOUT(IA)=EVRH*WFIN(NGRID2)+FSXH*WFIN(NGRID2-1)
     *           -DR2H*WFIN(NGRID2-2)
     *           +EHH*WFIN(NGRID2-3)-EZWH*WFIN(NGRID2-4)
      DWFOUT(NGRID2)=FZWTH*WFIN(NGRID2)-VRH*WFIN(NGRID2-1)
     *              +DRH*WFIN(NGRID2-2)-VDRH*WFIN(NGRID2-3)
     *              +EVRH*WFIN(NGRID2-4)
C
      RETURN
      END
C
C---   ENERGY   --------------------------------------------------------
C
      SUBROUTINE ENERGY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'CONDIT.INC'
      INCLUDE 'CONSTA.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'RESULT.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'YMSTAR.INC'
      INCLUDE 'YLAMDA.INC'
      INCLUDE 'LAMDA.INC'
      INCLUDE 'YSSMES.INC'
C
C
      SUMSPN=0.0
      DO 10 N=1,NSHLPR+NSHLNE
      SUMSPN=2.0*SPENRG(N)*WEIGHT(N)*XMULT(N)+SUMSPN
   10 CONTINUE
C
      SUMPON=0.0
      DO 20 I=2,NGRID
      RX=RSTEP*(I-1)
      SUMPON=-WGTINT(I)*RX**2
     *      *(SPOTD(I)*SRHO(I)+VPOT(I)*VRHO(I)+DPOT(I)*DRHO(I)
     *       +RPOT(I)*RRHO(I)+COULPT(I)*CRHO(I))+SUMPON
   20 CONTINUE
      SUMPON=12.5660*SUMPON*0.5D0*HBARC
C
      ENREAR=0.0
      DO 30 I=2,NGRID
       RX=RSTEP*(I-1)
       SFIELD=SPOTD(I)/SMESCPD(I)
       VFIELD=VPOT(I)/VMESCP
       IF(DMESCP.EQ.0.0) THEN
          DFIELD=0.D0
       ELSE       
          DFIELD=DPOT(I)/DMESCP
       ENDIF
       IF(RMESCP.EQ.0.0) THEN
          RFIELD=0.0
        ELSE
          RFIELD=RPOT(I)/RMESCP
        ENDIF
       ENREAR=WGTINT(I)*RX**2
     *       *(-SMSB2C*SFIELD**3/6.D0-SMSB3C*SFIELD**4/4.D0
     *         +VMCPC3*VFIELD**4/4.D0
     *         -ISOVEC*DMCPD3*DFIELD**4/4.D0
     *         +ISOVEC*RMCPE3*RFIELD**4/4.D0)
     *         +ENREAR
   30 CONTINUE
      ENREAR=12.5660*HBARC*ENREAR
C
      IF(IZPE.EQ.0) THEN
       ZPE=3.D0/4.D0*41.D0*(ABARY+NLAMDA)**(-1.D0/3.D0)
      ENDIF
      IF(IZPE.EQ.1) THEN
       CALL ZPECOR
      ENDIF
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C     ADD LAMDA
C
      SUMPONL=0.0
      SUMPONSS=0.0
C
      IF (NLAMDA.EQ.0) GOTO 500
C
      DO 25 I=2,NGRID
      RX=RSTEP*(I-1)
      SUMPONL=-WGTINT(I)*RX**2 *(
     *       SPOTDL(I)*SRHOL(I)+VPOTL(I)*VRHOL(I)+COULPT(I)*CRHOL(I)
     *       +(VPOTL(I)/VMESCPL)*VMESTCL/YMSL(1)*TERHOL(I)
     *       ) +SUMPONL
   25 CONTINUE
      SUMPONL=12.5660*SUMPONL*0.5D0*HBARC
C
C     ADD strange mesons
C
      DO 50 I=2,NGRID
      RX=RSTEP*(I-1)
      SUMPONSS=-WGTINT(I)*RX**2*(SSPOTL(I)*SRHOL(I)+VSPOTL(I)*VRHOL(I))
     *	     +SUMPONSS
   50 CONTINUE
      SUMPONSS=12.5660*SUMPONSS*0.5D0*HBARC
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  500 CONTINUE
C
      SUMPON=SUMPON+SUMPONL
C
      ESHELL=SUMSPN+SUMPON
      ENRG=ESHELL+EPAIR+ENREAR-ZPE+NLAMDA*SPENRGL+SUMPONSS
C
      RETURN
      END
C
C---   HPSI  -----------------------------------------------------------
C
      SUBROUTINE HPSI(WFIN,HWFOUT,XKAP,PARITY,ISO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'CONDIT.INC'
      INCLUDE 'CONSTA.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'YMSTAR.INC'
      DIMENSION WFIN(NGRDM2),HWFOUT(NGRDM2),AUX1(NGRDM2)
C
C
      CALL DERIV1(WFIN,AUX1,PARITY)
C
      IF(ISO.EQ.1) THEN
      HWFOUT(1)=(XKAP-1.D0)*AUX1(NGRIDP)+
     *          (BARMSE(1)+SPOTM(1)+VPOT(1)
     *          -DPOT(1)-RPOT(1)+COULPT(1))*WFIN(1)
     *          -D1VPOT(1)*WFIN(NGRIDP)
      HWFOUT(NGRIDP)=(XKAP+1.D0)*AUX1(1)+
     *               (-BARMSE(1)-SPOTM(1)+VPOT(1)
     *                +DPOT(1)-RPOT(1)+COULPT(1))*WFIN(NGRIDP)
     *                -D1VPOT(1)*WFIN(1)
C
      DO 10 I=2,NGRID
      XKPHBR=XKAP/(I*RSTEP-RSTEP)
      HWFOUT(I)=(XKPHBR-D1VPOT(I))*WFIN(I+NGRID)-AUX1(I+NGRID)+
     *          (BARMSE(1)+SPOTM(I)+VPOT(I)
     *          -DPOT(I)-RPOT(I)+COULPT(I))*WFIN(I)
      HWFOUT(I+NGRID)=(XKPHBR-D1VPOT(I))*WFIN(I)+AUX1(I)+
     *                (-BARMSE(1)-SPOTM(I)+VPOT(I)
     *                +DPOT(I)-RPOT(I)+COULPT(I))*WFIN(I+NGRID)
   10 CONTINUE
C
      ELSE
      HWFOUT(1)=(XKAP-1.00)*AUX1(NGRIDP)+
     *          (BARMSE(2)+SPOTM(1)+VPOT(1)
     *          +DPOT(1)+RPOT(1))*WFIN(1)
     *          -D1VPOT(1)*WFIN(NGRIDP)
      HWFOUT(NGRIDP)=(XKAP+1.00)*AUX1(1)+
     *               (-BARMSE(2)-SPOTM(1)+VPOT(1)
     *                -DPOT(1)+RPOT(1))*WFIN(NGRIDP)
     *                -D1VPOT(1)*WFIN(1)
C
      DO 20 I=2,NGRID
      XKPHBR=XKAP/(I*RSTEP-RSTEP)
      HWFOUT(I)=(XKPHBR-D1VPOT(I))*WFIN(I+NGRID)-AUX1(I+NGRID)+
     *          (BARMSE(2)+SPOTM(I)+VPOT(I)
     *          +DPOT(I)+RPOT(I))*WFIN(I)
      HWFOUT(I+NGRID)=(XKPHBR-D1VPOT(I))*WFIN(I)+AUX1(I)+
     *                (-BARMSE(2)-SPOTM(I)+VPOT(I)
     *                 -DPOT(I)+RPOT(I))*WFIN(I+NGRID)
   20 CONTINUE
      ENDIF
C
      RETURN
      END
C
CSSS   HPSIL  -----------------------------------------------------------
C
      SUBROUTINE HPSIL(WFIN,HWFOUT,XKAP,PARITY,ISO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'CONDIT.INC'
      INCLUDE 'CONSTA.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'YMSTAR.INC'
      INCLUDE 'YLAMDA.INC'
      INCLUDE 'LAMDA.INC'
      INCLUDE 'YSSMES.INC'
      DIMENSION WFIN(NGRDM2),HWFOUT(NGRDM2),AUX1(NGRDM2)
C
C
      CALL DERIV1(WFIN,AUX1,PARITY)
C
      HWFOUT(1)=(XKAP-1.00)*AUX1(NGRIDP)+
     *          (YMSL(1)+SPOTML(1)+VPOTL(1)+SSPOTL(1)+VSPOTL(1)
     *           +COULPT(1))*WFIN(1)
     *          -D1VPOTL(1)*WFIN(NGRIDP)
      HWFOUT(NGRIDP)=(XKAP+1.00)*AUX1(1)+
     *    (-YMSL(1)-SPOTML(1)+VPOTL(1)-SSPOTL(1)+VSPOTL(1)
     *    +COULPT(1))*WFIN(NGRIDP)
     *    -D1VPOTL(1)*WFIN(1)
C
      DO 20 I=2,NGRID
      XKPHBR=XKAP/(I*RSTEP-RSTEP)
       HWFOUT(I)=(XKPHBR-D1VPOTL(I))*WFIN(I+NGRID)-AUX1(I+NGRID)+
     *          (YMSL(1)+SPOTML(I)+VPOTL(I)+SSPOTL(I)+VSPOTL(I)
     *           +COULPT(I))*WFIN(I)
      HWFOUT(I+NGRID)=(XKPHBR-D1VPOTL(I))*WFIN(I)+AUX1(I)+
     *          (-YMSL(1)-SPOTML(I)+VPOTL(I)-SSPOTL(I)+VSPOTL(I)
     *           +COULPT(I))*WFIN(I+NGRID)
   20 CONTINUE
C
      RETURN
      END
C      
C---   INITGR   --------------------------------------------------------
C
      SUBROUTINE INITGR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C      
      INCLUDE 'CONDIT.INC'
      INCLUDE 'CONSTA.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'SPECIE.INC'
      INCLUDE 'WAVEFN.INC'
      DIMENSION NSP(005,0019),C(4),CAN(13),NMAXFL(2,11)
C
      DATA EZW,VDR,FSX,FZWT,ZDR/.83333333E-1,1.3333333,.83333333,
     *                          2.0833333,.66666667/
      DATA NMAXFL/1,2, 3,6, 3,8, 6,20, 7,28, 11,50, 16,82, 22,126,
     *             29,184, 37,258, 46,350/
      DATA NSP/1,6,15,28,45,           3,10,21,36,99,
     *         2,9,20,35,99,           5,14,27,44,99,
     *         4,13,26,43,99,          8,19,34,99,99,
     *         7,18,33,99,99,          12,25,42,99,99,
     *         11,24,41,99,99,            17,32, 3*99,
     *         16,31,    3*99,           23,40,    3*99,
     *         22,39,    3*99,           30,       4*99,
     *         29,       4*99,           38,       4*99,
     *         37,       4*99,                     5*99,
     *         46,       4*99/
      DATA C,CAN/.62222222,1.4222222,.53333333,1.422222222,
     *           .68611111,1.125,1.125,.375,.64444444,1.33333333,
     *           .70833333,1.125,1.125,.375,.64444444,1.33333333,
     *           .33333333/
C
C
      IF(IPAIR.EQ.1) THEN
       DO 25 I=1,11
       ISAV=I
       NSHLPR=NMAXFL(1,I)
       IF(NPROT.LE.NMAXFL(2,I)) GO TO 30
   25  CONTINUE
C
   30  CONTINUE
       DO 35 I=1,11
       ISAV=I
       NSHLNE=NMAXFL(1,I)
       IF(NNEUT.LE.NMAXFL(2,I)) GO TO 100
   35  CONTINUE
C
      ELSE
       NSHLPR=INSHEL(ITURN,1)
       NSHLNE=INSHEL(ITURN,2)
      ENDIF
C
  100 CONTINUE
C 
      WGTINT(1)=.5*C(1)*RSTEP
      DO 110 I=2,NGRID
      N=MOD(I-1,4)+1
      WGTINT(I)=C(N)*RSTEP
  110 CONTINUE
      NSW=MOD(NGRID,4)
      IF(NSW.EQ.1) THEN
       WGTINT(NGRID)=.5*WGTINT(NGRID)
      ELSEIF(NSW.EQ.2) THEN
       IA=NGRID-6
       DO 120 I=IA,NGRID
       WGTINT(I)=CAN(I-IA+5)*RSTEP
  120  CONTINUE
      ELSEIF(NSW.EQ.3) THEN
       IA=NGRID-3
       DO 130 I=IA,NGRID
       WGTINT(I)=CAN(I-IA+11)*RSTEP
  130  CONTINUE
      ELSE
       IA=NGRID-4
       DO 140 I=IA,NGRID
       WGTINT(I)=CAN(I-IA+1)*RSTEP
  140  CONTINUE
      ENDIF
      DO 150 I=1,NGRID
      WGTINT(I+NGRID)=WGTINT(I)
  150 CONTINUE
C
      NGRID2=NGRID+NGRID
      NGRIDP=NGRID+1
C
      EZWH=EZW/RSTEP
      ESXH=EZWH+EZWH
      VDRH=VDR/RSTEP
      FSXH=FSX/RSTEP
      FZWTH=FZWT/RSTEP
      ZDRH=ZDR/RSTEP
      EVRH=0.25D0/RSTEP
      VRH=4.D0/RSTEP
      DRH=3.D0/RSTEP
      DR2H=1.5D0/RSTEP
      EHH=0.5D0/RSTEP
C
      APROT=FLOAT(NPROT)
      ANEUT=FLOAT(NNEUT)
      ABARY=FLOAT(NPROT+NNEUT)
C
      NSHLTO=NSHLPR+NSHLNE
      ENINV0=ENINV/HBARC
C
       ZENTRF(1)=0.0
       DO 210 I=2,NGRID
       RA=I*RSTEP-RSTEP
       ZENTRF(I)=1.0/(RA*RA)
  210  CONTINUE
C
C
      DO 410 ISO=1,2
      DO 410 JS=1,JSMXMX
      DO 410 NSHL=1,NMAXMX
  410 NPLACE(NSHL,JS,ISO)=0.
C
      DO 420 ISO=1,2
      DO 420 JS=1,JSMXMX
  420 NSMAX(JS,ISO)=0
C
      DO 430 ISO=1,2
      IF(ISO.EQ.1) THEN
       NSHLCT=NSHLPR
      ELSE
       NSHLCT=NSHLNE
      ENDIF
      DO 440 JS=1,19
      DO 450 N=1,5
      IF(NSHLCT-NSP(N,JS).LT.0) GOTO 460
  450 CONTINUE
  460 CONTINUE
      NSMAX(JS,ISO)=N-1
  440 CONTINUE
      DO 470 JS=1,19
      IF(NSMAX(JS,ISO).NE.0) JSMAX(ISO)=JS
  470 CONTINUE
  430 CONTINUE
C
      NPL=1
      DO 480 ISO=1,2
      DO 480 JS=1,JSMAX(ISO)
      DO 480 NSHL=1,NSMAX(JS,ISO)
      NPLACE(NSHL,JS,ISO)=NPL
      NPL=NGRID2+NPL
  480 CONTINUE
C
      NPL=1
      DO 510 ISO=1,2
      DO 510 JS=1,JSMAX(ISO)
      DO 510 NSHL=1,NSMAX(JS,ISO)
      XMULT(NPL)=(XJ(JS)+0.5D0)
      NPL=1+NPL
  510 CONTINUE
C
      RETURN
      END
C      
C---   INITWF   --------------------------------------------------------
C
      SUBROUTINE INITWF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C      
      INCLUDE 'CONDIT.INC'
      INCLUDE 'CONSTA.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'SPECIE.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'YMSTAR.INC'
      INCLUDE 'YLAMDA.INC'
      INCLUDE 'LAMDA.INC'
      INCLUDE 'YSSMES.INC'
      DIMENSION NSP(5,19),NZAHL(46),NMAG(46),NPLORD(NSHLMX,2)
      DIMENSION IIPAIR(2)
C
      DATA NZAHL/2,6,8,14,18,20,28,34,38,40,50,58,64,68,70,82,92,100,
     *       106,110,112,126,138,148,156,162,166,168,184,
     *       198,210,220,228,234,238,240,258,274,288,300,310,318,
     *       324,328,330,350/
      DATA NMAG/101,202,303,3*406,707,4*811,5*1216,6*1722,7*2329,
     *          8*3037, 9*3846/
      DATA NSP/1,6,15,28,45,           3,10,21,36,99,
     *         2,9,20,35,99,           5,14,27,44,99,
     *         4,13,26,43,99,          8,19,34,99,99,
     *         7,18,33,99,99,         12,25,42,99,99,
     *         11,24,41,99,99,         17,32,    3*99,
     *         16,31,    3*99,         23,40,    3*99,
     *         22,39,    3*99,         30,       4*99,
     *         29,       4*99,         38,       4*99,
     *         37,       4*99,                   5*99,
     *         46,       4*99/
C
C
      DO 10 I=1,NGRID
      SRHO(I)=0.D0
      CRHO(I)=0.D0
      VRHO(I)=0.D0
      DRHO(I)=0.D0
      RRHO(I)=0.D0
      SPOTM(I)=0.D0
      SPOTD(I)=0.D0
      VPOT(I)=0.D0
      DPOT(I)=0.D0
      RPOT(I)=0.D0
      COULPT(I)=0.D0
      TERHO(I)=0.D0
      D1VPOT(I)=0.D0
C
      SRHOL(I)=0.D0
      VRHOL(I)=0.D0
      CRHOL(I)=0.D0
      SPOTML(I)=0.D0
      SPOTDL(I)=0.D0
      VPOTL(I)=0.D0
	TERHOL(I)=0.D0
	D1VPOTL(I)=0.D0
   10 CONTINUE
C
      DO 20 I=1,NGRID
      SMESON(I)=0.D0
      VMESON(I)=0.D0
      DMESON(I)=0.D0
      RMESON(I)=0.D0
   20 CONTINUE
C
      CALL DCOUPL
C
      DO 30 N=1,NSHLTO
      SPENRG(N)=0.D0
   30 CONTINUE
C
C     Initial wave functions
C
      DMPWID=0.8D0
      WIDTH=DMPWID*0.5*HMC*OMEG/ABARY**(1.D0/3.D0)
      DO 100 ISO=1,2
      DO 100 JS=1,JSMAX(ISO)
      IF(NSMAX(JS,ISO).LE.0) GO TO 100
      DO 110 NSHL=1,NSMAX(JS,ISO)
      NPLM=NPLACE(NSHL,JS,ISO)-1
C
      DO 120 I=NGRIDP,NGRID2
      WFN(I+NPLM)=0.D0
  120 CONTINUE
C
      IF(NSHL.EQ.1) THEN
       IF(JS.EQ.1) THEN
        DO 130 I=1,NGRID
        RX=RSTEP*(I-1)
        WFN(I+NPLM)=EXP(-WIDTH*RX**2)*RX
  130   CONTINUE
       ELSE
        XLA=XL(JS)
        DO 140 I=1,NGRID
        RX=RSTEP*(I-1)
        WFN(I+NPLM)=WFN(I)*RX**XLA
  140   CONTINUE
       ENDIF
C
      ELSE
       DO 150 I=1,NGRID
       RX=RSTEP*(I-1)
       WFN(I+NPLM)=WFN(I+NPLM-NGRID2)*RX**2
  150  CONTINUE
      ENDIF
  110 CONTINUE
C
      CALL ORTHOG(NPLACE(1,JS,ISO)-1,NSMAX(JS,ISO))
  100 CONTINUE
C
C     Initial occupation probabilities
C
      IF(IPAIR.NE.0) THEN
      IF(IPAIR.EQ.3) THEN
       IIPAIR(1)=IFPAIR(1)
       IIPAIR(2)=IFPAIR(2)
      ENDIF
      NPL=0
      DO 200 ISO=1,2
      DO 200 JS=1,JSMAX(ISO)
      DO 200 NSHL=1,NSMAX(JS,ISO)
      NORD=NSP(NSHL,JS)
      NPL=1+NPL
      NPLORD(NORD,ISO)=NPL
  200 CONTINUE
C
      DO 210 N=1,NSHLTO
      WEIGHT(N)=0.D0
  210 CONTINUE
C
      DO 220 ISO=1,2
      IF(ISO.EQ.1) THEN
       NZ=NPROT
       IE=NSHLPR
      ELSE
       NZ=NNEUT
       IE=NSHLNE
      ENDIF
C
      DO 230 N=1,IE
      IF(NZAHL(N).GE.NZ) GO TO 240
      NPL=NPLORD(N,ISO)
      WEIGHT(NPL)=1.0
  230 CONTINUE
      STOP '  NZ WAS TOO LARGE !'
C
  240 CONTINUE
      IA=NMAG(N)
      IB=MOD(IA,100)
      IA=IA/100
CSSSSSSSSSSSSSS
      IF (IA.EQ.1) THEN
	NOLD=0
	ELSE
      NOLD=NZAHL(IA-1)
	ENDIF
CSSSSSSSSSSSSSS
      X=FLOAT(NZ-NOLD)/FLOAT(NZAHL(IB)-NOLD)
      DO 250 N=IA,IB
      NPL=NPLORD(N,ISO)
      IF(X.GT.1.0) X=1.0
      IF(X.LT.0.D0) X=0.D0 
  250 WEIGHT(NPL)=X
      IF(X.GT. 0.99999 .AND. IB.EQ.IE) THEN
       IFPAIR(ISO)=0
      ELSE
       IFPAIR(ISO)=1
      ENDIF
  220 CONTINUE
C
      ELSE
      DO 300 I=1,NSHLPR
      WEIGHT(I)=WEIGIN(ITURN,I)
  300 CONTINUE      
      DO 310 I=NSHLPR+1,NSHLPR+NSHLNE
      WEIGHT(I)=WEIGIN(ITURN,I)
  310 CONTINUE
      ENDIF
C
      IF(IPAIR.EQ.3  .AND.  IIPAIR(1).EQ.0) THEN
      DO 350 I=1,NSHLPR
      WEIGHT(I)=WEIGIN(ITURN,I)
  350 CONTINUE
      IFPAIR(1)=IIPAIR(1)
      ENDIF
      IF(IPAIR.EQ.3  .AND.  IIPAIR(2).EQ.0) THEN
      DO 360 I=NSHLPR+1,NSHLPR+NSHLNE
      WEIGHT(I)=WEIGIN(ITURN,I)
  360 CONTINUE
      IFPAIR(2)=IIPAIR(2)
      ENDIF
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C     ADD LAMDA
C
      SPENRGL=0.D0
      JS=JSL
      ISO=2
      NSHL=1
      NPLM=NPLACE(NSHL,JS,ISO)-1
C
      WFNL(1)=0.0
      DO 510 I=2,NGRID2
      WFNL(I)=WFN(I+NPLM)
  510 CONTINUE
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C     ADD strange mesons
C
      DO 610 I=1,NGRID
      SSMESON(I)=0.D0
      VSMESON(I)=0.D0
C
      SSPOTL(I)=0.D0
      VSPOTL(I)=0.D0
  610 CONTINUE
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      RETURN
      END
C      
C---   LAPLAB   --------------------------------------------------------
C
      SUBROUTINE LAPLAB(R,RHS,ESUB,RSTEP,NMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'WRKINV.INC'
      DIMENSION R(NGRDMX),RHS(NGRDMX)
C      
      DATA FSX,EZW/.83333333,.83333333E-1/
C
C
       DO 17 I=1,NMAX
   17  R(I)   = RHS(I)/SQB(I)
C
        OFFDIG = 1./(RSTEP*RSTEP)
        OFFDG2 = OFFDIG+OFFDIG
        NMAXM  = NMAX-1
C
        BINV(1) = 0.0
        AINV(1) = 1.0
        X      = 0.
        Y      = 1.0
        RIM1   = R(1)
C
        DO 13 I=2,NMAXM
C                                         FORWARD ELIM. FOR R.H.S.
        X      = EZW*(R(I+1)+RIM1)+FSX*R(I)- X*(BINV(I-1)/Y)
        RIM1   = R(I)
        R(I)   = X
C                                         FORWARD ELIM. FOR MATRIX
        UI     = (POTL(I)-ESUB)/(SQB(I)*SQB(I)) + D2SQB(I)
        Z      = -OFFDIG+EZW*UI
        BINV(I) = Z
        Y      = OFFDG2+FSX*UI-Z*(BINV(I-1)/Y)
        AINV(I) = Y
   13   CONTINUE
C
        R(NMAX) = EZW*RIM1+FSX*R(NMAX)-X*(BINV(NMAXM)/Y)
        UI     = (POTL(NMAX)-ESUB)/(SQB(NMAX)*SQB(NMAX)) + D2SQB(NMAX)
        Z      = -OFFDIG+EZW*UI
        BINV(NMAX) = Z
        Y      = OFFDG2+FSX*UI-Z*(BINV(NMAXM)/Y)
        AINV(NMAX) = Y
C
       X      = R(NMAX)/AINV(NMAX)
       R(NMAX) = X/SQB(NMAX)
       DO 15 I=NMAXM,2,(-1)
       X      = (R(I)-X*BINV(I+1))/AINV(I)
       R(I)   = X/SQB(I)
   15  CONTINUE
       R(1)   = 0.
C
       RETURN
       END
C       
C---   LAPLAC   --------------------------------------------------------
C
      SUBROUTINE LAPLAC(R,RHS,XMAS2,RSTEP,NMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'WRKINV.INC'
      DIMENSION R(NGRDMX),RHS(NGRDMX)
C      
      DATA FSX,EZW/.83333333,.83333333E-1/
C      
C
        RINV2  = 1./(RSTEP*RSTEP)
        OFFDIG = -RINV2+EZW*XMAS2
        OFFDG2 = OFFDIG*OFFDIG
        ONDIAG = RINV2+RINV2+FSX*XMAS2
        NMAXM  = NMAX-1
C
        AINV(1) = 1.0
        X      = 0.
        Y      = 1.0
C
        X      = EZW*(RHS(3)+RHS(1))+FSX*RHS(2)
        R(2)   = X
        Y      = ONDIAG
        AINV(2)   = Y
C                                         NOW THE MOST OTHER I
        DO 13 I=3,NMAXM
C                                         FORWARD ELIM. FOR R.H.S.
        X      = EZW*(RHS(I+1)+RHS(I-1))+FSX*RHS(I)- X*(OFFDIG/Y)
        R(I)   = X
C                                         FORWARD ELIM. FOR MATRIX
        Y      = ONDIAG-OFFDG2/Y
        AINV(I) = Y
   13   CONTINUE
C                                         LAST ELIM. STEP EXPLICITELY
        R(NMAX) = EZW*RHS(NMAXM)+FSX*RHS(NMAX)-X*(OFFDIG/Y)
        Y      = ONDIAG-OFFDG2/Y
        AINV(NMAX) = Y
C
       X      = R(NMAX)/AINV(NMAX)
       R(NMAX) = X
       DO 15 I=NMAXM,2,(-1)
       X      = (R(I)-X*OFFDIG)/AINV(I)
       R(I)   = X
   15  CONTINUE
       R(1)   = 0.
C
       RETURN
       END
C       
C---   LOWCMP   --------------------------------------------------------
C
      SUBROUTINE LOWCMP(WFIN,DWFOUT,XKAP,PARITY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C      
      INCLUDE 'CONSTA.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'WAVEFN.INC'
      DIMENSION WFIN(NGRDMX),DWFOUT(NGRDMX)
C
C
      XKR=XKAP/RSTEP
      IF(PARITY.GT.0.) THEN
        DWFOUT(1)=0.D0
      ELSE
        DWFOUT(1)=(VDRH*WFIN(2)-ESXH*WFIN(3))*(1.00+XKAP)
      ENDIF
      DWFOUT(2)=(ZDRH*(WFIN(3)-WFIN(1))
     *           -EZWH*(WFIN(4)-PARITY*WFIN(2))
     *           +(XKR-D1VPOT(2))*WFIN(2))
      DO 10 I=3,NGRID-2
      DWFOUT(I)=(EZWH*(WFIN(I-2)-WFIN(I+2))
     *          +ZDRH*(WFIN(I+1)-WFIN(I-1))
     *          +(XKR/(I-1)-D1VPOT(I))*WFIN(I) )
   10 CONTINUE
      IA=NGRID-1
      DWFOUT(IA)=(EVRH*WFIN(NGRID)+FSXH*WFIN(NGRID-1)-DR2H*WFIN(NGRID-2)
     *           +EHH*WFIN(NGRID-3)-EZWH*WFIN(NGRID-4)
     *           +(XKR/(IA-1)-D1VPOT(IA))*WFIN(IA) )
      DWFOUT(NGRID)=(FZWTH*WFIN(NGRID)-VRH*WFIN(NGRID-1)
     *              +DRH*WFIN(NGRID-2)-VDRH*WFIN(NGRID-3)
     *              +EVRH*WFIN(NGRID-4)
     *              +(XKR/(NGRID-1)-D1VPOT(NGRID))*WFIN(NGRID) )
C
      RETURN
      END
C
C---   MESON   ---------------------------------------------------------
C
      SUBROUTINE MESON
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'CONDIT.INC'      
      INCLUDE 'CONSTA.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'WORKSP.INC'
      INCLUDE 'YMSTAR.INC'
      INCLUDE 'YLAMDA.INC'
      INCLUDE 'LAMDA.INC'
      INCLUDE 'YSSMES.INC'
      DIMENSION RHS(NGRDMX)
C
C     For scalar meson
C
      RHS(1)=0.0
      DO 10 I=2,NGRID
      RX=RSTEP*(I-1)
      SMS=SMESON(I)
      RHS(I)=-SMESCPD(I)*RX*SRHO(I)
     *      -(SMSB2C/RX+SMSB3C/(RX*RX)*SMS)*SMS*SMS
     *      -SMESCPDL(I)*RX*SRHOL(I)
   10 CONTINUE
C
      CALL LAPLAC(SMESON,RHS,SMASS2,RSTEP,NGRID)
      CALL DCOUPL
C      
C      SPOTD(1)=SMESCPD(1)*(VDRH*SMESON(2)-ESXH*SMESON(3))
C      DO 15 I=2,NGRID
C      RX=RSTEP*(I-1)
C      SPOTD(I)=SMESCPD(I)*SMESON(I)/RX
C   15 CONTINUE
C
C     For vector meson
C
      RHS(1)=0.0
      DO 20 I=2,NGRID
      RX=RSTEP*(I-1)
      VFIELD=VMESON(I)
      RHS(I)=VMESCP*RX*VRHO(I)
     *      -VMCPC3/(RX*RX)*VFIELD**3
     *      +RX*VMESTC/BARMSE(2)*TERHO(I)
     *      +VMESCPL*RX*VRHOL(I)
     *      +RX*VMESTCL/YMSL(1)*TERHOL(I)
   20 CONTINUE
C
      CALL LAPLAC(VMESON,RHS,VMASS2,RSTEP,NGRID)
C      
      VPOT(1)=VMESCP*(VDRH*VMESON(2)-ESXH*VMESON(3))
      VPOTL(1)=VPOT(1)*VMESCPL/VMESCP
      DO 25 I=2,NGRID
      RX=RSTEP*(I-1)
      VPOT(I)=VMESCP*VMESON(I)/RX
      VPOTL(I)=VPOT(I)*VMESCPL/VMESCP
   25 CONTINUE
C
C     Tensor part
C
      FACT=VMESTC/RSTEP
      DO 27 I=2,NGRID
      AUX1(I)=FACT/(2.0*BARMSE(2))*VMESON(I)/(I-1)
   27 CONTINUE
      AUX1(1)=(4.0*AUX1(2)-AUX1(3))*0.33333330
      DO 28 I=NGRIDP,NGRID2
      AUX1(I)=0.0D0
   28 CONTINUE
      CALL DERIV0(AUX1,D1VPOT,ONE)
C
C     LAMBDA Tensor part
C  
      FACT=VMESTCL/RSTEP
      DO 22 I=2,NGRID
      AUX1(I)=FACT/(2.0*YMSL(1))*VMESON(I)/(I-1)
   22 CONTINUE
      AUX1(1)=(4.0*AUX1(2)-AUX1(3))*0.33333330
      DO 23 I=NGRIDP,NGRID2
      AUX1(I)=0.0D0
   23 CONTINUE
      CALL DERIV0(AUX1,D1VPOTL,ONE)
C
C     For iso-vector scalar meson
C
      RHS(1)=0.0
      DO 30 I=2,NGRID
      RX=RSTEP*(I-1)
      DFIELD=DMESON(I)
      RHS(I)=-ISOVEC*DMESCP*RX*DRHO(I)-DMCPD3/(RX*RX)*DFIELD**3
   30 CONTINUE
C
      CALL LAPLAC(DMESON,RHS,DMASS2,RSTEP,NGRID)
C      
      DPOT(1)=DMESCP*(VDRH*DMESON(2)-ESXH*DMESON(3))
      DO 35 I=2,NGRID
      RX=RSTEP*(I-1)
      DPOT(I)=DMESCP*DMESON(I)/RX
   35 CONTINUE
C
C     For iso-vector vector meson
C
      RHS(1)=0.0
      DO 40 I=2,NGRID
      RX=RSTEP*(I-1)
      RFIELD=RMESON(I)
      RHS(I)=ISOVEC*RMESCP*RX*RRHO(I)-RMCPE3/(RX*RX)*RFIELD**3
   40 CONTINUE
C     
      CALL LAPLAC(RMESON,RHS,RMASS2,RSTEP,NGRID)
      RPOT(1)=RMESCP*(VDRH*RMESON(2)-ESXH*RMESON(3))
      DO 45 I=2,NGRID
      RX=RSTEP*(I-1)
      RPOT(I)=RMESCP*RMESON(I)/RX
   45 CONTINUE
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C     ADD strange mesons
C
      IF (NLAMDA.LT.2) GOTO 700
C
C     For strange scalar meson
C
      RHS(1)=0.0
      DO 50 I=2,NGRID
      RX=RSTEP*(I-1)
      RHS(I)=-SSMESCPL*RX*SRHOL(I)
   50 CONTINUE
C
      CALL LAPLAC(SSMESON,RHS,SSMASS2,RSTEP,NGRID)
C      
      SSPOTL(1)=SSMESCPL*(VDRH*SSMESON(2)-ESXH*SSMESON(3))
      DO 55 I=2,NGRID
      RX=RSTEP*(I-1)
      SSPOTL(I)=SSMESCPL*SSMESON(I)/RX
   55 CONTINUE
C
C     For vector meson
C
      RHS(1)=0.0
      DO 60 I=2,NGRID
      RX=RSTEP*(I-1)
      RHS(I)=+VSMESCPL*RX*VRHOL(I)
   60 CONTINUE
C
      CALL LAPLAC(VSMESON,RHS,VSMASS2,RSTEP,NGRID)
C      
      VSPOTL(1)=VSMESCPL*(VDRH*VSMESON(2)-ESXH*VSMESON(3))
      DO 65 I=2,NGRID
      RX=RSTEP*(I-1)
      VSPOTL(I)=VSMESCPL*VSMESON(I)/RX
   65 CONTINUE
C
  700 CONTINUE
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      RETURN
      END
C
C---   ORTHOG  ---------------------------------------------------------
C
      SUBROUTINE ORTHOG(NBASE,NMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C        
      INCLUDE 'GRID.INC'
      INCLUDE 'WAVEFN.INC'
C      
C
      DO 10 N1=1,NMAX
      IPL1=N1*NGRID2-NGRID2+NBASE
C
      IF(N1.GT.1) THEN
       DO 20 N2=1,N1-1
       IPL2=N2*NGRID2-NGRID2+NBASE
       OVLP=0.0
       DO 30 I=2,NGRID2
   30  OVLP=WFN(I+IPL1)*WFN(I+IPL2)*WGTINT(I)+OVLP
       DO 40 I=1,NGRID2
   40  WFN(I+IPL1)=WFN(I+IPL1)-OVLP*WFN(I+IPL2)
   20  CONTINUE
      ENDIF
C
      OVLP=0.0
      DO 50 I=1,NGRID2
   50 OVLP=WFN(I+IPL1)*WFN(I+IPL1)*WGTINT(I)+OVLP
      OVLP=1./SQRT(OVLP)
      DO 60 I=1,NGRID2
   60 WFN(I+IPL1)=WFN(I+IPL1)*OVLP
   10 CONTINUE
C
      RETURN
      END
C      
C---   PAIR   ----------------------------------------------------------
C
      SUBROUTINE PAIR(E,GW,PH,NZ,NMAX,GP,DELTA,ITER,EPS,LAST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C       
       DIMENSION E(NSHLMX),GW(NSHLMX),PH(NSHLMX)
       DIMENSION DEFERM(3),DDELTA(3),DPARN(3),GR(2,2)
C       
       DATA DEFERM,DDELTA,EPSIL/.0,1.,.0,.0,.0,.10,.9/
       DATA GR,DPARN/7*0./
       DATA EFMAX,DETMIN,DEFMIN/3.,.1E-3,.1E-3/
C
C
       IWARNG=0
       DEFERM(1)=0.0D0
       DEFERM(2)=1.0D0
       DDELTA(3)=0.1D0
       ITERA=IABS(ITER)
C
      SLAM=0.0
      ELAM=0.0
      DO 110 IT=1,NMAX
      WGT=SQRT(GW(IT)-GW(IT)*GW(IT))*PH(IT)
      SLAM=WGT      +SLAM
      ELAM=WGT*E(IT)+ELAM
  110 CONTINUE
      ELAM=ELAM/(SLAM+1.E-20)
      EFERM=ELAM
C
      DO 310 IT=1,ITERA
      DO 320 KA=1,2
      K=3-KA
      ELAM=EFERM+DEFERM(K)
      PARNUM=0.0D0
      DO 330 I=1,NMAX
      EQUASI=SQRT((E(I)-ELAM)*(E(I)-ELAM)+DELTA *DELTA )
      GW(I)=.5-.5*(E(I)-ELAM)/EQUASI
      IF(E(I).GT.0.D0) GW(I)=0.D0
      IF(GW(I).LT.0.0) GW(I)=0.0  
  305 PARNUM=GW(I)*PH(I)+PARNUM
  330 CONTINUE
  331 CONTINUE
      DPARN(K)=PARNUM+PARNUM-NZ
  320 CONTINUE
C
      GR(1,1)=(DPARN(2)-DPARN(1))/DEFERM(2)
      DET=SIGN(MAX(ABS(GR(1,1)),DETMIN),GR(1,1))
      EFOLD=EFERM
      EFSTEP=EPSIL*DPARN(1)/DET
      EFERM=EFERM-SIGN(MIN(ABS(EFSTEP),EFMAX),EFSTEP)
      DEFERM(2)=EFOLD-EFERM
      DEFERM(2)=SIGN(MAX(ABS(DEFERM(2)),DEFMIN),DEFERM(2))
      IF(ABS(DPARN(1)).LT.EPS)  GOTO 350
  310 CONTINUE
      IF(ABS(DPARN(1)).GT.EPS*10.E0) IWARNG=1
C
  350 CONTINUE
      GAPEQ=0.0D0
      DO 360 I=1,NMAX
      IF(E(I).GT.0.D0) GO TO 360
      GAPEQ=PH(I)/SQRT((E(I)-ELAM)*(E(I)-ELAM)+DELTA*DELTA)
     *     +GAPEQ
  360 CONTINUE
      IF(ABS(GAPEQ).LT.1.D-20) THEN
       GP=0.D0
      ELSE
       GP=2.0D0/GAPEQ
      ENDIF
C
      IF(LAST.EQ.1.  AND.  IWARNG.EQ.1) THEN
        WRITE(10,*) ' GAP EQ. DOES NOT CONVERGE !!'
      ENDIF

C
       RETURN
       END
C
C---   RMF   -----------------------------------------------------------
C
      SUBROUTINE RMF(ITER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'CONDIT.INC'
      INCLUDE 'CONSTA.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'RESULT.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'YMSTAR.INC'
      INCLUDE 'YLAMDA.INC'
      INCLUDE 'LAMDA.INC'
      INCLUDE 'YSSMES.INC'
      DIMENSION CFRHO(NGRDMX)
      dimension derrho(ngrdmx)
      dimension dspot(ngrdmx),dvpot(ngrdmx),drpot(ngrdmx),dcoul(ngrdmx)
C
C
      IERROR=0
      NGRID=INT((RMAX+1.D-8)/RSTEP)+1
C
      CALL INITGR
      CALL INITWF
C
      EFERMP=0.0
      EFERMN=0.0
      GAPPR=0.0
      GAPNE=0.0
      DELTPR=DELPIN/SQRT(ABARY)
      DELTNE=DELNIN/SQRT(ABARY)
C
      RMSCON=0.D0
C
    5 CONTINUE
      DO 10 ITER=1,ITERMX
      CALL DENSIT(ITER)
      CALL MESON
      CALL COULMB
C
      ENERPR=SPENRG(1)
      ENERNE=SPENRG(NSHLPR+1)
      CALL WFSTEP
CSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      IF (NLAMDA.EQ.0) THEN
      ELSE
      CALL WFSTEPL
      ENDIF
CSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      IF(IERROR.EQ.1) THEN
       IF(IWRIT.GE.1) THEN
        WRITE(10,*) ' NEGATIVE MASS!!'
       ENDIF
       RETURN
      ENDIF
C
      EPAIR=0.0
      LAST=0
      IF(IPAIR.NE.0) THEN
       IF(IFPAIR(1).EQ.1) THEN
        CALL PAIR(SPENRG,WEIGHT,XMULT,NPROT,NSHLPR,GAPPR,DELTPR,
     *            1*ITPAIR,ENDCON,LAST)
        IF(ABS(GAPPR).GT.1.D-8) THEN
          EPAIR=EPAIR-DELTPR*DELTPR/GAPPR
        ENDIF
       ENDIF
C
       IF(IFPAIR(2).EQ.1) THEN
        CALL PAIR(SPENRG(NSHLPR+1),WEIGHT(NSHLPR+1),XMULT(NSHLPR+1),
     *            NNEUT,NSHLNE,GAPNE,DELTNE,1*ITPAIR,ENDCON,LAST)
        IF(ABS(GAPNE).GT.1.D-8) THEN
         EPAIR=EPAIR-DELTNE*DELTNE/GAPNE
        ENDIF
       ENDIF
      ENDIF
C
      ESHLOL=ENRG
      CALL ENERGY
C
      RMSOLD=RMSCON
      RMSCON=0.D0
      DO 20 I=1,NGRID
      RX=RSTEP*(I-1)
      RMSCON=RX**4*WGTINT(I)*CRHO(I)+RMSCON
   20 CONTINUE
C    
      IF(ABS(ENRG-ESHLOL).LT.ENDCON*ABS(ENRG)  .AND.
     *   ABS(RMSCON-RMSOLD).LE.ENDCON*ABS(RMSCON)  .AND.
     *   ABS(ENERPR-SPENRG(1)) .LE. ENDCON*ABS(ENERPR)  .AND.
     *   ABS(ENERNE-SPENRG(NSHLPR+1)) .LE.
     *   ENDCON*ABS(ENERNE)  ) GO TO 30
C
   10 CONTINUE
C
      IF(IWRIT.GE.1) THEN
       WRITE(10,'(1(A,I5))') ' ITERMX EXHAUSTED AT ITER=',ITER
      ENDIF
      IERROR=1
      RETURN
C
   30  CONTINUE
       IF(IWRIT.GE.1) THEN
        WRITE(10,'(A,F7.3)') 'RMAX=',RSTEP*(NGRID-1)
        WRITE(10,'(A,I5)') ' END OF ITERATION AT ITER=',ITER
        WRITE(10,'(A,F22.15)') ' ENERGY1=',ENRG
       ENDIF
C
      CALL DENSIT(1)
      CALL MESON
      CALL COULMB
      CALL WFSTEP
CSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      IF (NLAMDA.EQ.0) THEN
      ELSE
      CALL WFSTEPL
      ENDIF
CSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      EPAIR=0.0
      IF(IPAIR.NE.0) THEN
       IF(IFPAIR(1).EQ.1) THEN
        LAST=1
        CALL PAIR(SPENRG,WEIGHT,XMULT,NPROT,NSHLPR,GAPPR,DELTPR,
     *            1*ITPAIR,ENDCON,LAST)
         IF(ABS(GAPPR).LT.1.D-8) THEN
          WRITE(10,*) ' THERE IS SOME TROUBLE IN PARING CALCULATION!!'
         ENDIF
        EPAIR=EPAIR-DELTPR*DELTPR/GAPPR
       ENDIF
C    
       IF(IFPAIR(2).EQ.1) THEN
        LAST=1
        CALL PAIR(SPENRG(NSHLPR+1),WEIGHT(NSHLPR+1),XMULT(NSHLPR+1),
     *            NNEUT,NSHLNE,GAPNE,DELTNE,1*ITPAIR,ENDCON,LAST)
         IF(ABS(GAPNE).LT.1.D-8) THEN
          WRITE(10,*) ' THERE IS SOME TROUBLE IN PARING CALCULATION!!'
         ENDIF        
        EPAIR=EPAIR-DELTNE*DELTNE/GAPNE
       ENDIF
      ENDIF
C
      CALL ENERGY
C
      IF(IWRIT.GE.1) THEN
       WRITE(10,'(A,F22.15)') ' ENERGY2=',ENRG
      ENDIF
C
C   
      SUMP=0.D0
      SUMN=0.D0
      SUMM=0.D0
      SUML=0.D0
      X=0.D0
      DO 210 I=2,NGRID
      X=RSTEP+X
      XX=(X*X)*(X*X)*WGTINT(I)
      SUMP=XX*CRHO(I)+SUMP
      SUMN=XX*(VRHO(I)-CRHO(I))+SUMN
      SUMM=XX*(VRHO(I)+VRHOL(I))+SUMM
      SUML=XX*(VRHOL(I))+SUML
  210 CONTINUE
      RMSP=4.D0*PI*SUMP/APROT
      RMSN=4.D0*PI*SUMN/ANEUT
      RMSM=4.D0*PI*SUMM/(ABARY+NLAMDA)
      IF (NLAMDA.EQ.0) THEN
	RMSL=0.0D0
	ELSE
	RMSL=4.D0*PI*SUML/NLAMDA
	ENDIF
C
c      IF(IZPE.EQ.0) THEN
c       A1=1.519D0*ABARY**(-2.D0/3.D0)
c       RMSP=RMSP-A1
c       RMSN=RMSN-A1
c       RMSM=RMSM-A1
c      ENDIF
c      IF(IZPE.EQ.1) THEN
c       RMSP=RMSP-2.D0*RMSP/ABARY+RMSM/ABARY
c       RMSN=RMSN-2.D0*RMSN/ABARY+RMSM/ABARY
c       RMSM=RMSM-RMSM/ABARY
c      ENDIF
C
      RMSC=RMSP+0.862D0*0.862D0-0.3359D0*0.3359D0*ANEUT/APROT
      RMSC=DSQRT(RMSC)
      RMSP=DSQRT(RMSP)
      RMSN=DSQRT(RMSN)
      RMSM=DSQRT(RMSM)
      RMSL=DSQRT(RMSL)
C
C     Print the final results
C
      IF(IWRIT.GE.1) THEN
       WRITE(10,'(/,A,/,2(A,F11.4),/,3(A,F11.4))')
     *           '  ENERGIES (IN MEV):',
     *           '  TOTAL=',ENRG,'     PER PAR=',ENRG/(ABARY+NLAMDA),
     *           '    ZPE=',ZPE,'   PAIR-PART=',EPAIR
       WRITE(10,'(3(A,F11.4))')
     *           '  SUMSPN=',SUMSPN,'   SUMPON=',SUMPON,
     *           '  ENREAR=',ENREAR
       WRITE(10,'(2(A,F11.4))')
     *           '  SUMPONSS=',SUMPONSS,
     *           '  EN(no L)=',ENRG-NLAMDA*SPENRGL
       WRITE(10,'(/,14X,2A)') '         CHARG    PROTON   ',
     *                       ' NEUTRON   MATTER   '
       WRITE(10,'(A,4(F10.4))')
     *              'RMS-RADII (IN FM):',RMSC,RMSP,RMSN,RMSM
       WRITE(10,'(/,2(A,I6),2(A,F10.4))')
     *           '  LAMDA=',NLAMDA,'  JS=',JSL,
     *           '  SPE=',SPENRGL,'  RMS=',RMSL
       CALL PRISPE
      ENDIF
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUMCORE=4.D0*PI*(SUMP+SUMN)/ABARY
      RMSCORE=DSQRT(SUMCORE)
	ECORE  =ENRG-NLAMDA*SPENRGL
C
      WRITE(11,'(4(I5),1(F6.2),3(F12.4))') NLAMDA,NPROT,NNEUT,JSL
     *           ,XJ(JSL),ENRG,ECORE,SPENRGL
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
c
c      write(41,'(2i5,5f14.6)') nprot,nneut,-enrg/abary,
c     *                         rmsc,rmsp,rmsn,rmsm
c
      call deriv0(vrho,derrho,-one)
      srmax=rmsm*sqrt(5.0/3.0)
      nrmax=INT((srmax+1.d-8)/rstep)+1
      do 111 i=nrmax,1,-1
      if(derrho(i) .ge. 0.0) then
        kmax=i
        go to 112
      endif
  111 continue
  112 continue
      cdens=0.0
      do 113 i=1,kmax
      cdens=cdens+vrho(i)
  113 continue
      cdens=cdens/kmax
      write(41,'(2i5,f14.6,f8.3)') nprot,nneut,cdens,rstep*(kmax-1)
c
c
C
      IF(IWRIT.GE.2) THEN
       OPEN(UNIT=50,FILE='densit250.rmf',STATUS='UNKNOWN')
       WRITE(50,'(A,/,2(A,I3),/,A)')
     *    '++++++++++++++++++++++++++++++++++++++++++++++++++++++',
     *       'NPROT=',NPROT,'   NNEUT=',NNEUT,
     *    '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       UCONST=(0.71D0)**(0.5D0)*10.0D0**3/HBARC
       ACONST=-UCONST**2/4.0D0
       DO 310 J=2,NGRID
       AR=J*RSTEP-RSTEP
       CHARDE=0.0
       DO 320 I=1,NGRID
       FR=I*RSTEP-RSTEP
       SUMR=AR+FR
       SUBR=ABS(AR-FR)
       DIVR=FR/AR
       CHARDE=WGTINT(I)*ACONST*DIVR*(CRHO(I)+CRHOL(I))
     *       *((SUMR+1.0/UCONST)*EXP(-UCONST*SUMR)
     *       - (SUBR+1.0/UCONST)*EXP(-UCONST*SUBR))
     *       +CHARDE
  320  CONTINUE
       CFRHO(J)=CHARDE
  310  CONTINUE
       CFRHO(1)=0.0
       DO 330 I=1,NGRID
       FR=I*RSTEP-RSTEP
       CFRHO(1)=UCONST**3/2.0*WGTINT(I)*FR**2
     *         *(CRHO(I)+CRHOL(I))*EXP(-UCONST*FR)+CFRHO(1)
  330  CONTINUE
C
       WRITE(50,*)  '   R[fm]         SRHO            VROH 
     *         DRHO        RRHO [fm^-3]'
       DO 370 I=1,NGRID
       R=RSTEP*(I-1)
       WRITE(50,'(F8.3,4E16.6)') R,SRHO(I),VRHO(I),DRHO(I),RRHO(I)
  370  CONTINUE
       WRITE(50,*)  '   R[fm]
     *         PROTON    NEUTRON   CHARGE  DENSITY[fm^-3]'
       DO 380 I=1,NGRID
       R=RSTEP*(I-1)
       WRITE(50,'(F8.3,3E20.8)') R,CRHO(I),VRHO(I)-CRHO(I),CFRHO(I)
  380  CONTINUE
C
      ENDIF
C
      IF(IWRIT.GE.3) THEN
       OPEN(UNIT=51,FILE='wfnmes250.rmf',STATUS='UNKNOWN')
       WRITE(51,'(A,/,2(A,I3),/,A)')
     *    '++++++++++++++++++++++++++++++++++++++++++++++++++++++',
     *       'NPROT=',NPROT,'   NNEUT=',NNEUT,
     *    '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
       WRITE(51,*) '  R   SIGMA   OMEGA   DELTA    RHO    SIGMA*    PHI'
       DO 400 I=1,NGRID
C
       S1=SPOTD(I)/SMESCPD(I)*HBARC
       W1=VPOT(I)/VMESCP*HBARC
       IF(ABS(DMESCP).LT.1.E-12) THEN
        D1=0.D0
       ELSE
        D1=DPOT(I)/DMESCP*HBARC
       ENDIF
       IF(ABS(RMESCP).LT.1.E-12) THEN
        R1=0.D0
       ELSE
        R1=RPOT(I)/RMESCP*HBARC
       ENDIF
       IF(ABS(SSMESCPL).LT.1.E-12) THEN
        SS1=0.D0
       ELSE
        SS1=SSPOTL(I)/SSMESCPL*HBARC
       ENDIF
       IF(ABS(VSMESCPL).LT.1.E-12) THEN
        VS1=0.D0
       ELSE
        VS1=VSPOTL(I)/VSMESCPL*HBARC
       ENDIF
C
       R=RSTEP*(I-1)
       WRITE(51,'(F8.3,6E12.3)') R,S1,W1,D1,R1,SS1,VS1
  400  CONTINUE
C
       OPEN(UNIT=52,FILE='wfnbar250.rmf',STATUS='UNKNOWN')
       WRITE(52,'(A)') 'BARYON - WAVEFUNCTION '
       WRITE(52,'(A,/,2(A,I3),/,A)')
     *    '++++++++++++++++++++++++++++++++++++++++++++++++++++++',
     *       'NPROT=',NPROT,'   NNEUT=',NNEUT,
     *    '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
C
       NPLORD=0
       DO 410 ISO=1,2
       DO 410 JS=1,JSMAX(ISO)
       DO 410 NS=1,NSMAX(JS,ISO)
         NPL=NPLACE(NS,JS,ISO)
         WRITE(52,'(A,2I5,2F10.3)')
     *         'ISO  N L J=',ISO,NS,XL(JS),XJ(JS)
         WRITE(52,*) 'UPPER COMPONENT'
         WRITE(52,'(F7.4,E20.8)')
     *        (RSTEP*(I-NPL),WFN(I),I=NPL,NPL+NGRID-1)
         WRITE(52,*) 'LOWER COMPONENT'
         WRITE(52,'(F7.4,E20.8)')
     *        (RSTEP*(I-NPL-NGRID),WFN(I),
     *         I=NPL+NGRID,NPL+NGRID2-1)
         NPLORD=NPLORD+1
  410  CONTINUE
      ENDIF
C
C     Check the convergence
C
      TPROT=0.D0
      TNEUT=0.D0
      TNUCL=0.D0
      DO 610 I=2,NGRID
      R=RSTEP*(I-1)
      TPROT=WGTINT(I)*R**2*CRHO(I)+TPROT
      TNEUT=WGTINT(I)*R**2*(VRHO(I)-CRHO(I))+TNEUT
      TNUCL=WGTINT(I)*R**2*VRHO(I)+TNUCL
  610 CONTINUE
      TPROT=4.D0*PI*TPROT
      TNEUT=4.D0*PI*TNEUT
      TNUCL=4.D0*PI*TNUCL
      IF(IWRIT.GE.1) THEN
       WRITE(10,'(3(A,F10.4))') 'PROTON=',TPROT,'    NEUTRON=',TNEUT,
     *                          '    NUCLEON=',TNUCL
      ENDIF
C
c
       write(42,'(a,/,2(a,i3),/,a)')
     *    '++++++++++++++++++++++++++++++++++++++++++++++++++++++',
     *       'nprot=',nprot,'   nneut=',nneut,
     *    '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      call deriv0(SPOTM,dspot,one)
      call deriv0(vpot,dvpot,one)
      call deriv0(rpot,drpot,one)
      call deriv0(coulpt,dcoul,one)
      do 333 i=1,ngrid
      rx=rstep*(i-1)
      write(42,'(f8.2,4e17.6)')
     *  rx,dspot(i)*hbarc,dvpot(i)*hbarc,drpot(i)*hbarc,
     *  dcoul(i)*hbarc
  333 continue
      write(42,*)
      do 334 i=1,ngrid
      rx=rstep*(i-1)
      write(42,'(f8.2,4e17.6)')
     *  rx,SPOTM(i)*hbarc,vpot(i)*hbarc,rpot(i)*hbarc,
     *  coulpt(i)*hbarc
  334 continue
c
      RETURN
      END
C
C---   RMFMN   ---------------------------------------------------------
C
      SUBROUTINE RMFMN(KERROR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C      
      INCLUDE 'BBCOUP.INC'
      INCLUDE 'CONDIT.INC'
      INCLUDE 'CONSTA.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'RESULT.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'SPECIE.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'WORKSP.INC'
      INCLUDE 'YLAMDA.INC'
      INCLUDE 'LAMDA.INC'
C
C
      BARMSE(1)=BARYMS(1)/HBARC
      BARMSE(2)=BARYMS(2)/HBARC
      SMESMS=SMESMS/HBARC
      VMESMS=VMESMS/HBARC
      RMESMS=RMESMS/HBARC
      DMESMS=DMESMS/HBARC
      SMASS2=SMESMS**2
      VMASS2=VMESMS**2
      RMASS2=RMESMS**2
      DMASS2=DMESMS**2
C
      DO 1000 ITURN=1,NNUCL
      NPROT=INPROT(ITURN)
      NNEUT=INNEUT(ITURN)
      RMAX=RMAXIN(ITURN)
	NLAMDA=INLAMDA(ITURN)
	JSL=IJSL(ITURN)
C
      SMESCP=BBCOUP(1,1)+BBCOUP(1,2)/FLOAT(NPROT+NNEUT)**POWER
      VMESCP=BBCOUP(2,1)+BBCOUP(2,2)/FLOAT(NPROT+NNEUT)**POWER
      RMESCP=BBCOUP(3,1)+BBCOUP(3,2)/FLOAT(NPROT+NNEUT)**POWER
      SMSB2C=BBCOUP(4,1)+BBCOUP(4,2)/FLOAT(NPROT+NNEUT)**POWER
      SMSB3C=BBCOUP(5,1)+BBCOUP(5,2)/FLOAT(NPROT+NNEUT)**POWER
      VMCPC3=BBCOUP(6,1)+BBCOUP(6,2)/FLOAT(NPROT+NNEUT)**POWER
      DMESCP=0.D0
      DMCPD3=0.D0
      RMCPE3=0.D0
cc-suga
c      vmestc=-vmescp
      vmestc=0.d0
      VMESCPL=VMESCP*0.7938D0*RGW(ITURN)
c      VMESTCL= 0.0*VMESCPL
      VMESTCL=-1.0*VMESCPL
cc-suga
C
      IF(IWRIT.GE.1) THEN
       WRITE(10,*)
       WRITE(10,'(A)') '-----------------------------------------------'      
       WRITE(10,'(2(A,I3))') 
     *     'NPROT=',NPROT,'      NNEUT=',NNEUT
       WRITE(10,'(A)') '-----------------------------------------------'      
      ENDIF
C
      ITER=1
      CALL RMF(ITER)
      IF(IERROR.EQ.1) THEN
       KERROR=1
      ENDIF
C
 1000 CONTINUE
C
      RETURN
      END
C
C---   SCRINV   --------------------------------------------------------
C
      SUBROUTINE SCRINV(WFNA,ISO,JS,SPE0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'CONDIT.INC'
      INCLUDE 'CONSTA.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'WRKINV.INC'
      DIMENSION WFNA(NGRDM2)
C
C
      XISOFC=2.0-ISO
      DO 10 I=1,NGRID
      SQB(I)=1.0/(SPE0+XMSISO(I))
      IF(SQB(I).LE.0.D0) THEN
        IERROR=1
        RETURN
      ENDIF
   10 CONTINUE
C
      CALL DERIV0(SQB,D1B,ONE)
C
      POTMIN=1.0D10
      XK=XKAPPA(JS)
      XK2=(XK*XK+XK)
      SQB(1)=SQRT(SQB(1))
C
      DO 20 I=2,NGRID
      RX=RSTEP*(I-1)
      D1A=D1VPOT(I)
      POTL(I)=POTISO(I)+XK2*SQB(I)/(RX*RX)-D1B(I)*(XK/RX-D1A)
     *       +SQB(I)*(D2VPOT(I)+D1A*D1A-2.0*XK/RX*D1A)
      POTMIN=MIN(POTL(I),POTMIN)
      SQB(I)=SQRT(SQB(I))
   20 CONTINUE
C
      CALL DERIV0(SQB,D1B,ONE)
      CALL DERIV0(D1B,D2SQB,-ONE)
C
      DO 30 I=1,NGRID
      D2SQB(I)=D2SQB(I)/SQB(I)
   30 CONTINUE
C
      CALL LAPLAB(WFNA,WFNA,POTMIN-ENINV0,RSTEP,NGRID)
C
      CALL LOWCMP(WFNA,WFNA(NGRIDP),XKAPPA(JS),GPARIT(JS))
      DO 40 I=1,NGRID
      WFNA(I+NGRID)=SQB(I)*SQB(I)*WFNA(I+NGRID)
   40 CONTINUE
cc-suga
c      kki=kki+1
c      if(kki.gt.3000) then
c      do 50 i=1,ngrid
c      wfna(i+ngrid)=0.d0
c   50 continue
c      endif
cc-suga
C
      RETURN
      END
C
C---   SCRINVL   --------------------------------------------------------
C
      SUBROUTINE SCRINVL(WFNA,ISO,JS,SPE0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'CONDIT.INC'
      INCLUDE 'CONSTA.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'WRKINV.INC'
      INCLUDE 'YLAMDA.INC'
      DIMENSION WFNA(NGRDM2)
C
C
      XISOFC=2.0-ISO
      DO 10 I=1,NGRID
      SQB(I)=1.0/(SPE0+XMSISO(I))
      IF(SQB(I).LE.0.D0) THEN
        IERROR=1
        RETURN
      ENDIF
   10 CONTINUE
C
      CALL DERIV0(SQB,D1B,ONE)
C
      POTMIN=1.0D10
      XK=XKAPPA(JS)
      XK2=(XK*XK+XK)
      SQB(1)=SQRT(SQB(1))
C
      DO 20 I=2,NGRID
      RX=RSTEP*(I-1)
      D1A=D1VPOTL(I)
      POTL(I)=POTISO(I)+XK2*SQB(I)/(RX*RX)-D1B(I)*(XK/RX-D1A)
     *       +SQB(I)*(D2VPOTL(I)+D1A*D1A-2.0*XK/RX*D1A)
      POTMIN=MIN(POTL(I),POTMIN)
      SQB(I)=SQRT(SQB(I))
   20 CONTINUE
C
      CALL DERIV0(SQB,D1B,ONE)
      CALL DERIV0(D1B,D2SQB,-ONE)
C
      DO 30 I=1,NGRID
      D2SQB(I)=D2SQB(I)/SQB(I)
   30 CONTINUE
C
      CALL LAPLAB(WFNA,WFNA,POTMIN-ENINV0,RSTEP,NGRID)
C
      CALL LOWCMP(WFNA,WFNA(NGRIDP),XKAPPA(JS),GPARIT(JS))
      DO 40 I=1,NGRID
      WFNA(I+NGRID)=SQB(I)*SQB(I)*WFNA(I+NGRID)
   40 CONTINUE
cc-suga
c      kki=kki+1
c      if(kki.gt.3000) then
c      do 50 i=1,ngrid
c      wfna(i+ngrid)=0.d0
c   50 continue
c      endif
cc-suga
C
      RETURN
      END
C
C---   WFSTEP   --------------------------------------------------------
C
      SUBROUTINE WFSTEP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'CONDIT.INC'
      INCLUDE 'CONSTA.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'WORKSP.INC'
      INCLUDE 'WRKINV.INC'
      INCLUDE 'YMSTAR.INC'
C      
C
      CALL DERIV0(D1VPOT,D2VPOT,-ONE)
C
      DO 10 ISO=1,2
      RPFAC=2.0*ISO-3.0
      XISOFC=2.0-ISO
      DO 50 I=1,NGRID
      XMSISO(I)=BARMSE(ISO)+SPOTM(I)-VPOT(I)
     *         +RPFAC*(DPOT(I)-RPOT(I))
     *         -XISOFC*COULPT(I)
      POTISO(I)=SPOTM(I)+VPOT(I)+RPFAC*(DPOT(I)+RPOT(I))
     *         +XISOFC*COULPT(I)
   50 CONTINUE
      DO 10 JS=1,JSMAX(ISO)
      DO 20 NSHL=1,NSMAX(JS,ISO)
      NPL=NPLACE(NSHL,JS,ISO)
      NPLM=NPL-1
C
      CALL HPSI(WFN(NPL),AUX1,XKAPPA(JS),GPARIT(JS),ISO)
      SPE=WGTINT(1)*WFN(NPL)*AUX1(1)
      DO 30 I=2,NGRID2
      SPE=SPE+WGTINT(I)*AUX1(I)*WFN(I+NPLM)
   30 CONTINUE
      SPENRG(NPL/NGRID2+1)=(SPE-BARMSE(ISO))*HBARC
C
      CALL SCRINV(WFN(NPL),ISO,JS,SPE)
C
  20  CONTINUE
C
      CALL ORTHOG(NPLACE(1,JS,ISO)-1,NSMAX(JS,ISO))
C
  10  CONTINUE
C
      RETURN
      END
C
CSSS   WFSTEPL   --------------------------------------------------------
C
      SUBROUTINE WFSTEPL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'CONDIT.INC'
      INCLUDE 'CONSTA.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'WORKSP.INC'
      INCLUDE 'WRKINV.INC'
      INCLUDE 'YMSTAR.INC'
      INCLUDE 'YLAMDA.INC'
      INCLUDE 'LAMDA.INC'
      INCLUDE 'YSSMES.INC'
C
      CALL DERIV0(D1VPOTL,D2VPOTL,-ONE)
C
      DO 10 I=1,NGRID
      XMSISO(I)=YMSL(1)+SPOTML(I)-VPOTL(I)-COULPT(I)
     *                 +SSPOTL(I)-VSPOTL(I)
      POTISO(I)=SPOTML(I)+VPOTL(I)
     *                 +SSPOTL(I)+VSPOTL(I)+COULPT(I)
   10 CONTINUE
C
      JS=JSL
      ISO=2
      NPL=1
      NPLM=NPL-1
C
      CALL HPSIL(WFNL(NPL),AUX1,XKAPPA(JS),GPARIT(JS),ISO)
      SPE=WGTINT(1)*WFNL(NPL)*AUX1(1)
      DO 20 I=2,NGRID2
      SPE=SPE+WGTINT(I)*AUX1(I)*WFNL(I+NPLM)
   20 CONTINUE
      SPENRGL=(SPE-YMSL(1))*HBARC
C
      CALL SCRINVL(WFNL(NPL),ISO,JS,SPE)
C
      IPL1=NPLM
C
      OVLP=0.0
      DO 30 I=1,NGRID2
   30 OVLP=WFNL(I+IPL1)*WFNL(I+IPL1)*WGTINT(I)+OVLP
      OVLP=1./SQRT(OVLP)
      DO 60 I=1,NGRID2
   60 WFNL(I+IPL1)=WFNL(I+IPL1)*OVLP
C
      RETURN
      END
C
C---   ZPECOR  ---------------------------------------------------------
C
      SUBROUTINE ZPECOR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'CONSTA.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'RESULT.INC'
      INCLUDE 'SHELLS.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'WORKSP.INC'
      INCLUDE 'YLAMDA.INC'
      INCLUDE 'LAMDA.INC'
C
      DIMENSION WFAUXL(NGRDMX)
C
       DIAG   = 0.0
       OFFDIG = 0.0
       DO 10 ISO=1,2
C
       DO 20 JSA=1,JSMAX(ISO)
       IF(NSMAX(JSA,ISO).LE.0) GOTO 20
       XKA    = XKAPPA(JSA)
       XJA    = XJ(JSA)
       XLA    = XL(JSA)
       LA     = NINT(XLA)
       XLA2   = XLA*XLA+XLA
       YLA    = XJA+XJA-XLA
       LAL    = NINT(YLA)
       YLA2   = YLA+YLA*YLA
C
       DO 30 NSHLA=1,NSMAX(JSA,ISO)
       NPLA   = NPLACE(NSHLA,JSA,ISO)
       NPLMA  = NPLA-1
       NPLOLM = NGRID2*NSHLA-NGRID2
       NPLOL  = NPLOLM + 1
       CALL DERIV1(WFN(NPLA),WFAUX(NPLOL),GPARIT(JSA))
       SUM    = 0.0
       DO 35 I=2,NGRID
       J      = I+NGRID
       SUM    = WGTINT(I)*( WFAUX(I+NPLOLM)*WFAUX(I+NPLOLM)
     *                     +WFAUX(J+NPLOLM)*WFAUX(J+NPLOLM)
     *                     +( WFN(J+NPLMA)*WFN(J+NPLMA)*YLA2
     *                       +WFN(I+NPLMA)*WFN(I+NPLMA)*XLA2)
     *                        *ZENTRF(I))     + SUM
   35  CONTINUE
       DIAG   = WEIGHT(NPLMA/NGRID2+1)*XMULT(NPLMA/NGRID2+1)*(SUM+SUM)
     *          + DIAG
   30  CONTINUE
C
       DO 40 JSB=1,JSA
       IF(NSMAX(JSB,ISO).LE.0) GOTO 40
       XKB    = XKAPPA(JSB)
       XJB    = XJ(JSB)
       XLB    = XL(JSB)
       LB     = NINT(XLB)
       YLB    = XJB+XJB-XLB
       LBL    = NINT(YLB)
C
       IF(MOD(LA+LB,2).NE.1) GOTO 40
       IF(IABS(LA-LB).GT.1 .AND. IABS(LAL-LBL).GT.1) GOTO 40
       IF(ABS(XJA-XJB).GT.1.00) GOTO 40
C
       IF(IABS(LA-LB).EQ.1) THEN
        IF(XJA-XLA.EQ.XJB-XLB) THEN
         SIXFCU = SQRT((XLA+XJB+2.5)*(XLA+XJB-0.5)
     *                 /((XLA+XLA+1.00)*(XLB+XLB+1.00)))
        ELSE
         SIXFCU = SQRT((XJB-XLA+1.5)*(XLA-XJB+1.5)
     *                 /((XLA+XLA+1.00)*(XLB+XLB+1.00)))
     *            *2.0*(XJB-XLB)
        ENDIF
        IF(XLA.GT.XLB) THEN
         FACDRU = SQRT(XLA)
         FACCTU = -XLA
        ELSE
         FACDRU = -SQRT(XLB)
         FACCTU = XLB
        ENDIF
       ELSE
        SIXFCU = 0.0
        FACDRU = 0.0
        FACCTU = 0.0
       ENDIF
C
       IF(IABS(LAL-LBL).EQ.1) THEN
        IF(XJA-YLA.EQ.XJB-YLB) THEN
         SIXFCL = SQRT((YLA+XJB+2.5)*(YLA+XJB-0.5)
     *                 /((YLA+YLA+1.00)*(YLB+YLB+1.00)))
        ELSE
         SIXFCL = SQRT((XJB-YLA+1.5)*(YLA-XJB+1.5)
     *                 /((YLA+YLA+1.00)*(YLB+YLB+1.00)))
     *            *2.0*(XJB-YLB)
        ENDIF
        IF(YLA.GT.YLB) THEN
         FACDRL = SQRT(YLA)
         FACCTL = -YLA
        ELSE
         FACDRL = -SQRT(YLB)
         FACCTL = YLB
        ENDIF
       ELSE
        SIXFCL = 0.0
        FACDRL = 0.0
        FACCTL = 0.0
       ENDIF
C
       DO 50 NSHLA=1,NSMAX(JSA,ISO)
       NPLA   = NPLACE(NSHLA,JSA,ISO)
       NPLMA  = NPLA-1
       V2A    = WEIGHT(NPLA/NGRID2+1)
       UVA    = SQRT(V2A-V2A*V2A)
       NPLOLM = NGRID2*NSHLA-NGRID2
       DO 50 NSHLB=1,NSMAX(JSB,ISO)
       NPLB   = NPLACE(NSHLB,JSB,ISO)
       NPLMB  = NPLB-1
       V2B    = WEIGHT(NPLB/NGRID2+1)
       UVB    = SQRT(V2B-V2B*V2B)
       RADMTU = 0.0
       RADMTL = 0.0
       DO 60 I=2,NGRID
       RA     = I*RSTEP-RSTEP
       RADMTU = WGTINT(I)*WFN(I+NPLMB)
     *                   *(-WFAUX(I+NPLOLM)+FACCTU*WFN(I+NPLMA)/RA)
     *          + RADMTU
       J      = I+NGRID
       RADMTL = WGTINT(I)*WFN(J+NPLMB)
     *                   *(-WFAUX(J+NPLOLM)+FACCTL*WFN(J+NPLMA)/RA)
     *          + RADMTL
   60  CONTINUE
       OFFDIG = (RADMTU*FACDRU*SIXFCU+RADMTL*FACDRL*SIXFCL)**2
     *           *(V2A*V2B+UVA*UVB)  + OFFDIG
C
   50  CONTINUE
   40  CONTINUE
   20  CONTINUE
   10  CONTINUE
C
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C     ADD LAMDA
C
      IF (NLAMDA.EQ.0) GOTO 500
C
       JSA=JSL
C
       XKA    = XKAPPA(JSA)
       XJA    = XJ(JSA)
       XLA    = XL(JSA)
       LA     = NINT(XLA)
       XLA2   = XLA*XLA+XLA
       YLA    = XJA+XJA-XLA
       LAL    = NINT(YLA)
       YLA2   = YLA+YLA*YLA
C
       NPLA   = 1
       NPLMA  = NPLA-1
       NPLOLM = 0
       NPLOL  = NPLOLM + 1
C 
       CALL DERIV1(WFNL(NPLA),WFAUXL(NPLOL),GPARIT(JSA))
       SUM    = 0.0
       DO 535 I=2,NGRID
       J      = I+NGRID
       SUM    = WGTINT(I)*( WFAUXL(I+NPLOLM)*WFAUXL(I+NPLOLM)
     *                     +WFAUXL(J+NPLOLM)*WFAUXL(J+NPLOLM)
     *                     +( WFNL(J+NPLMA)*WFNL(J+NPLMA)*YLA2
     *                       +WFNL(I+NPLMA)*WFNL(I+NPLMA)*XLA2)
     *                        *ZENTRF(I))     + SUM
  535  CONTINUE
       DIAG   = NLAMDA*SUM + DIAG
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
 500   CONTINUE
C
       PWIDTH = (DIAG-OFFDIG-OFFDIG)
       ZPE    = HBARC*0.5*PWIDTH/(APROT*BARMSE(1)+ANEUT*BARMSE(2)
     *                           +NLAMDA*YMSL(1))
C
      RETURN
      END
C      
C---   PRISPE   --------------------------------------------------------
C
      SUBROUTINE PRISPE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C      
      INCLUDE 'SHELLS.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'WORKSP.INC'
      CHARACTER*3 TX
      CHARACTER STATE(0:8)*2
      INTEGER NSAV(NSHLMX),LSAV(NSHLMX),JSAV(NSHLMX)
C       
      DATA TX/'/2='/
      DATA STATE/ ' s',' p',' d',' f',' g',' h',' i',' j',' k'/
      EQUIVALENCE (AUX1(1),NSAV(1)),(AUX1(NSHLMX+1),LSAV(1)),
     *            (AUX2(1),JSAV(1))
C
C     Print single particle energies with shell nummber
C
      WRITE(10,'(2(/),A)') 'S.P. ENERGIES (IN MEV) 
     * AND WEIGHT WITH SHELL NUMBER '
      NPL=0
      DO 10 ISO=1,2
      NPLA=NPL+1
      IF(ISO.EQ.1) THEN
        WRITE(10,'(A)') '  FOR PROTONS:'
        WRITE(10,'(A)') '  N   L   J        ENERGY       WEIGHT'
      ENDIF
      IF(ISO.EQ.2) THEN
        WRITE(10,'(A)') '  FOR NEUTRONS:'
        WRITE(10,'(A)') '  N   L   J        ENERGY       WEIGHT'
      ENDIF
c
      DO 20 JS=1,JSMAX(ISO)
      DO 20 NSHL=1,NSMAX(JS,ISO)
      NPL=NPL + 1
      JSAV(NPL)=XJ(JS)*2
      LSAV(NPL)=INT(XL(JS)*1)
      NSAV(NPL)=NSHL
   20 CONTINUE
      NPLE=NPL
C
      ABSMIN=-1.E+10
      DO 30 JJ=NPLA,NPLE
      AMIN=1.E+10
      DO 40 KK=NPLA,NPLE
      IF(SPENRG(KK).GT.ABSMIN  .AND.  SPENRG(KK).LT.AMIN) THEN
       AMIN=SPENRG(KK)
       KMIN=KK
      ENDIF
   40 CONTINUE
      ABSMIN=AMIN
      N=KMIN
      WRITE(10,1529)
     *       NSAV(N),STATE(LSAV(N)),JSAV(N),TX,SPENRG(N),WEIGHT(N)
 1529 FORMAT(1X,I2,2X,A2,2X,I2,A3,2F12.5)
   30 CONTINUE 
C
   10 CONTINUE
C
      RETURN
	END
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
C     EDITED BY K. SUMIYOSHI
C     ORIGINAL TAKEN FROM SUBSKY8.FOR
C
C     FIVE POINT DIFFERENTATION FORMULA (ABRAMOWITZ, PAG 914)
C
      SUBROUTINE DERIV(F,DF,H)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(5,5)
      DIMENSION F(5),DF(5)
C
      DATA A(1,1),A(1,2),A(1,3),A(1,4),A(1,5)
     +                  /-50.,96.,-72.,32.,-6./
      DATA A(2,1),A(2,2),A(2,3),A(2,4),A(2,5)
     +                  /-6.,-20.,36.,-12.,2./
      DATA A(3,1),A(3,2),A(3,3),A(3,4),A(3,5)
     +                  /2.,-16.,0.,16.,-2./
      DATA A(4,1),A(4,2),A(4,3),A(4,4),A(4,5)
     +                     /-2.,12.,-36.,20.,6./
      DATA A(5,1),A(5,2),A(5,3),A(5,4),A(5,5)
     +                     /6.,-32.,72.,-96.,50./
      DATA FACT/24./
C
      DO 2000 J=1,5
C
      SUM=0.0D0
      DO 1000 I=1,5
      SUM=SUM+A(J,I)*F(I)
 1000 CONTINUE
      DF(J)=SUM/(H*FACT)
C
 2000 CONTINUE
C
      RETURN
      END
CSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
      SUBROUTINE DCOUPL
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PAR000.INC'
C
      INCLUDE 'CONSTA.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'FORCES.INC'
      INCLUDE 'WAVEFN.INC'
      INCLUDE 'YMSTAR.INC'
      INCLUDE 'YLAMDA.INC'
      INCLUDE 'LAMDA.INC'
C

C
      SIG0=VDRH*SMESON(2)-ESXH*SMESON(3)
      CALL SPLINE(SIG0,XMS,GS,XMSL,GSL)
C
      SMESCPD(1)=GS
      SPOTM(1)=XMS-BARMSE(1)
      SPOTD(1)=GS*SIG0
C
      SMESCPDL(1)=GSL
      SPOTML(1)=XMSL-YMSL(1)
      SPOTDL(1)=GSL*SIG0
C
      DO 10 I=2,NGRID
      RX=RSTEP*(I-1)
      SIG0=SMESON(I)/RX
      CALL SPLINE(SIG0,XMS,GS,XMSL,GSL)
C
      SMESCPD(I)=GS
      SPOTM(I)=XMS-BARMSE(1)
      SPOTD(I)=GS*SIG0
C
      SMESCPDL(I)=GSL
      SPOTML(I)=XMSL-YMSL(1)
      SPOTDL(I)=GSL*SIG0
   10 CONTINUE
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C---   SPLINE   --------------------------------------------------------
C
      SUBROUTINE SPLINE(X,YM,YG,YML,YGL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'PAR000.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'YMSTAR.INC'
      INCLUDE 'YLAMDA.INC'
      INCLUDE 'LAMDA.INC'
C
      I1=1
      I2=2
C
      DO 1000 I=2,NNY
      IF (X .LT. YSIG(I)) THEN
      GOTO 1000
      ELSE
      I1=I-1
      I2=I
      GOTO 1100
      ENDIF
 1000 CONTINUE
      I1=NNY-1
      I2=NNY
C
 1100 CONTINUE
C
      DX=YSIG(I2)-YSIG(I1)
C            
      DM=(YMS(I2)-YMS(I1))/DX
      YM=DM*(X-YSIG(I1))+YMS(I1)
C
      DG=(YDMS(I2)-YDMS(I1))/DX
      YG=DG*(X-YSIG(I1))+YDMS(I1)
C
      DML=(YMSL(I2)-YMSL(I1))/DX
      YML=DML*(X-YSIG(I1))+YMSL(I1)
C
      DGL=(YDMSL(I2)-YDMSL(I1))/DX
      YGL=DGL*(X-YSIG(I1))+YDMSL(I1)
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC