C
C        MAIN PROGRAM FOR THE TEST RUN
C
      CALL MEDOC
      WRITE(6,1)
    1 FORMAT(1H1)
      STOP
      END
      SUBROUTINE MEDOC
C
C        CALCULATION OF THE NORM AND COLLISION MATRIX ELEMENTS FOR OEDM.
C
      DOUBLE PRECISION AO,CRAD,DELTA,DFLOAT,E,EP,FNOR,FNORP,GAMRAD,
     1OMEZUT,PE,PEP,QU,R,RP,RR,SIGMA,SIGMAP,XF,XFP,XG,XGP
      INTEGER INDEX1,INDEX2
      CHARACTER*260 FILENAME
      DIMENSION RR(300)
      DIMENSION N1(10),L1(10),M1(10)
C
C        COMMON/LAG/ CONTAINS THE ABCISSA AND WEIGHTS FOR GAUSS LAGUERRE
C     QUADRATURES.
C
      COMMON/LAG/XI(50),OMEGA(50),NXI
C
C        COMMON BRA CONTAINS THE COEFFICIENTS OF THE SEMI-ANALYTICAL
C     EXPANSION OF THE BRA VECTOR IN THE MATRIX ELEMENT.
C
      COMMON/BRA/XFP(40),XGP(40),NFP,NGP
C
C        COMMON KET CONTAINS THE COEFFICIENTS OF THE SEMI-ANALYTICAL
C     EXPANSION OF THE KET VECTOR IN THE MATRIX ELEMENT.
C
      COMMON/KET/XF(40),XG(40),NF,NG
C
C        THE FOLLOWING COMMONS ARE USED FOR THE TRANSMISSION OF
C     VARIABLES TO THE SUBROUTINES
C
      COMMON//R,PE,PEP,SIGMA,SIGMAP,DELTA,L,ME
      COMMON/ORIGE/OMEZUT
      OPEN(77,FILE='MEDOC_INPUT')
      OPEN(99,FILE='IOEDM')
      OPEN(777,FILE='r.txt')
      OPEN(105,FILE='rot.txt')
      OPEN(106,FILE='rad.txt')

    1 FORMAT (4E18.8)
    3 FORMAT (D12.6,4I3/(D20.14))
    4 FORMAT (3I3)
    5 FORMAT (D20.10)
    6 FORMAT (44H MATRIX ELEMENTS ARE ZERO BETWEEN THE STATES,3(I3,2X),1ACXY0049
     1H/,3(I3,2X))
    9 FORMAT (D12.8)
   91 FORMAT (1X,6HZ2/Z1=,1PD20.10,5X,2HR=,D15.8,4HU.A.,3X,1H*,3I2,4H/DRACXY0052
     1/,3I2,2H*=,1PD15.8,3X,5HGAMA=,D15.8)
   92 FORMAT (1X,6HZ2/Z1=,1PD20.10,5X,2HR=,D15.8,4HU.A.,3X,1H*,3I2,4H/LYACXY0054
     1/,3I2,2H*=,1PD15.8,3X,5HGAMA=,D15.8)
   93 FORMAT (F6.3,6I2,4D15.8)
   94 FORMAT (1H1,37HDISTANCE FROM CENTER 1 TO THE ORIGIN=,D20.10,////)
   34 FORMAT (E9.4,E15.7)

C
C        INPUT FROM UNIT IN:
C
C        AO= DISTANCE OF THE ORIGIN OF THE ELECTRONIC COORDINATE SYSTEM
C     TO THE CHARGE Z1 IN UNIT OF R. (SEE EQ. 7 OF L.W.UP).
C        NXI=NUMBER OF PIVOTS IN THE GAUSS-LAGUERRE INTEGRATION.
C        XI,OMEGA= ABCISSAS AND WEIGHTS FOR THE GAUSS-LAGUERRE INTEGRA-
C     TION ( NXI SINGLE PRECISION VALUES).
C        JPERF= IF JPERF IS NON-ZERO, RESULTS ARE WRITTEN ON UNIT 7.
C        NR= TWO POSSIBILITIES:
C
C     *****NR NON ZERO:
C        THEN NR IS THE NUMBER OF VALUES OF THE INTERNUCLEAR DISTANCE
C     AND THE FOLLOWING QUANTITIES ARE ALSO READ ON UNIT IN:
C        RR= ARRAY OF VALUES OF THE INTERNUCLEAR DISTANCE ( NR VALUES).
C        QU= RATIO OF NUCLEAR CHARGES.
C        KN= NUMBER OF STATES ( KN.GE.2 ).
C        N1,L1,M1= ARRAYS OF QUANTUM NUMBERS OF THE MOLECULAR STATES IN
C     THE UNITED ATOM LIMIT ( KN VALUES).
C
C     *****NR IS ZERO:
C        THEN THE ONLY FURTHER QUANTITY THAT SHOULD BE GIVEN ON UNIT IN
C     IS:
C        KN= NUMBER OF STATES (DEFAULT VALUE: 2).
C     AND THE PROGRAM READS THE OTHER DATA ON UNIT IOEDM FOR EACH STATE
C     (KN SETS OF DATA).
C        QU= RATION OF NUCLEAR CHARGES
C        N1,L1,M1= ARRAYS OF QUANTUM NUMBERS OF THE MOLECULAR STATES IN
C     THE UNITED ATOM LIMIT.
C        NR= NUMBER OF VALUES OF R.
C        RR= ARRAY OF VALUES OF R ( INTERNUCLEAR DISTANCE).
C
      READ(77,*) AO
      READ(77,*) NXI
      READ(77,*) (XI(N),OMEGA(N),N=1,NXI)
      READ(77,*) JPERF
      READ(777,*) NR
      READ(777,*) (RR(N),N=1,NR)
      READ(77,*) QU
      READ(77,*) KN
      READ(77,*) (N1(K),L1(K),M1(K),K=1,KN)
      READ(77,*) FILENAME
      GO TO 13
   10 READ(77,*) KN
      IF(KN.EQ.0) KN=2
      DO 12 K=1,KN
      READ(99,*) QU,N1(K),L1(K),M1(K),NR,(RR(N),N=1,NR)
      write(*,*) 'called?'
   12 CONTINUE
   13 CONTINUE
C
      OMEZUT=1.D0-2.D0*AO
      WRITE(6,94) AO
      IF(JPERF.NE.0) WRITE(7,94) AO
      DO 200 K=1,KN-1
      DO 300 KK=K+1,KN
      IF(IABS(M1(K)-M1(KK))-1) 15,18,11
   18 CONTINUE
      IF(M1(K).GT.M1(KK)) GO TO 16
   15 CONTINUE
      NN=N1(K)
      L=L1(K)
      ME=M1(K)
      INDEX1 = 100*NN+10*L+ME
      NNP=N1(KK)
      LP=L1(KK)
      MEP=M1(KK)
      INDEX2 = 100*NNP+10*LP+MEP
      GO TO 17
   16 CONTINUE
      NN=N1(KK)
      L=L1(KK)
      ME=M1(KK)
      INDEX1 = 100*NN+10*L+ME
      NNP=N1(K)
      LP=L1(K)
      MEP=M1(K)
      INDEX2 = 100*NNP+10*LP+MEP
   17 CONTINUE

C
C        LOOP ON INTERNUCLEAR DISTANCES
C
      OPEN(111,FILE=FILENAME)

      DO 100 JN=1,NR
      R=RR(JN)
      RP=R
c     KET, 210
      CALL ENTREE(NN,L,ME,R,E,PE,NG,XG,NF,XF,QU,INDEX1)
c     BRA, 211
      CALL ENTREE(NNP,LP,MEP,R,EP,PEP,NGP,XGP,NFP,XFP,QU,INDEX2)
      SIGMA=R*(1.D0+QU)/(2.D0*PE)-1.D0-DFLOAT(ME)
      SIGMAP=RP*(1.D0+QU)/(2.D0*PEP)-1.D0-DFLOAT(MEP)
      IF(MEP.NE.ME) GO TO 14
C
C        RADIAL COUPLING CASE
C

      CALL NORDIF(ME,PE,SIGMA,NG,XG,NF,XF,FNOR,R)

      CALL NORDIF(MEP,PEP,SIGMAP,NGP,XGP,NFP,XFP,FNORP,RP)
      WRITE(111,*) FNOR*FNORP
      DELTA=E-EP
      CALL RADAS(CRAD,QU,GAMRAD)
      GAMRAD=GAMRAD*FNOR*FNORP
      CRAD=CRAD*FNOR*FNORP
      WRITE(6,91) QU,R,NNP,LP,MEP,NN,L,ME,CRAD,GAMRAD
      WRITE(106,34) R, CRAD
      IF(JPERF.NE.0) WRITE(7,93) QU,NNP,LP,MEP,NN,L,ME,R,CRAD,E,EP
      GO TO 100
C
C        ROTATIONAL COUPLING CASE.
C
   14 CONTINUE

      CALL NORDIF(ME,PE,SIGMA,NG,XG,NF,XF,FNOR,R)
      CALL NORDIF(MEP,PEP,SIGMAP,NGP,XGP,NFP,XFP,FNORP,R)
      WRITE(111,*) FNOR*FNORP
      CALL ROTDIS(CRAD,GAMRAD)

      GAMRAD=GAMRAD*FNOR*FNORP
      CRAD=CRAD*FNOR*FNORP
      WRITE(6,92) QU,R,NNP,LP,MEP,NN,L,ME,CRAD,GAMRAD
      WRITE(105,34) R,CRAD
      IF(JPERF.NE.0) WRITE(7,93) QU,NNP,LP,MEP,NN,L,ME,R,CRAD,E,EP
C
  100 CONTINUE
C
C        OUTPUT ON UNIT 6:
C        QU= RATIO OF NUCLEAR CHARGES.
C        NNP,LP,MEP= QUANTUM NUMBERS DEFINING THE "BRA" STATE.
C        NN,L,ME= QUANTUM NUMBERS DEFINING THE "KET" STATE.
C        R= INTERNUCLEAR DISTANCE.
C        CRAD= COUPLING MATRIX ELEMENT ( EQ. 7 OF THE L.W.UP).         0
C        GAMRAD= MATRIX ELEMENT OF THE OPERATOR D2 OR L2 (EQ. 9 OR 11 OF
C     THE L.W.UP).
C
C
C        OUTPUT ON UNIT 7 ( IF JPERF IS NON-ZERO).
C        QU,NNP,LP,MEP,NN,L,ME,R,CRAD= AS ON UNIT 6.
C        E= ENERGY FOR THE "KET" STATE.
C        EP= ENERGY FOR THE "BRA" STATE.
C
      GO TO 300
   11 WRITE(6,6) N1(K),L1(K),M1(K),N1(KK),L1(KK),M1(KK)
  300 CONTINUE
  200 CONTINUE

      RETURN
      END
      SUBROUTINE ENTREE(N,L,ME,R,E,PE,NG,XG,NF,XF,QU,INDEX)
C
C        SUBROUTINE USED TO READ THE CHARACTERISTICS OF ONE STATE
C     THROUGH SUBROUTINE LECISA.
C
      DOUBLE PRECISION A,AOUT,E,EE,PE,PEPE,PU,QU,R,RR,XF,XG,XXF,XXG
      INTEGER INDEX
      DIMENSION XG(1),XF(1)
      DIMENSION AOUT(88)
      EQUIVALENCE (PU,AOUT)
      COMMON/TRAP/PU,RR,MME,NN,LL,IVER,EE,PEPE,A,NNG,NNF,XXG(40),XXF(40)
      DATA IVER/2/
      PU=QU
      NN=N
      LL=L
      MME=ME
      RR=R
      CALL LECISA(INDEX)
      E=EE
      PE=PEPE
      NG=NNG
      NF=NNF
      DO 10 I=1,NG
   10 XG(I)=XXG(I)
      DO 20 I=1,NF
   20 XF(I)=XXF(I)
      RETURN
      END
      SUBROUTINE LECISA(INDEX)
C
C
C        USER SUPPLIED SUBROUTINE USED TO PICK UP THE CHARACTERISTICS OF
C     THE STATE CONSIDERED AS CALCULATED IN THE PROGRAM GRAVE. SHOULD BE
C      WRITTEN IN MACHINE LANGUAGE.
C
C
C        THE VERSION GIVEN HERE IS ONLY SUITABLE FOR THE TEST RUN AND
C     SHOULD NOT BE CONSIDERED AS AN INDICATION ON THE ROLE OF THE
C     SUBROUTINE. SEE LONG WRITE UP FOR DETAILS ON THE CONCEPTION OF
C     THE SUBROUTINE.
C
      DOUBLE PRECISION AA,B,C,E,PE,A,      XG,XF
      INTEGER INDEX
      DIMENSION AA(1)
      COMMON/TRAP/B,C,M1,M2,M3,M4,E,PE,A,NG,NF,XG(40),XF(40)
      READ(INDEX,*) E,PE,A,NG,NF
      READ(INDEX,*) (XG(I),I=1,NG)
      READ(INDEX,*) (XF(I),I=1,NF)
    1 FORMAT(3D20.14,2I3,/,(4D18.12))
      RETURN
      END

      SUBROUTINE NORDIF(ME,PE,SIGMA,NG,XG,NF,XF,FNOR,R)
C
C        CALCULATES THE NORM OF THE FUNCTION.
C
C
C     IMPORTANT REMARK: IN THE PROGRAM GRAVE, THE LAST TERM IN THE SEMI
C     ANALYTICAL EXPANSION IS SET EQUAL TO ONE. THE SUBROUTINE NORDIF
C     FIRST TRANSFORMS THE COEFFICIENTS SO THAT THE FIRST NON-ZERO TERM
C     IS NOW EQUAL TO ONE. THE NORM IS CALCULATED FOR THIS NEW FORM
C     OF EXPANSION.
      DOUBLE PRECISION ALAM,ALAMXI,ALAM1,ALAM2,ALPHA,BETA,DEXP,DFLOAT,
     1DLOG,DSQRT,FNOR,FONC,PE,PP,R,SAM,SG,SIGMA,SIM,SPG,SUM,TROC,TRUC,
     2XF,XF2,XG,XN,Y
      INTEGER NG,NF
      DIMENSION XG(NG),XF(NF)

      COMMON/LAG/XI(50),OMEGA(50),NXI
      XN=XG(1)
      DO 16 K=1,NG
   16 XG(K)=XG(K)/XN
      XN=XF(1)
      IF(XN.EQ.0.D0) XN=XF(2)
      DO 15 K=1,NF
   15 XF(K)=XF(K)/XN
C
C        INTEGRALES SUR MU
C
      IALPHA=1
      ME2=2*ME
      DO 10 K=1,ME2
   10 IALPHA=IALPHA*K
      ALPHA=DFLOAT(IALPHA)
      SAM=0.D0
      SUM=0.D0
      DO 11 IT=1,NF
      IS=IT-1
      XF2=XF(IT)*XF(IT)
      IM=ME2+IS
      IME=IM+IS+1
      SUM=SUM+ALPHA*2.D0*XF2/DFLOAT(IME)
      TROC=DFLOAT((IT)*(IM+1))/DFLOAT(IME+2)
      TROC=TROC+DFLOAT(IS*IM)/DFLOAT(IME-2)
      TROC=TROC/DFLOAT(IME*IME)
      SAM=SAM+ALPHA*2.D0*XF2*TROC
      ALPHA=ALPHA*DFLOAT(IM+1)/DFLOAT(IT)
   11 CONTINUE
      BETA=DFLOAT(IALPHA*(ME2+1)*(ME2+2))
      SIM=0.D0
      IF(NF.LE.2) GO TO 21
      NFP=NF-2
      DO 20 IT=1,NFP
      IS=IT-1
      IME=ME2+IS+IS+1
      SIM=SIM+4.D0*XF(IT)*XF(IT+2)*BETA/DFLOAT(IME*(IME+2)*(IME+4))
      BETA=BETA*DFLOAT(ME2+IS+3)/DFLOAT(IT)
   20 CONTINUE
   21 SAM=SAM+SIM
C
C        INTEGRALES SUR LAMBDA
C
      SG=0.D0
      SPG=0.D0
      PP=2.D0*PE
      DO 30 N=1,NXI
c new variable as lambda
      ALAM=XI(N)/PP+1.D0
      ALAM1=ALAM-1.D0
      ALAM2=ALAM+1.D0
      ALAMXI=ALAM1/ALAM2
      FONC=DEXP(-PE*ALAM+SIGMA*DLOG(ALAM2))
c (FONC = (lam+1)^sig*exp(-p*lam))
      FONC=FONC*FONC
      TRUC=ALAM1*ALAM2
c TRUC = (lam+1)*(lam-1)
   35 TRUC=TRUC**ME
c ^m/2*^m/2 = ^m
      FONC=FONC*TRUC
c rwang Calculate terms before sum
c      Y=XG(NG)
c      DO 36 JN=2,NG
c      NN=NG-JN+1
c   36 Y=Y*ALAMXI+XG(NN)
      Y = 0.d0
      DO 36 JN =1,NG
   36 Y = Y+XG(JN)*ALAMXI**(JN-1)

      FONC=FONC*Y*Y*OMEGA(N)
      SG=SG+FONC
      SPG=SPG+FONC*ALAM*ALAM
   30 CONTINUE

      SG=SG/PP
      SPG=SPG/PP
C

      FNOR=R*R*R*(SPG*SUM-SG*SAM)/8.D0


      FNOR=DSQRT(1.D0/FNOR)
      RETURN
      END




      SUBROUTINE RADAS(CRAD,QU,GAMRAD)
C
C        CALCULATES THE RADIAL MATRIX ELEMENTS.
C
      DOUBLE PRECISION A,AI,ALAM,ALAMXI,ALAM1,ALAM2,ALPHA,B,BI,C,CI,
     1COEFF,CRAD,DELTA,DEXP,DFLOAT,DLOG,DY,FISC,FONC,GAMRAD,OMEZUT,
     2PE,PEP,QU,R,SAM,SG,SIGMA,SIGMAP,SIGPE,SIM,SJM,SPG,SPH,SPI,SUM,
     3TRIC,TRUC,XF,XFP,XG,XGP,Y,YP,SAG3,SAG2,SPK,SPK2

      A(I,I2)=(-DFLOAT((I+1)*(ME+I)*(ME2+I+1))/DFLOAT(ME2+I2+3)
     1+DFLOAT((ME2+I)*(ME+I+1)*I)/DFLOAT(ME2+I2-1))/DFLOAT(ME2+I2+1)
      B(I,I2)=DFLOAT((ME+I+3)*(ME2+I+2)*(ME2+I+1))/DFLOAT((ME2+I2+5)*
     1(ME2+I2+3))
      C(I,I2)=-DFLOAT((I-1)*I*(ME+I-2))/DFLOAT((ME2+I2-3)*(ME2+I2-1))
      COMMON/LAG/XI(50),OMEGA(50),NXI
      COMMON/BRA/XFP(40),XGP(40),NFP,NGP
      COMMON/KET/XF(40),XG(40),NF,NG
      COMMON//R,PE,PEP,SIGMA,SIGMAP,DELTA,L,ME
      COMMON/ORIGE/OMEZUT
      ME2=2*ME
C
C        INTEGRALS OVER MU
C
      IALPHA=1
      DO 11 K=1,ME2
   11 IALPHA=IALPHA*K
C
C        INTEGRAL E
C

      ALPHA=DFLOAT(IALPHA)
      NMAX=MIN0(NFP,NF)
      SUM=0.D0
      DO 10 IT=1,NMAX
      IS=IT-1
      IM=ME2+IS+1
      SUM=SUM+2.D0*ALPHA*XF(IT)*XFP(IT)/DFLOAT(IM+IS)
      ALPHA=ALPHA*DFLOAT(IM)/DFLOAT(IT)
   10 CONTINUE
C
C        INTEGRAL  F
C
      ALPHA=DFLOAT(IALPHA)
      NMAX=MIN0(NFP,NF-2)
      IF(NF.LT.3) XF(3)=0.D0
      IS=0
      IS2=0
      COEFF=2.D0/DFLOAT(ME2+1)
      SAM=XFP(1)*(A(IS,IS2)*XF(1)+B(IS,IS2)*XF(3))*COEFF*ALPHA
      ALPHA=ALPHA*DFLOAT(ME2+1)
      IF(NF.EQ.1) GO TO 28
      IF(NFP.EQ.1) GO TO 21
      IF(NF.LT.4) XF(4)=0.D0
      IS=1
      IS2=2
      COEFF=2.D0/DFLOAT(ME2+3)
      SAM=SAM+XFP(2)*(A(IS,IS2)*XF(2)+B(IS,IS2)*XF(4))*COEFF*ALPHA
      ALPHA=ALPHA*DFLOAT(ME2+2)/2.D0
      IF(NFP.EQ.2) GO TO 21
      IF(NF.EQ.2) GO TO 27
      IF(NF.EQ.3) GO TO 23
      IF(NF.EQ.4) GO TO 19
      DO 20 IT=3,NMAX
      IS=IT-1
      MIM=ME+IS
      IM=MIM+ME+1
      IME=IM+IS
      COEFF=2.D0/DFLOAT(IME)
      AI=(-DFLOAT(IT*MIM*IM)/DFLOAT(IME+2)+DFLOAT((IM-1)*(MIM+1)*IS)/
     1DFLOAT(IME-2))/DFLOAT(IME)
      BI=DFLOAT((MIM+3)*(IM+1)*IM)/DFLOAT((IME+4)*(IME+2))
      CI=-DFLOAT(IS*(IS-1)*(MIM-2))/DFLOAT((IME-4)*(IME-2))
      TRUC=XFP(IT)*(AI*XF(IT)+BI*XF(IT+2)+CI*XF(IT-2))*ALPHA*COEFF
      SAM=SAM+TRUC
      ALPHA=ALPHA*DFLOAT(IM)/DFLOAT(IT)
   20 CONTINUE
   19 IF(NMAX.EQ.NFP) GO TO 21
      IT=NMAX+1
      IF(IT.LE.2) GO TO 23
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0/DFLOAT(ME2+IS2+1)
      TRUC=XFP(IT)*(A(IS,IS2)*XF(IT)+C(IS,IS2)*XF(IT-2))*ALPHA*COEFF
      SAM=SAM+TRUC
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)
   23 IF(NFP.EQ.(NF-1)) GO TO 21
      IF(NFP-NF-1) 24,25,26
   24 IT=NF
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0/DFLOAT(ME2+IS2+1)
      TRUC=XFP(IT)*(A(IS,IS2)*XF(IT)+C(IS,IS2)*XF(IT-2))*ALPHA*COEFF
      SAM=SAM+TRUC
      GO TO 21
   28 IF(NFP.LT.3) GO TO 21
      IT=3
      IS=IT-1
      IS2=2*IS
      ALPHA=ALPHA*DFLOAT(ME2+2)/2.D0
      COEFF=2.D0/DFLOAT(ME2+IS2+1)
      TRUC=XFP(IT)*C(IS,IS2)*XF(IT-2)*ALPHA*COEFF
      SAM=SAM+TRUC
      GO TO 21
   27 IT=3
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0/DFLOAT(ME2+IS2+1)
      TRUC=XFP(IT)*C(IS,IS2)*XF(IT-2)*ALPHA*COEFF
      SAM=SAM+TRUC
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)
      IF(NFP.EQ.3) GO TO 21
      IT=4
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0/DFLOAT(ME2+IS2+1)
      TRUC=XFP(IT)*C(IS,IS2)*XF(IT-2)*ALPHA*COEFF
      SAM=SAM+TRUC
      GO TO 21
   25 IT=NF
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0/DFLOAT(ME2+IS2+1)
      TRUC=XFP(IT)*(A(IS,IS2)*XF(IT)+C(IS,IS2)*XF(IT-2))*ALPHA*COEFF
      SAM=SAM+TRUC
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)
      IT=NF+1
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0/DFLOAT(ME2+IS2+1)
      TRUC=XFP(IT)*C(IS,IS2)*XF(IT-2)*ALPHA*COEFF
      SAM=SAM+TRUC
      GO TO 21
   26 IT=NF
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0/DFLOAT(ME2+IS2+1)
      TRUC=XFP(IT)*(A(IS,IS2)*XF(IT)+C(IS,IS2)*XF(IT-2))*ALPHA*COEFF
      SAM=SAM+TRUC
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)
      IT=NF+1
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0/DFLOAT(ME2+IS2+1)
      TRUC=XFP(IT)*C(IS,IS2)*XF(IT-2)*ALPHA*COEFF
      SAM=SAM+TRUC
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)
      IT=NF+2
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0/DFLOAT(ME2+IS2+1)
      TRUC=XFP(IT)*C(IS,IS2)*XF(IT-2)*ALPHA*COEFF
      SAM=SAM+TRUC
   21 CONTINUE
C
C        INTEGRALS A AND B
C
      ALPHA=DFLOAT(IALPHA)
      SIM=2.D0*ALPHA*XFP(1)*XF(2)/DFLOAT(ME2+3)
      SJM=2.D0*ALPHA*XFP(1)*XF(2)*DFLOAT(ME+2)/DFLOAT(ME2+3)
      ALPHA=ALPHA*DFLOAT(ME2+1)
      NMAX=MIN0(NFP,NF-1)
      IF(NMAX.EQ.1) GO TO 300
      DO 30 IT=2,NMAX
      IS=IT-1
      IM=ME2+IS+1
      IME=IM+IS
      COEFF=2.D0*ALPHA/DFLOAT(IME)
      SIM=SIM+XFP(IT)*COEFF*(XF(IT-1)*DFLOAT(IS)/DFLOAT(IME-2)+
     1XF(IT+1)*DFLOAT(IM)/DFLOAT(IME+2))
      SJM=SJM+XFP(IT)*COEFF*(-XF(IT-1)*DFLOAT(IS*(ME+IS-1))/
     1DFLOAT(IME-2)+XF(IT+1)*DFLOAT((ME+IS+2)*IM)/DFLOAT(IME+2))
      ALPHA=ALPHA*DFLOAT(IM)/DFLOAT(IT)
   30 CONTINUE
  300 IF(NFP-NF) 31,32,33
   32 IT=NF
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0*ALPHA/DFLOAT(ME2+IS2+1)
      SIM=SIM+XFP(IT)*COEFF*(XF(IT-1)*DFLOAT(IS)/DFLOAT(ME2+IS2-1))
      SJM=SJM+XFP(IT)*COEFF*(-XF(IT-1)*DFLOAT(IS*(ME+IS-1))/
     1DFLOAT(ME2+IS2-1))
      GO TO 31
   33 IT=NF
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0*ALPHA/DFLOAT(ME2+IS2+1)
      SIM=SIM+XFP(IT)*COEFF*(XF(IT-1)*DFLOAT(IS)/DFLOAT(ME2+IS2-1))
      SJM=SJM+XFP(IT)*COEFF*(-XF(IT-1)*DFLOAT(IS*(ME+IS-1))/
     1DFLOAT(ME2+IS2-1))
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)
      IT=NF+1
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0*ALPHA/DFLOAT(ME2+IS2+1)
      SIM=SIM+XFP(IT)*COEFF*(XF(IT-1)*DFLOAT(IS)/DFLOAT(ME2+IS2-1))
      SJM=SJM+XFP(IT)*COEFF*(-XF(IT-1)*DFLOAT(IS*(ME+IS-1))/
     1DFLOAT(ME2+IS2-1))
   31 CONTINUE
C
C        INTEGRALS OVER LAMBDA
C
      SIGPE=PE+PEP
      SG=0.D0
      SPG=0.D0
      SPH=0.D0
      SPI=0.D0
      SAG3 = 0.D0
      SAG2 = 0.D0
      DO 70 N=1,NXI
      ALAM=XI(N)/SIGPE+1.D0
      ALAM1=ALAM-1.D0
      ALAM2=ALAM+1.D0
      ALAMXI=ALAM1/ALAM2
      FONC=DEXP(-SIGPE*ALAM+(SIGMA+SIGMAP)*DLOG(ALAM2))
      TRUC=ALAM1*ALAM2
      TRIC=TRUC**ME
      TRUC=TRUC*TRIC
      FISC=FONC*TRIC
      FONC=FONC*TRUC
      YP=XGP(NGP)
      DO 71 JN=2,NGP
      NN=NGP-JN+1
   71 YP=YP*ALAMXI+XGP(NN)
      Y=XG(NG)
      DY=0.D0
      DO 72 JN=2,NG
      NN=NG-JN+1
      Y=Y*ALAMXI+XG(NN)
      DY=DY*ALAMXI+XG(NN+1)*DFLOAT(NN)
   72 CONTINUE
      DY=2.D0*DY
      DY=DY+Y*(ME*ALAM/ALAMXI+ALAM2*(SIGMA-PE*ALAM2))
      DY=DY/(ALAM2*ALAM2)
      FISC=FISC*Y*YP*OMEGA(N)
      FONC=FONC*YP*DY*OMEGA(N)
      SG=SG+FISC
      SPG=SPG+FISC*ALAM
      SAG2 =SAG2+FISC*ALAM*ALAM
      SAG3 =SAG3+FISC*ALAM*ALAM*ALAM
      SPI=SPI+FONC*ALAM
      SPH=SPH+FONC
   70 CONTINUE
      SG=SG/SIGPE
      SPG=SPG/SIGPE
      SPH=SPH/SIGPE
      SPI=SPI/SIGPE
      SAG2 = SAG2/SIGPE
      SAG3 = SAG3/SIGPE
C


      SIGPE=PE+PEP
      SPK=0.D0
      SPK2 = 0.D0
      DO 700 N=1,NXI
      ALAM=XI(N)/SIGPE+1.D0
      ALAM1=ALAM-1.D0
      ALAM2=ALAM+1.D0
      ALAMXI=ALAM1/ALAM2
      FONC=DEXP(-SIGPE*ALAM+(SIGMA+SIGMAP)*DLOG(ALAM2))
      TRUC=ALAM1*ALAM2
      TRIC=TRUC**ME
      TRUC=TRIC*TRUC
      FISC=FONC*TRIC*OMEGA(N)
      FONC=FONC*TRUC*OMEGA(N)
      YP=XGP(NGP)
      DO 710 JN=2,NGP
      NN=NGP-JN+1
  710 YP=YP*ALAMXI+XGP(NN)
      Y=XG(NG)
      DY=0.D0
      DO 720 JN=2,NG
      NN=NG-JN+1
      Y=Y*ALAMXI+XG(NN)
      DY=DY*ALAMXI+XG(NN+1)*DFLOAT(NN)
  720 CONTINUE
      DY=2.D0*DY
      DY=DY+Y*(ME*ALAM/ALAMXI+ALAM2*(SIGMA-PE*ALAM2))
      DY=DY/(ALAM2*ALAM2)
      FISC=FISC*Y*YP
      SPK=SPK+FONC*Y*YP
      SPK2 = SPK2+FONC*Y*YP*ALAM*ALAM
  700 CONTINUE
      SPK=SPK/SIGPE
      SPK2 = SPK2/SIGPE

      CRAD=-R*((1.D0+QU)*SUM*SPG+(QU-1.D0)*SIM*SG)/(4.D0*DELTA)
     1-R*R*(SUM*SPI+SG*SAM)/8.D0
      GAMRAD=-R*R*(SIM*SPH+SJM*SPG)/8.D0
      WRITE(111,66) SPG,SAG3
   66 FORMAT(4E20.9)
      CRAD=CRAD+OMEZUT*GAMRAD
      RETURN
      END



      SUBROUTINE ROTDIS(CROT,GAMROT)
C
C        CALCULATES THE ROTATIONAL MATRIX ELEMENTS.
C
      DOUBLE PRECISION A,AI,ALAM,ALAMXI,ALAM1,ALAM2,ALPHA,B,BI,C,CI,COEFACXY0572
     1,COEFF,CROT,DELTA,DEXP,DFLOAT,DLOG,DSQRT,DY,FISC,FONC,GAMROT,GDEU,
     2GPDEU,GPUN,GTROI,GUN,OMEZUT,PE,PEP,R,SAM,SG,SIGMA,SIGMAP,SIGPE,
     3SIM,SJM,SKM,SP,SPG,SPI,SPJ,SPK,SPN,SUM,SUMP,SUMP1,SUMP2,SUM1,SUM2,
     4TRIC,TRUC,XF,XFP,XG,XGP,Y,YP,SPK2
      A(I,I2)=-2.D0*DFLOAT(I*(ME2+I+2))/DFLOAT((ME2+I2-1)*(ME2+I2+1)*
     1(ME2+I2+3))
      B(I,I2)=+2.D0*DFLOAT((ME2+I+2)*(-4*ME*ME-4*ME*(I+2)-2*I-3))/
     1DFLOAT((ME2+I2+3)**2*(ME2+I2+5)*(ME2+I2+1))
      C(I,I2)=2.D0*DFLOAT((ME2+I+3)*(ME2+I+2))/DFLOAT((ME2+I2+7)
     1*(ME2+I2+3)*(ME2+I2+5))
      COMMON/LAG/XI(50),OMEGA(50),NXI
      COMMON/BRA/XFP(40),XGP(40),NFP,NGP
      COMMON/KET/XF(40),XG(40),NF,NG
      COMMON//R,PE,PEP,SIGMA,SIGMAP,DELTA,L,ME
      COMMON/ORIGE/OMEZUT
      ME2=2*ME
C
C        INTEGRAL OVER MU
C
      IALPHA=1
      DO 11 K=1,ME2
   11 IALPHA=IALPHA*K
C
C        INTEGRALS B AND B"
C
      ALPHA=DFLOAT(IALPHA)
      GDEU=-2.D0*ALPHA*DFLOAT((ME +2)*(ME2+1))/DFLOAT(ME2+3)
      GPDEU=-2.D0*ALPHA*DFLOAT(ME2+1)/DFLOAT(ME2+3)
      SP=0.D0
      IBMAX=(NFP-1)/2+1
      DO 25 IB=1,IBMAX
      IN=2*IB-1
   25 SP=SP+XFP(IN)
      SUM2=XF(2)*GDEU*SP
      SUMP2=XF(2)*GPDEU*SP
      ALPHA=ALPHA*DFLOAT(ME2+1)
      SUM1=0.D0
      SUMP1=0.D0
      NMAX=MIN0(NFP,NF-1)
      IF(NMAX.LT.2) GO TO 21
      DO 20 IT=2,NMAX
      IS=IT-1
      IS2=2*IS
      MIM=ME+IS
      IM=MIM+ME+1
      IME=IM+IS
      GUN=ALPHA*2.D0*DFLOAT((MIM-1)*IS)/DFLOAT(IME-2)
      GPUN=-ALPHA*DFLOAT(IS2)/DFLOAT(IME-2)
      GDEU=-ALPHA*2.D0*DFLOAT((MIM+2)*IM)/DFLOAT(IME+2)
      GPDEU=-2.D0*ALPHA*DFLOAT(IM)/DFLOAT(IME+2)
      SP=0.D0
      IBMAX=(NFP-IT)/2+1
      DO 22 IB=1,IBMAX
      IN=2*IB-2+IT
   22 SP=SP+XFP(IN)
      SUM1=SUM1+XF(IT-1)*GUN*SP
      SUMP1=SUMP1+XF(IT-1)*GPUN*SP
      SUM2=SUM2+XF(IT+1)*GDEU*SP
      SUMP2=SUMP2+XF(IT+1)*GPDEU*SP
      ALPHA=ALPHA*DFLOAT(IM)/DFLOAT(IT)
   20 CONTINUE
   21 IF(NFP-NF) 24,27,26
   27 IT=NF
      IS=IT-1
      IS2=2*IS
      GUN=ALPHA*2.D0*DFLOAT((ME+IS-1)*IS)/DFLOAT(ME2+IS2-1)
      GPUN=-ALPHA*DFLOAT(IS2)/DFLOAT(ME2+IS2-1)
      SP=XFP(IT)
      SUM1=SUM1+XF(IT-1)*GUN*SP
      SUMP1=SUMP1+XF(IT-1)*GPUN*SP
      GO TO 24
   26 IT=NF
      IS=IT-1
      IS2=2*IS
      GUN=ALPHA*2.D0*DFLOAT((ME+IS-1)*IS)/DFLOAT(ME2+IS2-1)
      GPUN=-ALPHA*DFLOAT(IS2)/DFLOAT(ME2+IS2-1)
      SP=0.D0
      IBMAX=(NFP-IT)/2+1
      DO 32 IB=1,IBMAX
      IN=2*IB-2+IT
   32 SP=SP+XFP(IN)
      SUM1=SUM1+XF(IT-1)*GUN*SP
      SUMP1=SUMP1+XF(IT-1)*GPUN*SP
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)
      IT=NF+1
      IS=IT-1
      IS2=2*IS
      GUN=ALPHA*2.D0*DFLOAT((ME+IS-1)*IS)/DFLOAT(ME2+IS2-1)
      GPUN=-ALPHA*DFLOAT(IS2)/DFLOAT(ME2+IS2-1)
      SP=0.D0
      IBMAX=(NFP-IT)/2+1
      DO 33 IB=1,IBMAX
      IN=2*IB-2+IT
   33 SP=SP+XFP(IN)
      SUM1=SUM1+XF(IT-1)*GUN*SP
      SUMP1=SUMP1+XF(IT-1)*GPUN*SP
   24 CONTINUE
      SUM=SUM1+SUM2
      SUMP=SUMP1+SUMP2
C
C        INTEGRAL B'
C
      ALPHA=DFLOAT(IALPHA)
      IF(NF.LT.4) XF(4)=0.D0
      IF(NF.LT.2) XF(2)=0.D0
      IF(NFP.LT.2) XFP(2)=0.D0
      COEF=DFLOAT(ME2+1)
      SAM=COEF*XFP(1)*(B(0,0)*XF(2)+C(0,0)*XF(4))*ALPHA
      ALPHA=ALPHA*DFLOAT(ME2+1)
      NMAX=MIN0(NFP,NF-3)
      IF(NMAX.LT.2) GO TO 59
      DO 60 IT=2,NMAX
      IS=IT-1
      MIM=ME+IS
      IM=ME+MIM+1
      IME=IM+IS
      AI=-2.D0*DFLOAT(IS*(IM+1))/DFLOAT((IME-2)*IME*(IME+2))
      BI=2.D0*DFLOAT((IM+1)*(-4*ME*ME-4*ME*(IS+2)-2*IS-3))/DFLOAT(
     1(IME+2)*(IME+2)*(IME+4)*IME)
      CI=2.D0*DFLOAT((IM+2)*(IM+1))/DFLOAT((IME+6)*(IME+4)*(IME+2))
      COEF=DFLOAT(IM)
      SAM=SAM+COEF*XFP(IT)*(AI*XF(IT-1)+BI*XF(IT+1)+CI*XF(IT+3))*ALPHA
      ALPHA=ALPHA*DFLOAT(IM)/DFLOAT(IT)
   60 CONTINUE
   59 IF(NF-3) 64,67,68
   68 IF(NMAX.EQ.NFP) GO TO 61
      IF(NFP-NF+1) 62,63,63
   62 CONTINUE
      IT=NF-2
      IS=IT-1
      IS2=2*IS
      COEF=DFLOAT(ME2+IS+1)
      SAM=SAM+COEF*XFP(IT)*(A(IS,IS2)*XF(IT-1)+B(IS,IS2)*XF(IT+1))
      GO TO 61
   63 IT=NF-2
      IS=IT-1
      IS2=2*IS
      COEF=DFLOAT(ME2+IS+1)
      SAM=SAM+COEF*XFP(IT)*(A(IS,IS2)*XF(IT-1)+B(IS,IS2)*XF(IT+1))
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)
   67 IT=NF-1
      IS=IT-1
      IS2=2*IS
      COEF=DFLOAT(ME2+IS+1)
      SAM=SAM+COEF*XFP(IT)*(A(IS,IS2)*XF(IT-1)+B(IS,IS2)*XF(IT+1))
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)
   64 IF(NFP-NF) 61,65,66
   65 IT=NF
      IS=IT-1
      IS2=2*IS
      COEF=DFLOAT(ME2+IS+1)
      SAM=SAM+COEF*XFP(IT)*(A(IS,IS2)*XF(IT-1))
      GO TO 61
   66 IT=NF
      IS=IT-1
      IS2=2*IS
      COEF=DFLOAT(ME2+IS+1)
      SAM=SAM+COEF*XFP(IT)*(A(IS,IS2)*XF(IT-1))
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)
      IT=NF+1
      IS=IT-1
      IS2=2*IS
      COEF=DFLOAT(ME2+IS+1)
      SAM=SAM+COEF*XFP(IT)*(A(IS,IS2)*XF(IT-1))
   61 CONTINUE
C
C        INTEGRAL B*
C
      SIM=0.D0
      ALPHA=DFLOAT(IALPHA*(ME2+1)*(ME2+2))
      NMAX=MIN0(NFP,NF-2)
      IF(NMAX.LT.1) GO TO 12
      DO 10 IT=1,NMAX
      IS=IT-1
      IME=ME2+IS+IS+1
      COEFF=2.D0/DFLOAT(IME+2)
      SIM=SIM+COEFF*XFP(IT)*ALPHA*(XF(IT+2)/DFLOAT(IME+4)
     1-XF(IT)/DFLOAT(IME))
      ALPHA=ALPHA*DFLOAT(ME2+IS+3)/DFLOAT(IT)
   10 CONTINUE
   12 IF(NFP-NF+1) 231,232,233
  232 IT=NF-1
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0/DFLOAT(ME2+IS2+3)
      SIM=SIM-COEFF*XFP(IT)*ALPHA*XF(IT)/DFLOAT(ME2+IS2+1)
      GO TO 231
  233 IT=NF-1
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0/DFLOAT(ME2+IS2+3)
      SIM=SIM-COEFF*XFP(IT)*ALPHA*XF(IT)/DFLOAT(ME2+IS2+1)
      ALPHA=ALPHA*DFLOAT(ME2+IS+3)/DFLOAT(IS+1)
      IT=NF
      IS=IT-1
      IS2=2*IS
      COEFF=2.D0/DFLOAT(ME2+IS2+3)
      SIM=SIM-COEFF*XFP(IT)*ALPHA*XF(IT)/DFLOAT(ME2+IS2+1)
  231 CONTINUE
C
C        INTEGRAL B"*
C
      SJM=0.D0
      ALPHA=DFLOAT(IALPHA)
      NMAX=MIN0(NFP,NF)
      DO 110 IT=1,NMAX
      IBMAX=(NFP-IT)/2+1
      SP=0.D0
      DO 111 IB=1,IBMAX
      IN=2*IB-2+IT
  111 SP=SP+XFP(IN)
      SJM=SJM-SP*2.D0*ALPHA*XF(IT)
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)
  110 CONTINUE
C
C        INTEGRAL B'*
C
      SKM=0.D0
      ALPHA=DFLOAT(IALPHA)
      IF(NF.LT.4) XF(4)=0.D0
      IF(NF.LT.3) XF(3)=0.D0
      DO 121 IT=1,2
      IS=IT-1
      IS2=2*IS
      COEFF=-2.D0*ALPHA
      GUN=     (DFLOAT((ME2+IS)*(ME+IS+1)*IS)/DFLOAT(ME2+IS2-1)
     1-DFLOAT((ME+IS)*(ME2+IS+1)*(IS+1))/DFLOAT(ME2+IS2+3))/DFLOAT(ME2
     2+IS2+1)
      GDEU=DFLOAT((ME2+IS+2)*(ME+IS+3)*(ME2+IS+1))/DFLOAT((ME2+IS2+3)
     1*(ME2+IS2+5))
      SP=0.D0
      IBMAX=(NFP-IT)/2+1
      DO 122 IB=1,IBMAX
      IN=2*IB-2+IT
  122 SP=SP+XFP(IN)
      SKM=SKM+COEFF*SP*(XF(IT)*GUN+XF(IT+2)*GDEU)
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)
  121 CONTINUE
      NMAX=MIN0(NFP,NF-2)
      IF(NMAX.LT.3) GO TO 125
      DO 123 IT=3,NMAX
      IS=IT-1
      COEFF=-2.D0*ALPHA
      MIM=IS+ME
      IM=MIM+ME+1
      IME=IM+IS
      GUN=(DFLOAT((IM-1)*(MIM+1)*IS)/DFLOAT(IME-2)
     1-DFLOAT(MIM*IM*IT)/DFLOAT(IME+2))/DFLOAT(IME)
      GDEU=DFLOAT((IM+1)*(MIM+3)*IM)/DFLOAT((IME+2)*(IME+4))
      GTROI=-DFLOAT((MIM-2)*(IS-1)*IS)/DFLOAT((IME-2)*(IME-4))
      SP=0.D0
      IBMAX=(NFP-IT)/2+1
      DO 124 IB=1,IBMAX
      IN=2*IB-2+IT
  124 SP=SP+XFP(IN)
      SKM=SKM+COEFF*SP*(XF(IT)*GUN+XF(IT+2)*GDEU+XF(IT-2)*GTROI)
      ALPHA=ALPHA*DFLOAT(IM)/DFLOAT(IT)
  123 CONTINUE
  125 IF(NMAX.EQ.NFP) GO TO 130
      NMAX=NMAX+1
      IF(NMAX.LT.3) NMAX=3
      NPMAX=MIN0(NFP,NF)
      IF(NPMAX.LT.NMAX) GO TO 132
      DO 131 IT=NMAX,NPMAX
      IS=IT-1
      IS2=2*IS
      COEFF=-2.D0*ALPHA
      GUN=     (DFLOAT((ME2+IS)*(ME+IS+1)*IS)/DFLOAT(ME2+IS2-1)
     1-DFLOAT((ME+IS)*(ME2+IS+1)*(IS+1))/DFLOAT(ME2+IS2+3))/DFLOAT(ME2
     2+IS2+1)
      GTROI=-DFLOAT((ME+IS-2)*(IS-1)*IS)/DFLOAT((ME2+IS2-1)*(ME2+IS2-3))
      SP=0.D0
      IBMAX=(NFP-IT)/2+1
      DO 134 IB=1,IBMAX
      IN=2*IB-2+IT
  134 SP=SP+XFP(IN)
      SKM=SKM+COEFF*SP*(XF(IT)*GUN+XF(IT-2)*GTROI)
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)
  131 CONTINUE
  132 IF(NPMAX.EQ.NFP) GO TO 130
      NMAX=MIN0(NFP,NF+2)
      NPMAX=NPMAX+1
      IF(NMAX.LT.NPMAX) GO TO 130
      DO 141  IT=NPMAX,NMAX
      IS=IT-1
      IS2=2*IS
      COEFF=-2.D0*ALPHA
      GTROI=-DFLOAT((ME+IS-2)*(IS-1)*IS)/DFLOAT((ME2+IS2-1)*(ME2+IS2-3))
      SP=0.D0
      IBMAX=(NFP-IT)/2+1
      DO 144 IB=1,IBMAX
      IN=2*IB-2+IT
  144 SP=SP+XFP(IN)
      SKM=SKM+COEFF*SP*XF(IT-2)*GTROI
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)
  141 CONTINUE
  130 CONTINUE
C
C        INTEGRALS OVER LAMBDA
C
      SIGPE=PE+PEP
      SG=0.D0
      SPG=0.D0
      SPI=0.D0
      SPJ=0.D0
      SPK=0.D0
      SPN=0.D0
      SPK2=0.D0

      DO 70 N=1,NXI
      ALAM=XI(N)/SIGPE+1.D0
      ALAM1=ALAM-1.D0
      ALAM2=ALAM+1.D0
      ALAMXI=ALAM1/ALAM2
      FONC=DEXP(-SIGPE*ALAM+(SIGMA+SIGMAP)*DLOG(ALAM2))
      TRUC=ALAM1*ALAM2
      TRIC=TRUC**ME
      TRUC=TRIC*TRUC
      FISC=FONC*TRIC*OMEGA(N)
      FONC=FONC*TRUC*OMEGA(N)
      YP=XGP(NGP)
      DO 71 JN=2,NGP
      NN=NGP-JN+1
   71 YP=YP*ALAMXI+XGP(NN)
      Y=XG(NG)
      DY=0.D0
      DO 72 JN=2,NG
      NN=NG-JN+1
      Y=Y*ALAMXI+XG(NN)
      DY=DY*ALAMXI+XG(NN+1)*DFLOAT(NN)
   72 CONTINUE
      DY=2.D0*DY
      DY=DY+Y*(ME*ALAM/ALAMXI+ALAM2*(SIGMA-PE*ALAM2))
      DY=DY/(ALAM2*ALAM2)
      FISC=FISC*Y*YP
      SG=SG+FONC*DY*YP
      SPG=SPG+FONC*Y*YP*ALAM
      SPI=SPI+FISC*ALAM
      SPJ=SPJ+FONC*DY*YP*ALAM
      SPK=SPK+FONC*Y*YP
      SPK2=SPK2+FONC*Y*YP*ALAM*ALAM
      SPN=SPN+FISC
   70 CONTINUE
      SG=SG/SIGPE
      SPG=SPG/SIGPE
      SPI=SPI/SIGPE
      SPJ=SPJ/SIGPE
      SPK=SPK/SIGPE
      SPN=SPN/SIGPE
      SPK2=SPK2/SIGPE
      WRITE(111,*) SPK, SPK2
  677 FORMAT(4E20.9)
C
      CROT=+R*R*R*(SG*SAM-SPG*SUM+ME*(SPI*SAM+SUMP*SPG))/16.D0
      GAMROT=+R*R*R*(SPJ*SIM-SPK*SKM+ME*(SPK*SJM+SIM*SPN))/16.D0
      CROT=CROT+OMEZUT*GAMROT
      IF(ME.NE.0) RETURN
      CROT=CROT*DSQRT(2.D0)
      GAMROT=GAMROT*DSQRT(2.D0)
      RETURN
      END
