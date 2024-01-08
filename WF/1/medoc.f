C                                                                       ACXY0004
C        MAIN PROGRAM FOR THE TEST RUN                                  ACXY0005
C                                                                       ACXY0006
      CALL MEDOC                                                        ACXY0007
      WRITE(6,1)                                                        ACXY0008
    1 FORMAT(1H1)                                                       ACXY0009
      STOP                                                              ACXY0010
      END                                                               ACXY0011
      SUBROUTINE MEDOC                                                  ACXY0012
C                                                                       ACXY0013
C        CALCULATION OF THE NORM AND COLLISION MATRIX ELEMENTS FOR OEDM.ACXY0014
C                                                                       ACXY0015
      DOUBLE PRECISION AO,CRAD,DELTA,DFLOAT,E,EP,FNOR,FNORP,GAMRAD,     ACXY0016
     1OMEZUT,PE,PEP,QU,R,RP,RR,SIGMA,SIGMAP,XF,XFP,XG,XGP               ACXY0017
      INTEGER INDEX1,INDEX2
      CHARACTER*260 FILENAME
      DIMENSION RR(100)                                                 ACXY0018
      DIMENSION N1(10),L1(10),M1(10)                                    ACXY0019
C                                                                       ACXY0020
C        COMMON/LAG/ CONTAINS THE ABCISSA AND WEIGHTS FOR GAUSS LAGUERREACXY0021
C     QUADRATURES.                                                      ACXY0022
C                                                                       ACXY0023
      COMMON/LAG/XI(50),OMEGA(50),NXI                                   ACXY0024
C                                                                       ACXY0025
C        COMMON BRA CONTAINS THE COEFFICIENTS OF THE SEMI-ANALYTICAL    ACXY0026
C     EXPANSION OF THE BRA VECTOR IN THE MATRIX ELEMENT.                ACXY0027
C                                                                       ACXY0028
      COMMON/BRA/XFP(40),XGP(40),NFP,NGP                                ACXY0029
C                                                                       ACXY0030
C        COMMON KET CONTAINS THE COEFFICIENTS OF THE SEMI-ANALYTICAL    ACXY0031
C     EXPANSION OF THE KET VECTOR IN THE MATRIX ELEMENT.                ACXY0032
C                                                                       ACXY0033
      COMMON/KET/XF(40),XG(40),NF,NG                                    ACXY0034
C                                                                       ACXY0035
C        THE FOLLOWING COMMONS ARE USED FOR THE TRANSMISSION OF         ACXY0036
C     VARIABLES TO THE SUBROUTINES                                      ACXY0037
C                                                                       ACXY0038
      COMMON//R,PE,PEP,SIGMA,SIGMAP,DELTA,L,ME                          ACXY0039
      COMMON/ORIGE/OMEZUT                                               ACXY0040
      OPEN(77,FILE='MEDOC_INPUT')
c      OPEN(99,FILE='IOEDM')
      OPEN(777,FILE='r.txt')
      OPEN(105,FILE='rot.txt')
      OPEN(106,FILE='rad.txt')


C                                                                       ACXY0041
C        THE UNIT IOEDM SHOULD BE THE SAME AS IN THE PROGRAM GRAVE.     ACXY0042
C                                                                       ACXY0043
    1 FORMAT (4E18.8)
    3 FORMAT (D12.6,4I3/(D20.14))                                       ACXY0046
    4 FORMAT (3I3)                                                       ACXY0047
    5 FORMAT (D20.10)                                                   ACXY0048
    6 FORMAT (44H MATRIX ELEMENTS ARE ZERO BETWEEN THE STATES,3(I3,2X),1ACXY0049
     1H/,3(I3,2X))                                                      ACXY0050
    9 FORMAT (D12.8)                                                    ACXY0051
   91 FORMAT (1X,6HZ2/Z1=,1PD20.10,5X,2HR=,D15.8,4HU.A.,3X,1H*,3I2,4H/DRACXY0052
     1/,3I2,2H*=,1PD15.8,3X,5HGAMA=,D15.8)                              ACXY0053
   92 FORMAT (1X,6HZ2/Z1=,1PD20.10,5X,2HR=,D15.8,4HU.A.,3X,1H*,3I2,4H/LYACXY0054
     1/,3I2,2H*=,1PD15.8,3X,5HGAMA=,D15.8)                              ACXY0055
   93 FORMAT (F6.3,6I2,4D15.8)                                          ACXY0056
   94 FORMAT (1H1,37HDISTANCE FROM CENTER 1 TO THE ORIGIN=,D20.10,////) ACXY0057
   34 FORMAT (E9.4,E15.7)

C                                                                       ACXY0058
C        INPUT FROM UNIT IN:                                            ACXY0059
C                                                                       ACXY0060
C        AO= DISTANCE OF THE ORIGIN OF THE ELECTRONIC COORDINATE SYSTEM ACXY0061
C     TO THE CHARGE Z1 IN UNIT OF R. (SEE EQ. 7 OF L.W.UP).             ACXY0062
C        NXI=NUMBER OF PIVOTS IN THE GAUSS-LAGUERRE INTEGRATION.        ACXY0063
C        XI,OMEGA= ABCISSAS AND WEIGHTS FOR THE GAUSS-LAGUERRE INTEGRA- ACXY0064
C     TION ( NXI SINGLE PRECISION VALUES).                              ACXY0065
C        JPERF= IF JPERF IS NON-ZERO, RESULTS ARE WRITTEN ON UNIT 7.    ACXY0066
C        NR= TWO POSSIBILITIES:                                         ACXY0067
C                                                                       ACXY0068
C     *****NR NON ZERO:                                                 ACXY0069
C        THEN NR IS THE NUMBER OF VALUES OF THE INTERNUCLEAR DISTANCE   ACXY0070
C     AND THE FOLLOWING QUANTITIES ARE ALSO READ ON UNIT IN:            ACXY0071
C        RR= ARRAY OF VALUES OF THE INTERNUCLEAR DISTANCE ( NR VALUES). ACXY0072
C        QU= RATIO OF NUCLEAR CHARGES.                                  ACXY0073
C        KN= NUMBER OF STATES ( KN.GE.2 ).                              ACXY0074
C        N1,L1,M1= ARRAYS OF QUANTUM NUMBERS OF THE MOLECULAR STATES IN ACXY0075
C     THE UNITED ATOM LIMIT ( KN VALUES).                               ACXY0076
C                                                                       ACXY0077
C     *****NR IS ZERO:                                                  ACXY0078
C        THEN THE ONLY FURTHER QUANTITY THAT SHOULD BE GIVEN ON UNIT IN ACXY0079
C     IS:                                                               ACXY0080
C        KN= NUMBER OF STATES (DEFAULT VALUE: 2).                       ACXY0081
C     AND THE PROGRAM READS THE OTHER DATA ON UNIT IOEDM FOR EACH STATE ACXY0082
C     (KN SETS OF DATA).                                                ACXY0083
C        QU= RATION OF NUCLEAR CHARGES                                  ACXY0084
C        N1,L1,M1= ARRAYS OF QUANTUM NUMBERS OF THE MOLECULAR STATES IN ACXY0085
C     THE UNITED ATOM LIMIT.                                            ACXY0086
C        NR= NUMBER OF VALUES OF R.                                     ACXY0087
C        RR= ARRAY OF VALUES OF R ( INTERNUCLEAR DISTANCE).             ACXY0088
C                                                                       ACXY0089
      READ(77,*) AO                                                     ACXY0090
      READ(77,*) NXI
      READ(77,*) (XI(N),OMEGA(N),N=1,NXI)
      READ(77,*) JPERF                                                  ACXY0092
      READ(777,*) NR                                                     ACXY0093
      IF(NR.EQ.0) GO TO 10                                              ACXY0094
      READ(777,*) (RR(N),N=1,NR)                                         ACXY0095
      READ(77,*) QU                                                     ACXY0096
      READ(77,*) KN
      READ(77,*) (N1(K),L1(K),M1(K),K=1,KN)
      READ(77,*) FILENAME

      GO TO 13                                                          ACXY0098
   10 READ(77,*) KN                                                     ACXY0099
      IF(KN.EQ.0) KN=2                                                  ACXY0100
      DO 12 K=1,KN                                                      ACXY0101
      READ(99,*) QU,N1(K),L1(K),M1(K),NR,(RR(N),N=1,NR)
   12 CONTINUE                                                          ACXY0103
   13 CONTINUE                                                          ACXY0104
C                                                                       ACXY0105
      OMEZUT=1.D0-2.D0*AO                                               ACXY0106
      WRITE(6,94) AO                                                    ACXY0107
      IF(JPERF.NE.0) WRITE(7,94) AO                                     ACXY0108
      DO 200 K=1,KN-1                                                   ACXY0109
      DO 300 KK=K+1,KN                                                  ACXY0110
      IF(IABS(M1(K)-M1(KK))-1) 15,18,11                                 ACXY0111
   18 CONTINUE                                                          ACXY0112
      IF(M1(K).GT.M1(KK)) GO TO 16                                      ACXY0113
   15 CONTINUE                                                          ACXY0114
      NN=N1(K)                                                          ACXY0115
      L=L1(K)                                                           ACXY0116
      ME=M1(K)                                                          ACXY0117
      INDEX1 = 100*NN+10*L+ME
      NNP=N1(KK)                                                        ACXY0118
      LP=L1(KK)                                                         ACXY0119
      MEP=M1(KK)                                                        ACXY0120
      INDEX2 = 100*NNP+10*LP+MEP
      GO TO 17                                                          ACXY0121
   16 CONTINUE                                                          ACXY0122
      NN=N1(KK)                                                         ACXY0123
      L=L1(KK)                                                          ACXY0124
      ME=M1(KK)                                                         ACXY0125
      INDEX1 = 100*NN+10*L+ME
      NNP=N1(K)                                                         ACXY0126
      LP=L1(K)                                                          ACXY0127
      MEP=M1(K)                                                         ACXY0128
      INDEX2 = 100*NNP+10*LP+MEP
   17 CONTINUE                                                          ACXY0129

C                                                                       ACXY0130
C        LOOP ON INTERNUCLEAR DISTANCES                                 ACXY0131
C                                                                       ACXY0132
      OPEN(111,FILE=FILENAME)

      DO 100 JN=1,NR                                                    ACXY0133
      R=RR(JN)                                                          ACXY0134
      RP=R                                                              ACXY0135
c     KET, 210
      CALL ENTREE(NN,L,ME,R,E,PE,NG,XG,NF,XF,QU,INDEX1)
c     BRA, 211
      CALL ENTREE(NNP,LP,MEP,R,EP,PEP,NGP,XGP,NFP,XFP,QU,INDEX2)
      SIGMA=R*(1.D0+QU)/(2.D0*PE)-1.D0-DFLOAT(ME)                       ACXY0138
      SIGMAP=RP*(1.D0+QU)/(2.D0*PEP)-1.D0-DFLOAT(MEP)                   ACXY0139
      IF(MEP.NE.ME) GO TO 14                                            ACXY0140
C                                                                       ACXY0141
C        RADIAL COUPLING CASE                                           ACXY0142
C                                                                       ACXY0143

      CALL NORDIF(ME,PE,SIGMA,NG,XG,NF,XF,FNOR,R)                       ACXY0144

      CALL NORDIF(MEP,PEP,SIGMAP,NGP,XGP,NFP,XFP,FNORP,RP)              ACXY0145
      WRITE(111,*) FNOR*FNORP
      DELTA=E-EP                                                        ACXY0146
      CALL RADAS(CRAD,QU,GAMRAD)                                        ACXY0147
      GAMRAD=GAMRAD*FNOR*FNORP                                          ACXY0148
      CRAD=CRAD*FNOR*FNORP                                              ACXY0149
      WRITE(6,91) QU,R,NNP,LP,MEP,NN,L,ME,CRAD,GAMRAD                   ACXY0150
      WRITE(106,34) R, CRAD
      IF(JPERF.NE.0) WRITE(7,93) QU,NNP,LP,MEP,NN,L,ME,R,CRAD,E,EP      ACXY0151
      GO TO 100                                                         ACXY0152
C                                                                       ACXY0153
C        ROTATIONAL COUPLING CASE.                                      ACXY0154
C                                                                       ACXY0155
   14 CONTINUE                                                          ACXY0156

      CALL NORDIF(ME,PE,SIGMA,NG,XG,NF,XF,FNOR,R)                       ACXY0157
      CALL NORDIF(MEP,PEP,SIGMAP,NGP,XGP,NFP,XFP,FNORP,R)               ACXY0158
      WRITE(111,*) FNOR*FNORP
      CALL ROTDIS(CRAD,GAMRAD)                                          ACXY0159

      GAMRAD=GAMRAD*FNOR*FNORP                                          ACXY0160
      CRAD=CRAD*FNOR*FNORP                                              ACXY0161
      WRITE(6,92) QU,R,NNP,LP,MEP,NN,L,ME,CRAD,GAMRAD                   ACXY0162
      WRITE(105,34) R,CRAD
      IF(JPERF.NE.0) WRITE(7,93) QU,NNP,LP,MEP,NN,L,ME,R,CRAD,E,EP      ACXY0163
C                                                                       ACXY0164
  100 CONTINUE                                                          ACXY0165
C                                                                       ACXY0166
C        OUTPUT ON UNIT 6:                                              ACXY0167
C        QU= RATIO OF NUCLEAR CHARGES.                                  ACXY0168
C        NNP,LP,MEP= QUANTUM NUMBERS DEFINING THE "BRA" STATE.          ACXY0169
C        NN,L,ME= QUANTUM NUMBERS DEFINING THE "KET" STATE.             ACXY0170
C        R= INTERNUCLEAR DISTANCE.                                      ACXY0171
C        CRAD= COUPLING MATRIX ELEMENT ( EQ. 7 OF THE L.W.UP).         0ACXY0172
C        GAMRAD= MATRIX ELEMENT OF THE OPERATOR D2 OR L2 (EQ. 9 OR 11 OFACXY0173
C     THE L.W.UP).                                                      ACXY0174
C                                                                       ACXY0175
C                                                                       ACXY0176
C        OUTPUT ON UNIT 7 ( IF JPERF IS NON-ZERO).                      ACXY0177
C        QU,NNP,LP,MEP,NN,L,ME,R,CRAD= AS ON UNIT 6.                    ACXY0178
C        E= ENERGY FOR THE "KET" STATE.                                 ACXY0179
C        EP= ENERGY FOR THE "BRA" STATE.                                ACXY0180
C                                                                       ACXY0181
      GO TO 300                                                         ACXY0182
   11 WRITE(6,6) N1(K),L1(K),M1(K),N1(KK),L1(KK),M1(KK)                 ACXY0183
  300 CONTINUE                                                          ACXY0184
  200 CONTINUE                                                          ACXY0185

      RETURN                                                            ACXY0186
      END                                                               ACXY0187
      SUBROUTINE ENTREE(N,L,ME,R,E,PE,NG,XG,NF,XF,QU,INDEX)
C                                                                       ACXY0189
C        SUBROUTINE USED TO READ THE CHARACTERISTICS OF ONE STATE       ACXY0190
C     THROUGH SUBROUTINE LECISA.                                        ACXY0191
C                                                                       ACXY0192
      DOUBLE PRECISION A,AOUT,E,EE,PE,PEPE,PU,QU,R,RR,XF,XG,XXF,XXG     ACXY0193
      INTEGER INDEX
      DIMENSION XG(1),XF(1)                                             ACXY0194
      DIMENSION AOUT(88)                                                ACXY0195
      EQUIVALENCE (PU,AOUT)                                             ACXY0196
      COMMON/TRAP/PU,RR,MME,NN,LL,IVER,EE,PEPE,A,NNG,NNF,XXG(40),XXF(40)ACXY0197
      DATA IVER/2/                                                      ACXY0198
      PU=QU                                                             ACXY0199
      NN=N                                                              ACXY0200
      LL=L                                                              ACXY0201
      MME=ME                                                            ACXY0202
      RR=R                                                              ACXY0203
      CALL LECISA(INDEX)                                                ACXY0204
      E=EE                                                              ACXY0205
      PE=PEPE                                                           ACXY0206
      NG=NNG                                                            ACXY0207
      NF=NNF                                                            ACXY0208
      DO 10 I=1,NG                                                      ACXY0209
   10 XG(I)=XXG(I)                                                      ACXY0210
      DO 20 I=1,NF                                                      ACXY0211
   20 XF(I)=XXF(I)                                                      ACXY0212
      RETURN                                                            ACXY0213
      END                                                               ACXY0214
      SUBROUTINE LECISA(INDEX)
C                                                                       ACXY0216
C                                                                       ACXY0217
C        USER SUPPLIED SUBROUTINE USED TO PICK UP THE CHARACTERISTICS OFACXY0218
C     THE STATE CONSIDERED AS CALCULATED IN THE PROGRAM GRAVE. SHOULD BEACXY0219
C      WRITTEN IN MACHINE LANGUAGE.                                     ACXY0220
C                                                                       ACXY0221
C                                                                       ACXY0222
C        THE VERSION GIVEN HERE IS ONLY SUITABLE FOR THE TEST RUN AND   ACXY0223
C     SHOULD NOT BE CONSIDERED AS AN INDICATION ON THE ROLE OF THE      ACXY0224
C     SUBROUTINE. SEE LONG WRITE UP FOR DETAILS ON THE CONCEPTION OF    ACXY0225
C     THE SUBROUTINE.                                                   ACXY0226
C                                                                       ACXY0227
      DOUBLE PRECISION AA,B,C,E,PE,A,      XG,XF                        ACXY0228
      INTEGER INDEX
      DIMENSION AA(1)                                                   ACXY0229
      COMMON/TRAP/B,C,M1,M2,M3,M4,E,PE,A,NG,NF,XG(40),XF(40)            ACXY0230
      READ(INDEX,*) E,PE,A,NG,NF
      READ(INDEX,*) (XG(I),I=1,NG)
      READ(INDEX,*) (XF(I),I=1,NF)
    1 FORMAT(3D20.14,2I3,/,(4D18.12))                                   ACXY0232
      RETURN                                                            ACXY0233
      END                                                               ACXY0234

      SUBROUTINE NORDIF(ME,PE,SIGMA,NG,XG,NF,XF,FNOR,R)                 ACXY0235
C                                                                       ACXY0236
C        CALCULATES THE NORM OF THE FUNCTION.                           ACXY0237
C                                                                       ACXY0238
C                                                                       ACXY0239
C     IMPORTANT REMARK: IN THE PROGRAM GRAVE, THE LAST TERM IN THE SEMI ACXY0240
C     ANALYTICAL EXPANSION IS SET EQUAL TO ONE. THE SUBROUTINE NORDIF   ACXY0241
C     FIRST TRANSFORMS THE COEFFICIENTS SO THAT THE FIRST NON-ZERO TERM ACXY0242
C     IS NOW EQUAL TO ONE. THE NORM IS CALCULATED FOR THIS NEW FORM     ACXY0243
C     OF EXPANSION.                                                     ACXY0244
      DOUBLE PRECISION ALAM,ALAMXI,ALAM1,ALAM2,ALPHA,BETA,DEXP,DFLOAT,  ACXY0245
     1DLOG,DSQRT,FNOR,FONC,PE,PP,R,SAM,SG,SIGMA,SIM,SPG,SUM,TROC,TRUC,  ACXY0246
     2XF,XF2,XG,XN,Y                                                    ACXY0247
      INTEGER NG,NF
      DIMENSION XG(NG),XF(NF)

      COMMON/LAG/XI(50),OMEGA(50),NXI                                   ACXY0249
      XN=XG(1)                                                          ACXY0250
      DO 16 K=1,NG                                                      ACXY0251
   16 XG(K)=XG(K)/XN                                                    ACXY0252
      XN=XF(1)                                                          ACXY0253
      IF(XN.EQ.0.D0) XN=XF(2)                                           ACXY0254
      DO 15 K=1,NF                                                      ACXY0255
   15 XF(K)=XF(K)/XN                                                    ACXY0256
C                                                                       ACXY0257
C        INTEGRALES SUR MU                                              ACXY0258
C                                                                       ACXY0259
      IALPHA=1                                                          ACXY0260
      ME2=2*ME                                                          ACXY0261
      DO 10 K=1,ME2                                                     ACXY0262
   10 IALPHA=IALPHA*K                                                   ACXY0263
      ALPHA=DFLOAT(IALPHA)                                              ACXY0264
      SAM=0.D0                                                          ACXY0265
      SUM=0.D0                                                          ACXY0266
      DO 11 IT=1,NF                                                     ACXY0267
      IS=IT-1                                                           ACXY0268
      XF2=XF(IT)*XF(IT)                                                 ACXY0269
      IM=ME2+IS                                                         ACXY0270
      IME=IM+IS+1                                                       ACXY0271
      SUM=SUM+ALPHA*2.D0*XF2/DFLOAT(IME)                                ACXY0272
      TROC=DFLOAT((IT)*(IM+1))/DFLOAT(IME+2)                            ACXY0273
      TROC=TROC+DFLOAT(IS*IM)/DFLOAT(IME-2)                             ACXY0274
      TROC=TROC/DFLOAT(IME*IME)                                         ACXY0275
      SAM=SAM+ALPHA*2.D0*XF2*TROC                                       ACXY0276
      ALPHA=ALPHA*DFLOAT(IM+1)/DFLOAT(IT)                               ACXY0277
   11 CONTINUE                                                          ACXY0278
      BETA=DFLOAT(IALPHA*(ME2+1)*(ME2+2))                               ACXY0279
      SIM=0.D0                                                          ACXY0280
      IF(NF.LE.2) GO TO 21                                              ACXY0281
      NFP=NF-2                                                          ACXY0282
      DO 20 IT=1,NFP                                                    ACXY0283
      IS=IT-1                                                           ACXY0284
      IME=ME2+IS+IS+1                                                   ACXY0285
      SIM=SIM+4.D0*XF(IT)*XF(IT+2)*BETA/DFLOAT(IME*(IME+2)*(IME+4))     ACXY0286
      BETA=BETA*DFLOAT(ME2+IS+3)/DFLOAT(IT)                             ACXY0287
   20 CONTINUE                                                          ACXY0288
   21 SAM=SAM+SIM                                                       ACXY0289
C                                                                       ACXY0290
C        INTEGRALES SUR LAMBDA                                          ACXY0291
C                                                                       ACXY0292
      SG=0.D0                                                           ACXY0293
      SPG=0.D0                                                          ACXY0294
      PP=2.D0*PE                                                        ACXY0295
      DO 30 N=1,NXI                                                     ACXY0296
c new variable as lambda
      ALAM=XI(N)/PP+1.D0                                                ACXY0297
      ALAM1=ALAM-1.D0                                                   ACXY0298
      ALAM2=ALAM+1.D0                                                   ACXY0299
      ALAMXI=ALAM1/ALAM2                                                ACXY0300
      FONC=DEXP(-PE*ALAM+SIGMA*DLOG(ALAM2))                             ACXY0301
c (FONC = (lam+1)^sig*exp(-p*lam))
      FONC=FONC*FONC                                                    ACXY0302
      TRUC=ALAM1*ALAM2                                                  ACXY0303
c TRUC = (lam+1)*(lam-1)
   35 TRUC=TRUC**ME                                                     ACXY0304
c ^m/2*^m/2 = ^m
      FONC=FONC*TRUC                                                    ACXY0305
c rwang Calculate terms before sum
c      Y=XG(NG)                                                          ACXY0306
c      DO 36 JN=2,NG                                                     ACXY0307
c      NN=NG-JN+1                                                        ACXY0308
c   36 Y=Y*ALAMXI+XG(NN)                                                 ACXY0309
      Y = 0.d0
      DO 36 JN =1,NG
   36 Y = Y+XG(JN)*ALAMXI**(JN-1)

      FONC=FONC*Y*Y*OMEGA(N)                                            ACXY0310
      SG=SG+FONC                                                        ACXY0311
      SPG=SPG+FONC*ALAM*ALAM                                            ACXY0312
   30 CONTINUE                                                          ACXY0313

      SG=SG/PP                                                          ACXY0314
      SPG=SPG/PP                                                        ACXY0315
C                                                                       ACXY0316

      FNOR=R*R*R*(SPG*SUM-SG*SAM)/8.D0                                  ACXY0317


      FNOR=DSQRT(1.D0/FNOR)                                             ACXY0318
      RETURN                                                            ACXY0319
      END                                                               ACXY0320




















      SUBROUTINE RADAS(CRAD,QU,GAMRAD)                                  ACXY0321
C                                                                       ACXY0322
C        CALCULATES THE RADIAL MATRIX ELEMENTS.                         ACXY0323
C                                                                       ACXY0324
      DOUBLE PRECISION A,AI,ALAM,ALAMXI,ALAM1,ALAM2,ALPHA,B,BI,C,CI,    ACXY0325
     1COEFF,CRAD,DELTA,DEXP,DFLOAT,DLOG,DY,FISC,FONC,GAMRAD,OMEZUT,     ACXY0326
     2PE,PEP,QU,R,SAM,SG,SIGMA,SIGMAP,SIGPE,SIM,SJM,SPG,SPH,SPI,SUM,    ACXY0327
     3TRIC,TRUC,XF,XFP,XG,XGP,Y,YP,SAG3,SAG2,SPK,SPK2









      A(I,I2)=(-DFLOAT((I+1)*(ME+I)*(ME2+I+1))/DFLOAT(ME2+I2+3)         ACXY0329
     1+DFLOAT((ME2+I)*(ME+I+1)*I)/DFLOAT(ME2+I2-1))/DFLOAT(ME2+I2+1)    ACXY0330
      B(I,I2)=DFLOAT((ME+I+3)*(ME2+I+2)*(ME2+I+1))/DFLOAT((ME2+I2+5)*   ACXY0331
     1(ME2+I2+3))                                                       ACXY0332
      C(I,I2)=-DFLOAT((I-1)*I*(ME+I-2))/DFLOAT((ME2+I2-3)*(ME2+I2-1))   ACXY0333
      COMMON/LAG/XI(50),OMEGA(50),NXI                                   ACXY0334
      COMMON/BRA/XFP(40),XGP(40),NFP,NGP                                ACXY0335
      COMMON/KET/XF(40),XG(40),NF,NG                                    ACXY0336
      COMMON//R,PE,PEP,SIGMA,SIGMAP,DELTA,L,ME                          ACXY0337
      COMMON/ORIGE/OMEZUT                                               ACXY0338
      ME2=2*ME                                                          ACXY0339
C                                                                       ACXY0340
C        INTEGRALS OVER MU                                              ACXY0341
C                                                                       ACXY0342
      IALPHA=1                                                          ACXY0343
      DO 11 K=1,ME2                                                     ACXY0344
   11 IALPHA=IALPHA*K                                                   ACXY034
C                                                                       ACXY0346
C        INTEGRAL E                                                     ACXY0347
C

      ALPHA=DFLOAT(IALPHA)                                              ACXY0349
      NMAX=MIN0(NFP,NF)                                                 ACXY0350
      SUM=0.D0                                                          ACXY0351
      DO 10 IT=1,NMAX                                                   ACXY0352
      IS=IT-1                                                           ACXY0353
      IM=ME2+IS+1                                                       ACXY0354
      SUM=SUM+2.D0*ALPHA*XF(IT)*XFP(IT)/DFLOAT(IM+IS)                   ACXY0355
      ALPHA=ALPHA*DFLOAT(IM)/DFLOAT(IT)                                 ACXY0356
   10 CONTINUE                                                          ACXY0357
C                                                                       ACXY0358
C        INTEGRAL  F                                                    ACXY0359
C                                                                       ACXY0360
      ALPHA=DFLOAT(IALPHA)                                              ACXY0361
      NMAX=MIN0(NFP,NF-2)                                               ACXY0362
      IF(NF.LT.3) XF(3)=0.D0                                            ACXY0363
      IS=0                                                              ACXY0364
      IS2=0                                                             ACXY0365
      COEFF=2.D0/DFLOAT(ME2+1)                                          ACXY0366
      SAM=XFP(1)*(A(IS,IS2)*XF(1)+B(IS,IS2)*XF(3))*COEFF*ALPHA          ACXY0367
      ALPHA=ALPHA*DFLOAT(ME2+1)                                         ACXY0368
      IF(NF.EQ.1) GO TO 28                                              ACXY0369
      IF(NFP.EQ.1) GO TO 21                                             ACXY0370
      IF(NF.LT.4) XF(4)=0.D0                                            ACXY0371
      IS=1                                                              ACXY0372
      IS2=2                                                             ACXY0373
      COEFF=2.D0/DFLOAT(ME2+3)                                          ACXY0374
      SAM=SAM+XFP(2)*(A(IS,IS2)*XF(2)+B(IS,IS2)*XF(4))*COEFF*ALPHA      ACXY0375
      ALPHA=ALPHA*DFLOAT(ME2+2)/2.D0                                    ACXY0376
      IF(NFP.EQ.2) GO TO 21                                             ACXY0377
      IF(NF.EQ.2) GO TO 27                                              ACXY0378
      IF(NF.EQ.3) GO TO 23                                              ACXY0379
      IF(NF.EQ.4) GO TO 19                                              ACXY0380
      DO 20 IT=3,NMAX                                                   ACXY0381
      IS=IT-1                                                           ACXY0382
      MIM=ME+IS                                                         ACXY0383
      IM=MIM+ME+1                                                       ACXY0384
      IME=IM+IS                                                         ACXY0385
      COEFF=2.D0/DFLOAT(IME)                                            ACXY0386
      AI=(-DFLOAT(IT*MIM*IM)/DFLOAT(IME+2)+DFLOAT((IM-1)*(MIM+1)*IS)/   ACXY0387
     1DFLOAT(IME-2))/DFLOAT(IME)                                        ACXY0388
      BI=DFLOAT((MIM+3)*(IM+1)*IM)/DFLOAT((IME+4)*(IME+2))              ACXY0389
      CI=-DFLOAT(IS*(IS-1)*(MIM-2))/DFLOAT((IME-4)*(IME-2))             ACXY0390
      TRUC=XFP(IT)*(AI*XF(IT)+BI*XF(IT+2)+CI*XF(IT-2))*ALPHA*COEFF      ACXY0391
      SAM=SAM+TRUC                                                      ACXY0392
      ALPHA=ALPHA*DFLOAT(IM)/DFLOAT(IT)                                 ACXY0393
   20 CONTINUE                                                          ACXY0394
   19 IF(NMAX.EQ.NFP) GO TO 21                                          ACXY0395
      IT=NMAX+1                                                         ACXY0396
      IF(IT.LE.2) GO TO 23                                              ACXY0397
      IS=IT-1                                                           ACXY0398
      IS2=2*IS                                                          ACXY0399
      COEFF=2.D0/DFLOAT(ME2+IS2+1)                                      ACXY0400
      TRUC=XFP(IT)*(A(IS,IS2)*XF(IT)+C(IS,IS2)*XF(IT-2))*ALPHA*COEFF    ACXY0401
      SAM=SAM+TRUC                                                      ACXY0402
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)                         ACXY0403
   23 IF(NFP.EQ.(NF-1)) GO TO 21                                        ACXY0404
      IF(NFP-NF-1) 24,25,26                                             ACXY0405
   24 IT=NF                                                             ACXY0406
      IS=IT-1                                                           ACXY0407
      IS2=2*IS                                                          ACXY0408
      COEFF=2.D0/DFLOAT(ME2+IS2+1)                                      ACXY0409
      TRUC=XFP(IT)*(A(IS,IS2)*XF(IT)+C(IS,IS2)*XF(IT-2))*ALPHA*COEFF    ACXY0410
      SAM=SAM+TRUC                                                      ACXY0411
      GO TO 21                                                          ACXY0412
   28 IF(NFP.LT.3) GO TO 21                                             ACXY0413
      IT=3                                                              ACXY0414
      IS=IT-1                                                           ACXY0415
      IS2=2*IS                                                          ACXY0416
      ALPHA=ALPHA*DFLOAT(ME2+2)/2.D0                                    ACXY0417
      COEFF=2.D0/DFLOAT(ME2+IS2+1)                                      ACXY0418
      TRUC=XFP(IT)*C(IS,IS2)*XF(IT-2)*ALPHA*COEFF                       ACXY0419
      SAM=SAM+TRUC                                                      ACXY0420
      GO TO 21                                                          ACXY0421
   27 IT=3                                                              ACXY0422
      IS=IT-1                                                           ACXY0423
      IS2=2*IS                                                          ACXY0424
      COEFF=2.D0/DFLOAT(ME2+IS2+1)                                      ACXY0425
      TRUC=XFP(IT)*C(IS,IS2)*XF(IT-2)*ALPHA*COEFF                       ACXY0426
      SAM=SAM+TRUC                                                      ACXY0427
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)                         ACXY0428
      IF(NFP.EQ.3) GO TO 21                                             ACXY0429
      IT=4                                                              ACXY0430
      IS=IT-1                                                           ACXY0431
      IS2=2*IS                                                          ACXY0432
      COEFF=2.D0/DFLOAT(ME2+IS2+1)                                      ACXY0433
      TRUC=XFP(IT)*C(IS,IS2)*XF(IT-2)*ALPHA*COEFF                       ACXY0434
      SAM=SAM+TRUC                                                      ACXY0435
      GO TO 21                                                          ACXY0436
   25 IT=NF                                                             ACXY0437
      IS=IT-1                                                           ACXY0438
      IS2=2*IS                                                          ACXY0439
      COEFF=2.D0/DFLOAT(ME2+IS2+1)                                      ACXY0440
      TRUC=XFP(IT)*(A(IS,IS2)*XF(IT)+C(IS,IS2)*XF(IT-2))*ALPHA*COEFF    ACXY0441
      SAM=SAM+TRUC                                                      ACXY0442
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)                         ACXY0443
      IT=NF+1                                                           ACXY0444
      IS=IT-1                                                           ACXY0445
      IS2=2*IS                                                          ACXY0446
      COEFF=2.D0/DFLOAT(ME2+IS2+1)                                      ACXY0447
      TRUC=XFP(IT)*C(IS,IS2)*XF(IT-2)*ALPHA*COEFF                       ACXY0448
      SAM=SAM+TRUC                                                      ACXY0449
      GO TO 21                                                          ACXY0450
   26 IT=NF                                                             ACXY0451
      IS=IT-1                                                           ACXY0452
      IS2=2*IS                                                          ACXY0453
      COEFF=2.D0/DFLOAT(ME2+IS2+1)                                      ACXY0454
      TRUC=XFP(IT)*(A(IS,IS2)*XF(IT)+C(IS,IS2)*XF(IT-2))*ALPHA*COEFF    ACXY0455
      SAM=SAM+TRUC                                                      ACXY0456
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)                         ACXY0457
      IT=NF+1                                                           ACXY0458
      IS=IT-1                                                           ACXY0459
      IS2=2*IS                                                          ACXY0460
      COEFF=2.D0/DFLOAT(ME2+IS2+1)                                      ACXY0461
      TRUC=XFP(IT)*C(IS,IS2)*XF(IT-2)*ALPHA*COEFF                       ACXY0462
      SAM=SAM+TRUC                                                      ACXY0463
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)                         ACXY0464
      IT=NF+2                                                           ACXY0465
      IS=IT-1                                                           ACXY0466
      IS2=2*IS                                                          ACXY0467
      COEFF=2.D0/DFLOAT(ME2+IS2+1)                                      ACXY0468
      TRUC=XFP(IT)*C(IS,IS2)*XF(IT-2)*ALPHA*COEFF                       ACXY0469
      SAM=SAM+TRUC                                                      ACXY0470
   21 CONTINUE                                                          ACXY0471
C                                                                       ACXY0472
C        INTEGRALS A AND B                                              ACXY0473
C                                                                       ACXY0474
      ALPHA=DFLOAT(IALPHA)                                              ACXY0475
      SIM=2.D0*ALPHA*XFP(1)*XF(2)/DFLOAT(ME2+3)                         ACXY0476
      SJM=2.D0*ALPHA*XFP(1)*XF(2)*DFLOAT(ME+2)/DFLOAT(ME2+3)            ACXY0477
      ALPHA=ALPHA*DFLOAT(ME2+1)                                         ACXY0478
      NMAX=MIN0(NFP,NF-1)                                               ACXY047
      IF(NMAX.EQ.1) GO TO 300                                           ACXY0480
      DO 30 IT=2,NMAX                                                   ACXY0481
      IS=IT-1                                                           ACXY0482
      IM=ME2+IS+1                                                       ACXY0483
      IME=IM+IS                                                         ACXY0484
      COEFF=2.D0*ALPHA/DFLOAT(IME)                                      ACXY0485
      SIM=SIM+XFP(IT)*COEFF*(XF(IT-1)*DFLOAT(IS)/DFLOAT(IME-2)+         ACXY0486
     1XF(IT+1)*DFLOAT(IM)/DFLOAT(IME+2))                                ACXY0487
      SJM=SJM+XFP(IT)*COEFF*(-XF(IT-1)*DFLOAT(IS*(ME+IS-1))/            ACXY0488
     1DFLOAT(IME-2)+XF(IT+1)*DFLOAT((ME+IS+2)*IM)/DFLOAT(IME+2))        ACXY0489
      ALPHA=ALPHA*DFLOAT(IM)/DFLOAT(IT)                                 ACXY0490
   30 CONTINUE                                                          ACXY0491
  300 IF(NFP-NF) 31,32,33                                               ACXY0492
   32 IT=NF                                                             ACXY0493
      IS=IT-1                                                           ACXY0494
      IS2=2*IS                                                          ACXY0495
      COEFF=2.D0*ALPHA/DFLOAT(ME2+IS2+1)                                ACXY0496
      SIM=SIM+XFP(IT)*COEFF*(XF(IT-1)*DFLOAT(IS)/DFLOAT(ME2+IS2-1))     ACXY0497
      SJM=SJM+XFP(IT)*COEFF*(-XF(IT-1)*DFLOAT(IS*(ME+IS-1))/            ACXY0498
     1DFLOAT(ME2+IS2-1))                                                ACXY0499
      GO TO 31                                                          ACXY0500
   33 IT=NF                                                             ACXY0501
      IS=IT-1                                                           ACXY0502
      IS2=2*IS                                                          ACXY0503
      COEFF=2.D0*ALPHA/DFLOAT(ME2+IS2+1)                                ACXY0504
      SIM=SIM+XFP(IT)*COEFF*(XF(IT-1)*DFLOAT(IS)/DFLOAT(ME2+IS2-1))     ACXY0505
      SJM=SJM+XFP(IT)*COEFF*(-XF(IT-1)*DFLOAT(IS*(ME+IS-1))/            ACXY0506
     1DFLOAT(ME2+IS2-1))                                                ACXY0507
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)                         ACXY0508
      IT=NF+1                                                           ACXY0509
      IS=IT-1                                                           ACXY0510
      IS2=2*IS                                                          ACXY0511
      COEFF=2.D0*ALPHA/DFLOAT(ME2+IS2+1)                                ACXY0512
      SIM=SIM+XFP(IT)*COEFF*(XF(IT-1)*DFLOAT(IS)/DFLOAT(ME2+IS2-1))     ACXY0513
      SJM=SJM+XFP(IT)*COEFF*(-XF(IT-1)*DFLOAT(IS*(ME+IS-1))/            ACXY0514
     1DFLOAT(ME2+IS2-1))                                                ACXY0515
   31 CONTINUE                                                          ACXY0516
C                                                                       ACXY0517
C        INTEGRALS OVER LAMBDA                                          ACXY0518
C                                                                       ACXY0519
      SIGPE=PE+PEP                                                      ACXY0520
      SG=0.D0                                                           ACXY0521
      SPG=0.D0                                                          ACXY0522
      SPH=0.D0                                                          ACXY0523
      SPI=0.D0                                                          ACXY0524
      SAG3 = 0.D0
      SAG2 = 0.D0
      DO 70 N=1,NXI                                                     ACXY0525
      ALAM=XI(N)/SIGPE+1.D0                                             ACXY0526
      ALAM1=ALAM-1.D0                                                   ACXY0527
      ALAM2=ALAM+1.D0                                                   ACXY0528
      ALAMXI=ALAM1/ALAM2                                                ACXY0529
      FONC=DEXP(-SIGPE*ALAM+(SIGMA+SIGMAP)*DLOG(ALAM2))                 ACXY0530
      TRUC=ALAM1*ALAM2                                                  ACXY0531
      TRIC=TRUC**ME                                                     ACXY0532
      TRUC=TRUC*TRIC                                                    ACXY0533
      FISC=FONC*TRIC                                                    ACXY0534
      FONC=FONC*TRUC                                                    ACXY0535
      YP=XGP(NGP)                                                       ACXY0536
      DO 71 JN=2,NGP                                                    ACXY0537
      NN=NGP-JN+1                                                       ACXY0538
   71 YP=YP*ALAMXI+XGP(NN)                                              ACXY0539
      Y=XG(NG)                                                          ACXY0540
      DY=0.D0                                                           ACXY0541
      DO 72 JN=2,NG                                                     ACXY0542
      NN=NG-JN+1                                                        ACXY0543
      Y=Y*ALAMXI+XG(NN)                                                 ACXY0544
      DY=DY*ALAMXI+XG(NN+1)*DFLOAT(NN)                                  ACXY0545
   72 CONTINUE                                                          ACXY0546
      DY=2.D0*DY                                                        ACXY0547
      DY=DY+Y*(ME*ALAM/ALAMXI+ALAM2*(SIGMA-PE*ALAM2))                   ACXY0548
      DY=DY/(ALAM2*ALAM2)                                               ACXY0549
      FISC=FISC*Y*YP*OMEGA(N)                                           ACXY0550
      FONC=FONC*YP*DY*OMEGA(N)                                          ACXY0551
      SG=SG+FISC                                                        ACXY0552
      SPG=SPG+FISC*ALAM                                                 ACXY0553
      SAG2 =SAG2+FISC*ALAM*ALAM
      SAG3 =SAG3+FISC*ALAM*ALAM*ALAM
      SPI=SPI+FONC*ALAM                                                 ACXY0554
      SPH=SPH+FONC                                                      ACXY0555
   70 CONTINUE                                                          ACXY0556
      SG=SG/SIGPE                                                       ACXY0557
      SPG=SPG/SIGPE                                                     ACXY0558
      SPH=SPH/SIGPE                                                     ACXY0559
      SPI=SPI/SIGPE                                                     ACXY0560
      SAG2 = SAG2/SIGPE
      SAG3 = SAG3/SIGPE
C                                                                       ACXY0561


      SIGPE=PE+PEP                                                      ACXY0873
      SPK=0.D0                                                          ACXY0878
      SPK2 = 0.D0
      DO 700 N=1,NXI                                                    ACXY0880
      ALAM=XI(N)/SIGPE+1.D0                                             ACXY0881
      ALAM1=ALAM-1.D0                                                   ACXY0882
      ALAM2=ALAM+1.D0                                                   ACXY0883
      ALAMXI=ALAM1/ALAM2                                                ACXY0884
      FONC=DEXP(-SIGPE*ALAM+(SIGMA+SIGMAP)*DLOG(ALAM2))                 ACXY0885
      TRUC=ALAM1*ALAM2                                                  ACXY0886
      TRIC=TRUC**ME                                                     ACXY0887
      TRUC=TRIC*TRUC                                                    ACXY0888
      FISC=FONC*TRIC*OMEGA(N)                                           ACXY0889
      FONC=FONC*TRUC*OMEGA(N)                                           ACXY0890
      YP=XGP(NGP)                                                       ACXY0891
      DO 710 JN=2,NGP                                                    ACXY0892
      NN=NGP-JN+1                                                       ACXY0893
  710 YP=YP*ALAMXI+XGP(NN)                                              ACXY0894
      Y=XG(NG)                                                          ACXY0895
      DY=0.D0                                                           ACXY0896
      DO 720 JN=2,NG                                                     ACXY0897
      NN=NG-JN+1                                                        ACXY0898
      Y=Y*ALAMXI+XG(NN)                                                 ACXY0899
      DY=DY*ALAMXI+XG(NN+1)*DFLOAT(NN)                                  ACXY0900
  720 CONTINUE                                                          ACXY0901
      DY=2.D0*DY                                                        ACXY0902
      DY=DY+Y*(ME*ALAM/ALAMXI+ALAM2*(SIGMA-PE*ALAM2))                   ACXY0903
      DY=DY/(ALAM2*ALAM2)                                               ACXY0904
      FISC=FISC*Y*YP                                                    ACXY0905
      SPK=SPK+FONC*Y*YP                                                 ACXY0910
      SPK2 = SPK2+FONC*Y*YP*ALAM*ALAM
  700 CONTINUE                                                          ACXY0912
      SPK=SPK/SIGPE                                                     ACXY0917
      SPK2 = SPK2/SIGPE

      CRAD=-R*((1.D0+QU)*SUM*SPG+(QU-1.D0)*SIM*SG)/(4.D0*DELTA)         ACXY0562
     1-R*R*(SUM*SPI+SG*SAM)/8.D0                                        ACXY0563
      GAMRAD=-R*R*(SIM*SPH+SJM*SPG)/8.D0                                ACXY0564
      WRITE(111,66) SPG,SAG3
   66 FORMAT(4E20.9)
      CRAD=CRAD+OMEZUT*GAMRAD
      RETURN                                                            ACXY0566
      END                                                               ACXY0567




















      SUBROUTINE ROTDIS(CROT,GAMROT)                                    ACXY0568
C                                                                       ACXY0569
C        CALCULATES THE ROTATIONAL MATRIX ELEMENTS.                     ACXY0570
C                                                                       ACXY0571
      DOUBLE PRECISION A,AI,ALAM,ALAMXI,ALAM1,ALAM2,ALPHA,B,BI,C,CI,COEFACXY0572
     1,COEFF,CROT,DELTA,DEXP,DFLOAT,DLOG,DSQRT,DY,FISC,FONC,GAMROT,GDEU,ACXY0573
     2GPDEU,GPUN,GTROI,GUN,OMEZUT,PE,PEP,R,SAM,SG,SIGMA,SIGMAP,SIGPE,   ACXY0574
     3SIM,SJM,SKM,SP,SPG,SPI,SPJ,SPK,SPN,SUM,SUMP,SUMP1,SUMP2,SUM1,SUM2,ACXY0575
     4TRIC,TRUC,XF,XFP,XG,XGP,Y,YP,SPK2
      A(I,I2)=-2.D0*DFLOAT(I*(ME2+I+2))/DFLOAT((ME2+I2-1)*(ME2+I2+1)*   ACXY0577
     1(ME2+I2+3))                                                       ACXY0578
      B(I,I2)=+2.D0*DFLOAT((ME2+I+2)*(-4*ME*ME-4*ME*(I+2)-2*I-3))/      ACXY0579
     1DFLOAT((ME2+I2+3)**2*(ME2+I2+5)*(ME2+I2+1))                       ACXY0580
      C(I,I2)=2.D0*DFLOAT((ME2+I+3)*(ME2+I+2))/DFLOAT((ME2+I2+7)        ACXY0581
     1*(ME2+I2+3)*(ME2+I2+5))                                           ACXY0582
      COMMON/LAG/XI(50),OMEGA(50),NXI                                   ACXY0583
      COMMON/BRA/XFP(40),XGP(40),NFP,NGP                                ACXY0584
      COMMON/KET/XF(40),XG(40),NF,NG                                    ACXY0585
      COMMON//R,PE,PEP,SIGMA,SIGMAP,DELTA,L,ME                          ACXY0586
      COMMON/ORIGE/OMEZUT                                               ACXY0587
      ME2=2*ME                                                          ACXY0588
C                                                                       ACXY0589
C        INTEGRAL OVER MU                                               ACXY0590
C                                                                       ACXY0591
      IALPHA=1                                                          ACXY0592
      DO 11 K=1,ME2                                                     ACXY0593
   11 IALPHA=IALPHA*K                                                   ACXY0594
C                                                                       ACXY0595
C        INTEGRALS B AND B"                                             ACXY0596
C                                                                       ACXY0597
      ALPHA=DFLOAT(IALPHA)                                              ACXY0598
      GDEU=-2.D0*ALPHA*DFLOAT((ME +2)*(ME2+1))/DFLOAT(ME2+3)            ACXY0599
      GPDEU=-2.D0*ALPHA*DFLOAT(ME2+1)/DFLOAT(ME2+3)                     ACXY0600
      SP=0.D0                                                           ACXY0601
      IBMAX=(NFP-1)/2+1                                                 ACXY0602
      DO 25 IB=1,IBMAX                                                  ACXY0603
      IN=2*IB-1                                                         ACXY0604
   25 SP=SP+XFP(IN)                                                     ACXY0605
      SUM2=XF(2)*GDEU*SP                                                ACXY0606
      SUMP2=XF(2)*GPDEU*SP                                              ACXY0607
      ALPHA=ALPHA*DFLOAT(ME2+1)                                         ACXY0608
      SUM1=0.D0                                                         ACXY0609
      SUMP1=0.D0                                                        ACXY0610
      NMAX=MIN0(NFP,NF-1)                                               ACXY0611
      IF(NMAX.LT.2) GO TO 21                                            ACXY0612
      DO 20 IT=2,NMAX                                                   ACXY0613
      IS=IT-1                                                           ACXY0614
      IS2=2*IS                                                          ACXY0615
      MIM=ME+IS                                                         ACXY0616
      IM=MIM+ME+1                                                       ACXY0617
      IME=IM+IS                                                         ACXY0618
      GUN=ALPHA*2.D0*DFLOAT((MIM-1)*IS)/DFLOAT(IME-2)                   ACXY0619
      GPUN=-ALPHA*DFLOAT(IS2)/DFLOAT(IME-2)                             ACXY0620
      GDEU=-ALPHA*2.D0*DFLOAT((MIM+2)*IM)/DFLOAT(IME+2)                 ACXY0621
      GPDEU=-2.D0*ALPHA*DFLOAT(IM)/DFLOAT(IME+2)                        ACXY0622
      SP=0.D0                                                           ACXY0623
      IBMAX=(NFP-IT)/2+1                                                ACXY0624
      DO 22 IB=1,IBMAX                                                  ACXY0625
      IN=2*IB-2+IT                                                      ACXY0626
   22 SP=SP+XFP(IN)                                                     ACXY0627
      SUM1=SUM1+XF(IT-1)*GUN*SP                                         ACXY0628
      SUMP1=SUMP1+XF(IT-1)*GPUN*SP                                      ACXY0629
      SUM2=SUM2+XF(IT+1)*GDEU*SP                                        ACXY0630
      SUMP2=SUMP2+XF(IT+1)*GPDEU*SP                                     ACXY0631
      ALPHA=ALPHA*DFLOAT(IM)/DFLOAT(IT)                                 ACXY0632
   20 CONTINUE                                                          ACXY0633
   21 IF(NFP-NF) 24,27,26                                               ACXY0634
   27 IT=NF                                                             ACXY0635
      IS=IT-1                                                           ACXY0636
      IS2=2*IS                                                          ACXY0637
      GUN=ALPHA*2.D0*DFLOAT((ME+IS-1)*IS)/DFLOAT(ME2+IS2-1)             ACXY0638
      GPUN=-ALPHA*DFLOAT(IS2)/DFLOAT(ME2+IS2-1)                         ACXY0639
      SP=XFP(IT)                                                        ACXY0640
      SUM1=SUM1+XF(IT-1)*GUN*SP                                         ACXY0641
      SUMP1=SUMP1+XF(IT-1)*GPUN*SP                                      ACXY0642
      GO TO 24                                                          ACXY0643
   26 IT=NF                                                             ACXY0644
      IS=IT-1                                                           ACXY0645
      IS2=2*IS                                                          ACXY0646
      GUN=ALPHA*2.D0*DFLOAT((ME+IS-1)*IS)/DFLOAT(ME2+IS2-1)             ACXY0647
      GPUN=-ALPHA*DFLOAT(IS2)/DFLOAT(ME2+IS2-1)                         ACXY0648
      SP=0.D0                                                           ACXY0649
      IBMAX=(NFP-IT)/2+1                                                ACXY0650
      DO 32 IB=1,IBMAX                                                  ACXY0651
      IN=2*IB-2+IT                                                      ACXY0652
   32 SP=SP+XFP(IN)                                                     ACXY0653
      SUM1=SUM1+XF(IT-1)*GUN*SP                                         ACXY0654
      SUMP1=SUMP1+XF(IT-1)*GPUN*SP                                      ACXY0655
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)                         ACXY0656
      IT=NF+1                                                           ACXY0657
      IS=IT-1                                                           ACXY0658
      IS2=2*IS                                                          ACXY0659
      GUN=ALPHA*2.D0*DFLOAT((ME+IS-1)*IS)/DFLOAT(ME2+IS2-1)             ACXY0660
      GPUN=-ALPHA*DFLOAT(IS2)/DFLOAT(ME2+IS2-1)                         ACXY0661
      SP=0.D0                                                           ACXY0662
      IBMAX=(NFP-IT)/2+1                                                ACXY0663
      DO 33 IB=1,IBMAX                                                  ACXY0664
      IN=2*IB-2+IT                                                      ACXY0665
   33 SP=SP+XFP(IN)                                                     ACXY0666
      SUM1=SUM1+XF(IT-1)*GUN*SP                                         ACXY0667
      SUMP1=SUMP1+XF(IT-1)*GPUN*SP                                      ACXY0668
   24 CONTINUE                                                          ACXY0669
      SUM=SUM1+SUM2                                                     ACXY0670
      SUMP=SUMP1+SUMP2                                                  ACXY0671
C                                                                       ACXY0672
C        INTEGRAL B'                                                    ACXY0673
C                                                                       ACXY0674
      ALPHA=DFLOAT(IALPHA)                                              ACXY0675
      IF(NF.LT.4) XF(4)=0.D0                                            ACXY0676
      IF(NF.LT.2) XF(2)=0.D0                                            ACXY0677
      IF(NFP.LT.2) XFP(2)=0.D0                                          ACXY0678
      COEF=DFLOAT(ME2+1)                                                ACXY0679
      SAM=COEF*XFP(1)*(B(0,0)*XF(2)+C(0,0)*XF(4))*ALPHA                 ACXY0680
      ALPHA=ALPHA*DFLOAT(ME2+1)                                         ACXY0681
      NMAX=MIN0(NFP,NF-3)                                               ACXY0682
      IF(NMAX.LT.2) GO TO 59                                            ACXY0683
      DO 60 IT=2,NMAX                                                   ACXY0684
      IS=IT-1                                                           ACXY0685
      MIM=ME+IS                                                         ACXY0686
      IM=ME+MIM+1                                                       ACXY0687
      IME=IM+IS                                                         ACXY0688
      AI=-2.D0*DFLOAT(IS*(IM+1))/DFLOAT((IME-2)*IME*(IME+2))            ACXY0689
      BI=2.D0*DFLOAT((IM+1)*(-4*ME*ME-4*ME*(IS+2)-2*IS-3))/DFLOAT(      ACXY0690
     1(IME+2)*(IME+2)*(IME+4)*IME)                                      ACXY0691
      CI=2.D0*DFLOAT((IM+2)*(IM+1))/DFLOAT((IME+6)*(IME+4)*(IME+2))     ACXY0692
      COEF=DFLOAT(IM)                                                   ACXY0693
      SAM=SAM+COEF*XFP(IT)*(AI*XF(IT-1)+BI*XF(IT+1)+CI*XF(IT+3))*ALPHA  ACXY0694
      ALPHA=ALPHA*DFLOAT(IM)/DFLOAT(IT)                                 ACXY0695
   60 CONTINUE                                                          ACXY0696
   59 IF(NF-3) 64,67,68                                                 ACXY0697
   68 IF(NMAX.EQ.NFP) GO TO 61                                          ACXY0698
      IF(NFP-NF+1) 62,63,63                                             ACXY0699
   62 CONTINUE                                                          ACXY0700
      IT=NF-2                                                           ACXY0701
      IS=IT-1                                                           ACXY0702
      IS2=2*IS                                                          ACXY0703
      COEF=DFLOAT(ME2+IS+1)                                             ACXY0704
      SAM=SAM+COEF*XFP(IT)*(A(IS,IS2)*XF(IT-1)+B(IS,IS2)*XF(IT+1))      ACXY0705
      GO TO 61                                                          ACXY0706
   63 IT=NF-2                                                           ACXY0707
      IS=IT-1                                                           ACXY0708
      IS2=2*IS                                                          ACXY0709
      COEF=DFLOAT(ME2+IS+1)                                             ACXY0710
      SAM=SAM+COEF*XFP(IT)*(A(IS,IS2)*XF(IT-1)+B(IS,IS2)*XF(IT+1))      ACXY0711
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)                         ACXY0712
   67 IT=NF-1                                                           ACXY0713
      IS=IT-1                                                           ACXY0714
      IS2=2*IS                                                          ACXY0715
      COEF=DFLOAT(ME2+IS+1)                                             ACXY0716
      SAM=SAM+COEF*XFP(IT)*(A(IS,IS2)*XF(IT-1)+B(IS,IS2)*XF(IT+1))      ACXY0717
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)                         ACXY0718
   64 IF(NFP-NF) 61,65,66                                               ACXY0719
   65 IT=NF                                                             ACXY0720
      IS=IT-1                                                           ACXY0721
      IS2=2*IS                                                          ACXY0722
      COEF=DFLOAT(ME2+IS+1)                                             ACXY0723
      SAM=SAM+COEF*XFP(IT)*(A(IS,IS2)*XF(IT-1))                         ACXY0724
      GO TO 61                                                          ACXY0725
   66 IT=NF                                                             ACXY0726
      IS=IT-1                                                           ACXY0727
      IS2=2*IS                                                          ACXY0728
      COEF=DFLOAT(ME2+IS+1)                                             ACXY0729
      SAM=SAM+COEF*XFP(IT)*(A(IS,IS2)*XF(IT-1))                         ACXY0730
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)                         ACXY0731
      IT=NF+1                                                           ACXY0732
      IS=IT-1                                                           ACXY0733
      IS2=2*IS                                                          ACXY0734
      COEF=DFLOAT(ME2+IS+1)                                             ACXY0735
      SAM=SAM+COEF*XFP(IT)*(A(IS,IS2)*XF(IT-1))                         ACXY0736
   61 CONTINUE                                                          ACXY0737
C                                                                       ACXY0738
C        INTEGRAL B*                                                    ACXY0739
C                                                                       ACXY0740
      SIM=0.D0                                                          ACXY0741
      ALPHA=DFLOAT(IALPHA*(ME2+1)*(ME2+2))                              ACXY0742
      NMAX=MIN0(NFP,NF-2)                                               ACXY0743
      IF(NMAX.LT.1) GO TO 12                                            ACXY0744
      DO 10 IT=1,NMAX                                                   ACXY0745
      IS=IT-1                                                           ACXY0746
      IME=ME2+IS+IS+1                                                   ACXY0747
      COEFF=2.D0/DFLOAT(IME+2)                                          ACXY0748
      SIM=SIM+COEFF*XFP(IT)*ALPHA*(XF(IT+2)/DFLOAT(IME+4)               ACXY0749
     1-XF(IT)/DFLOAT(IME))                                              ACXY0750
      ALPHA=ALPHA*DFLOAT(ME2+IS+3)/DFLOAT(IT)                           ACXY0751
   10 CONTINUE                                                          ACXY0752
   12 IF(NFP-NF+1) 231,232,233                                          ACXY0753
  232 IT=NF-1                                                           ACXY0754
      IS=IT-1                                                           ACXY0755
      IS2=2*IS                                                          ACXY0756
      COEFF=2.D0/DFLOAT(ME2+IS2+3)                                      ACXY0757
      SIM=SIM-COEFF*XFP(IT)*ALPHA*XF(IT)/DFLOAT(ME2+IS2+1)              ACXY0758
      GO TO 231                                                         ACXY0759
  233 IT=NF-1                                                           ACXY0760
      IS=IT-1                                                           ACXY0761
      IS2=2*IS                                                          ACXY0762
      COEFF=2.D0/DFLOAT(ME2+IS2+3)                                      ACXY0763
      SIM=SIM-COEFF*XFP(IT)*ALPHA*XF(IT)/DFLOAT(ME2+IS2+1)              ACXY0764
      ALPHA=ALPHA*DFLOAT(ME2+IS+3)/DFLOAT(IS+1)                         ACXY0765
      IT=NF                                                             ACXY0766
      IS=IT-1                                                           ACXY0767
      IS2=2*IS                                                          ACXY0768
      COEFF=2.D0/DFLOAT(ME2+IS2+3)                                      ACXY0769
      SIM=SIM-COEFF*XFP(IT)*ALPHA*XF(IT)/DFLOAT(ME2+IS2+1)              ACXY0770
  231 CONTINUE                                                          ACXY0771
C                                                                       ACXY0772
C        INTEGRAL B"*                                                   ACXY0773
C                                                                       ACXY0774
      SJM=0.D0                                                          ACXY0775
      ALPHA=DFLOAT(IALPHA)                                              ACXY0776
      NMAX=MIN0(NFP,NF)                                                 ACXY0777
      DO 110 IT=1,NMAX                                                  ACXY0778
      IBMAX=(NFP-IT)/2+1                                                ACXY0779
      SP=0.D0                                                           ACXY0780
      DO 111 IB=1,IBMAX                                                 ACXY0781
      IN=2*IB-2+IT                                                      ACXY0782
  111 SP=SP+XFP(IN)                                                     ACXY0783
      SJM=SJM-SP*2.D0*ALPHA*XF(IT)                                      ACXY0784
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)                         ACXY0785
  110 CONTINUE                                                          ACXY0786
C                                                                       ACXY0787
C        INTEGRAL B'*                                                   ACXY0788
C                                                                       ACXY0789
      SKM=0.D0                                                          ACXY0790
      ALPHA=DFLOAT(IALPHA)                                              ACXY0791
      IF(NF.LT.4) XF(4)=0.D0                                            ACXY0792
      IF(NF.LT.3) XF(3)=0.D0                                            ACXY0793
      DO 121 IT=1,2                                                     ACXY0794
      IS=IT-1                                                           ACXY0795
      IS2=2*IS                                                          ACXY0796
      COEFF=-2.D0*ALPHA                                                 ACXY0797
      GUN=     (DFLOAT((ME2+IS)*(ME+IS+1)*IS)/DFLOAT(ME2+IS2-1)         ACXY0798
     1-DFLOAT((ME+IS)*(ME2+IS+1)*(IS+1))/DFLOAT(ME2+IS2+3))/DFLOAT(ME2  ACXY0799
     2+IS2+1)                                                           ACXY0800
      GDEU=DFLOAT((ME2+IS+2)*(ME+IS+3)*(ME2+IS+1))/DFLOAT((ME2+IS2+3)   ACXY0801
     1*(ME2+IS2+5))                                                     ACXY0802
      SP=0.D0                                                           ACXY0803
      IBMAX=(NFP-IT)/2+1                                                ACXY0804
      DO 122 IB=1,IBMAX                                                 ACXY0805
      IN=2*IB-2+IT                                                      ACXY0806
  122 SP=SP+XFP(IN)                                                     ACXY0807
      SKM=SKM+COEFF*SP*(XF(IT)*GUN+XF(IT+2)*GDEU)                       ACXY0808
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)                         ACXY0809
  121 CONTINUE                                                          ACXY0810
      NMAX=MIN0(NFP,NF-2)                                               ACXY0811
      IF(NMAX.LT.3) GO TO 125                                           ACXY0812
      DO 123 IT=3,NMAX                                                  ACXY0813
      IS=IT-1                                                           ACXY0814
      COEFF=-2.D0*ALPHA                                                 ACXY0815
      MIM=IS+ME                                                         ACXY0816
      IM=MIM+ME+1                                                       ACXY0817
      IME=IM+IS                                                         ACXY0818
      GUN=(DFLOAT((IM-1)*(MIM+1)*IS)/DFLOAT(IME-2)                      ACXY0819
     1-DFLOAT(MIM*IM*IT)/DFLOAT(IME+2))/DFLOAT(IME)                     ACXY0820
      GDEU=DFLOAT((IM+1)*(MIM+3)*IM)/DFLOAT((IME+2)*(IME+4))            ACXY0821
      GTROI=-DFLOAT((MIM-2)*(IS-1)*IS)/DFLOAT((IME-2)*(IME-4))          ACXY0822
      SP=0.D0                                                           ACXY0823
      IBMAX=(NFP-IT)/2+1                                                ACXY0824
      DO 124 IB=1,IBMAX                                                 ACXY0825
      IN=2*IB-2+IT                                                      ACXY0826
  124 SP=SP+XFP(IN)                                                     ACXY0827
      SKM=SKM+COEFF*SP*(XF(IT)*GUN+XF(IT+2)*GDEU+XF(IT-2)*GTROI)        ACXY0828
      ALPHA=ALPHA*DFLOAT(IM)/DFLOAT(IT)                                 ACXY0829
  123 CONTINUE                                                          ACXY0830
  125 IF(NMAX.EQ.NFP) GO TO 130                                         ACXY0831
      NMAX=NMAX+1                                                       ACXY0832
      IF(NMAX.LT.3) NMAX=3                                              ACXY0833
      NPMAX=MIN0(NFP,NF)                                                ACXY0834
      IF(NPMAX.LT.NMAX) GO TO 132                                       ACXY0835
      DO 131 IT=NMAX,NPMAX                                              ACXY0836
      IS=IT-1                                                           ACXY0837
      IS2=2*IS                                                          ACXY0838
      COEFF=-2.D0*ALPHA                                                 ACXY0839
      GUN=     (DFLOAT((ME2+IS)*(ME+IS+1)*IS)/DFLOAT(ME2+IS2-1)         ACXY0840
     1-DFLOAT((ME+IS)*(ME2+IS+1)*(IS+1))/DFLOAT(ME2+IS2+3))/DFLOAT(ME2  ACXY0841
     2+IS2+1)                                                           ACXY0842
      GTROI=-DFLOAT((ME+IS-2)*(IS-1)*IS)/DFLOAT((ME2+IS2-1)*(ME2+IS2-3))ACXY0843
      SP=0.D0                                                           ACXY0844
      IBMAX=(NFP-IT)/2+1                                                ACXY0845
      DO 134 IB=1,IBMAX                                                 ACXY0846
      IN=2*IB-2+IT                                                      ACXY0847
  134 SP=SP+XFP(IN)                                                     ACXY0848
      SKM=SKM+COEFF*SP*(XF(IT)*GUN+XF(IT-2)*GTROI)                      ACXY0849
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)                         ACXY0850
  131 CONTINUE                                                          ACXY0851
  132 IF(NPMAX.EQ.NFP) GO TO 130                                        ACXY0852
      NMAX=MIN0(NFP,NF+2)                                               ACXY0853
      NPMAX=NPMAX+1                                                     ACXY0854
      IF(NMAX.LT.NPMAX) GO TO 130                                       ACXY0855
      DO 141  IT=NPMAX,NMAX                                             ACXY0856
      IS=IT-1                                                           ACXY0857
      IS2=2*IS                                                          ACXY0858
      COEFF=-2.D0*ALPHA                                                 ACXY0859
      GTROI=-DFLOAT((ME+IS-2)*(IS-1)*IS)/DFLOAT((ME2+IS2-1)*(ME2+IS2-3))ACXY0860
      SP=0.D0                                                           ACXY0861
      IBMAX=(NFP-IT)/2+1                                                ACXY0862
      DO 144 IB=1,IBMAX                                                 ACXY0863
      IN=2*IB-2+IT                                                      ACXY0864
  144 SP=SP+XFP(IN)                                                     ACXY0865
      SKM=SKM+COEFF*SP*XF(IT-2)*GTROI                                   ACXY0866
      ALPHA=ALPHA*DFLOAT(ME2+IS+1)/DFLOAT(IS+1)                         ACXY0867
  141 CONTINUE                                                          ACXY0868
  130 CONTINUE                                                          ACXY0869
C                                                                       ACXY0870
C        INTEGRALS OVER LAMBDA                                          ACXY0871
C                                                                       ACXY0872
      SIGPE=PE+PEP                                                      ACXY0873
      SG=0.D0                                                           ACXY0874
      SPG=0.D0                                                          ACXY0875
      SPI=0.D0                                                          ACXY0876
      SPJ=0.D0                                                          ACXY0877
      SPK=0.D0                                                          ACXY0878
      SPN=0.D0                                                          ACXY0879
      SPK2=0.D0

      DO 70 N=1,NXI                                                     ACXY0880
      ALAM=XI(N)/SIGPE+1.D0                                             ACXY0881
      ALAM1=ALAM-1.D0                                                   ACXY0882
      ALAM2=ALAM+1.D0                                                   ACXY0883
      ALAMXI=ALAM1/ALAM2                                                ACXY0884
      FONC=DEXP(-SIGPE*ALAM+(SIGMA+SIGMAP)*DLOG(ALAM2))                 ACXY0885
      TRUC=ALAM1*ALAM2                                                  ACXY0886
      TRIC=TRUC**ME                                                     ACXY0887
      TRUC=TRIC*TRUC                                                    ACXY0888
      FISC=FONC*TRIC*OMEGA(N)                                           ACXY0889
      FONC=FONC*TRUC*OMEGA(N)                                           ACXY0890
      YP=XGP(NGP)                                                       ACXY0891
      DO 71 JN=2,NGP                                                    ACXY0892
      NN=NGP-JN+1                                                       ACXY0893
   71 YP=YP*ALAMXI+XGP(NN)                                              ACXY0894
      Y=XG(NG)                                                          ACXY0895
      DY=0.D0                                                           ACXY0896
      DO 72 JN=2,NG                                                     ACXY0897
      NN=NG-JN+1                                                        ACXY0898
      Y=Y*ALAMXI+XG(NN)                                                 ACXY0899
      DY=DY*ALAMXI+XG(NN+1)*DFLOAT(NN)                                  ACXY0900
   72 CONTINUE                                                          ACXY0901
      DY=2.D0*DY                                                        ACXY0902
      DY=DY+Y*(ME*ALAM/ALAMXI+ALAM2*(SIGMA-PE*ALAM2))                   ACXY0903
      DY=DY/(ALAM2*ALAM2)                                               ACXY0904
      FISC=FISC*Y*YP                                                    ACXY0905
      SG=SG+FONC*DY*YP                                                  ACXY0906
      SPG=SPG+FONC*Y*YP*ALAM                                            ACXY0907
      SPI=SPI+FISC*ALAM                                                 ACXY0908
      SPJ=SPJ+FONC*DY*YP*ALAM                                           ACXY0909
      SPK=SPK+FONC*Y*YP                                                 ACXY0910
      SPK2=SPK2+FONC*Y*YP*ALAM*ALAM
      SPN=SPN+FISC                                                      ACXY0911
   70 CONTINUE                                                          ACXY0912
      SG=SG/SIGPE                                                       ACXY0913
      SPG=SPG/SIGPE                                                     ACXY0914
      SPI=SPI/SIGPE                                                     ACXY0915
      SPJ=SPJ/SIGPE                                                     ACXY0916
      SPK=SPK/SIGPE                                                     ACXY0917
      SPN=SPN/SIGPE
      SPK2=SPK2/SIGPE
      WRITE(111,*) SPK, SPK2
  677 FORMAT(4E20.9)
C                                                                       ACXY0919
      CROT=+R*R*R*(SG*SAM-SPG*SUM+ME*(SPI*SAM+SUMP*SPG))/16.D0          ACXY0920
      GAMROT=+R*R*R*(SPJ*SIM-SPK*SKM+ME*(SPK*SJM+SIM*SPN))/16.D0        ACXY0921
      CROT=CROT+OMEZUT*GAMROT                                           ACXY0922
      IF(ME.NE.0) RETURN                                                ACXY0923
      CROT=CROT*DSQRT(2.D0)                                             ACXY0924
      GAMROT=GAMROT*DSQRT(2.D0)                                         ACXY0925
      RETURN                                                            ACXY0926
      END                                                               ACXY0927
