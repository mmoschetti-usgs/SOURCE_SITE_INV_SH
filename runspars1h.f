      PROGRAM RUNSPARS1H
C
C  THIS PROGRAM RUNS SPARS1 FOR A REAL INPUT MATRIX.
C  THIS VERSION OF THE PROGRAM RUNS BOTH HORIZONTAL COMPONENTS.
C
C   This program calls routines to do an iterative inversion of a
C  large matrix.  This version has the matrix stored in the vector
C  Q in an arbitrary order.  The vectors II and JJ keep track of the
C  indices.  This program should be linked with "spars1" routines
C  to have consistent indexing.
C
C  (also calculates elapsed CPU time)
C
C  dimensions:
C		nc   = maximum # of columns in A matrix
C		nr   = maximum # of rows in A matrix
C		nit  = maximum # of iterations
C		nbig = maximum # of nonzero elements of A 
C-----------------------------------------------------------------------
C
      parameter ( nc=3555, nr=34049, nit=512, nbig=75000 )
c
      real q(nbig), x(nc), dx(nr), sigma(nit), w(nit), res(nr), d(nr)
      integer ii(nbig), jj(nbig)
      PI = 3.14159265
      nc1 = nc
      nr1 = nr
      na = nbig
      handle = 0
      icnts = 0
C
      WRITE(6,600)
 600  FORMAT(1X,'ENTER THE NUMBER OF ITERATIONS ON EIGENVALUE')
      READ(5,*) NITEIG
      WRITE(6,601)
 601  FORMAT(1X,'ENTER THE NUMBER OF ITERATIONS ON INVERSE')
      READ(5,*) NITINV
      WRITE(6,603)
 603  FORMAT(1X,'ENTER THE NUMBER OF REPETITIVE CALLS TO LSEIG')
      READ(5,*) NREPS
      WRITE(6,602)
 602  FORMAT(1X,'ENTER RATIO OF LARGEST TO SMALLEST EIGENVALUE')
      READ(5,*) EIGRAT
      WRITE(6,607)
 607  FORMAT(1X,'DO YOU WANT TO READ IN A FILE OF SIGNAL TO',/,
     +1X,'NOISE RATIOS TO BE USED TO WEIGHT THE ROWS OF THE',/,
     +1X,'MATRIX  1=YES 0=NO')
      READ(5,*) IWT 
C
C   READ IN II, AND JJ VECTORS HERE
C
      OPEN(UNIT=15,FILE='II.DAT',STATUS='OLD',ACCESS='SEQUENTIAL',
     +FORM='UNFORMATTED')
      READ(15) NRA,NCA,NSOR,NSIT,ICREC,IFN,NII,NIIOLD
      READ(15) (II(I),I=1,NII)
      CLOSE(15)
      OPEN(UNIT=15,FILE='JJ.DAT',STATUS='OLD',ACCESS='SEQUENTIAL',
     +FORM='UNFORMATTED')
      READ(15) NRA,NCA,NSOR,NSIT,ICREC,IFN,NJJ,NJJOLD
      READ(15) (JJ(I),I=1,NJJ)
      CLOSE(15)
      NQ=NJJ
C
C READ IN DATA VECTOR, D, AND INITIALIZE SOLUTION VECTOR, X. 
C
      IF(IWT .EQ. 1) THEN
      OPEN(UNIT=14,FILE='WTrad.DAT',STATUS='OLD',ACCESS='SEQUENTIAL',
     +FORM='UNFORMATTED')
      OPEN(UNIT=16,FILE='WTtrans.DAT',STATUS='OLD',ACCESS='SEQUENTIAL',
     +FORM='UNFORMATTED')
      END IF
      OPEN(UNIT=15,FILE='D.DAT',STATUS='OLD',ACCESS='SEQUENTIAL',
     +FORM='UNFORMATTED')
      DO 40 I=1,ICREC
      K1=IFN*(I-1)+1
      K2=IFN*I
      READ(15) ILOOP,JLOOP,IFN,FMIN,FSTEP
      READ(15) (D(K),K=K1,K2)
      IF(IWT .EQ. 1) THEN
      IUN=14
      ITT=I/2
      TT2=FLOAT(ITT)
      TT1=FLOAT(I)/2.
      IF(TT1 .EQ. TT2) IUN=16
      READ(IUN) ILOOP,JLOOP,IFN,FMIN,FSTEP
      READ(IUN) (DX(K),K=K1,K2)
      DO 32 K=K1,K2
      SNWT=AMIN1(DX(K),5.0)
      SNWT=AMAX1(SNWT,1.0)
      SNWT=SNWT/5.0
      DX(K)=SNWT
  32  D(K)=D(K)*SNWT
      END IF
  40  CONTINUE
      IF(IWT .EQ. 1) THEN
      CLOSE(14)
      CLOSE(16)
      END IF
      K1=IFN*ICREC+1
      K2=IFN*(ICREC+1)
      READ(15) WTSIT,IAVG,IFN,FAVG
      READ(15) (D(K),K=K1,K2)
      CLOSE(15)
      DO 45 I=1,NCA
  45  X(I)=0.
C
C     SET UP Q ARRAY
C
      K=1
      ICC=0
      DO 30 I=1,NJJOLD
      SNWT=1.0
      IF(IWT .EQ. 1) THEN
      ICC=ICC+1
      IF(ICC .GT. 2) THEN
      ICC=1
      K=K+1
      END IF
      SNWT=DX(K)
      END IF
  30  Q(I)=1.0*SNWT
      DO 31 I=1,NRA
  31  DX(I)=0.      
      I1=NJJOLD+1
      DO 50 I=I1,NQ
  50  Q(I)=WTSIT
C
C  the first call to lseig initializes the largest eigenvalue and
C  eigenvector.  ncyc and nsig should be set to zero in this call.
C
      ncyc = 0
      nsig = 0
      call lseig(q,ii,jj,na,nq,nca,nra,ncyc,nsig,x,dx,sigma,w,
     .           smax,err,sup,res)
      OPEN(UNIT=10,FILE='RUNSPARS1.OUT',STATUS='NEW')
      rewind 10
      WRITE(10,608) WTSIT,IAVG,IFN,FAVG
 608  FORMAT(1X,'WTSIT= ',F12.6,' IAVG= ',I5,' IFN= ',I5,' FAVG= ',
     +F12.6)
      nz = 0
      write(10,*) 'Calculating largest eigenvalue of ',nca,' by ',
     .             nra,' array'
C
C  this call does "niteig" iterations on the eigenvalue, eigenvector
C  estimate.
C
      DO 35 IREPS=1,NREPS
      call lseig(q,ii,jj,na,nq,nca,nra,niteig,nsig,x,dx,sigma,w,
     .             smax,err,sup,res)
      WRITE(10,604) SMAX,ERR,SUP
 604  FORMAT(1X,'SMAX= ',E14.6,' ERR= ',E14.6,' SUP= ',E14.6)
  35  CONTINUE
C
C calculate the weights for each iteration based on the largest and
C smallest eigenvalues over which uniform convergence is desired.
C
      slo = 0
      if ( eigrat.ne.0.0 ) slo = smax/eigrat
      print *, 'uniform convergence over the range: ',smax,' to ',slo
      call chebyu(sigma,nitinv,smax,slo,w)
      ERRLSQ=ERRLIM(SIGMA,NITINV,SMAX,SLO)
      WRITE(10,605) ERRLSQ
 605  FORMAT(1X,'ERRLSQ= ',E14.6)
      do 80 j=1,nca
 80     x(j) = 0.0
C
C calculate the new estimate of "x"
C
      call lsquc (q,ii,jj,na,nq,nca,nra,x,dx,d,nitinv,sigma,res)
      write(10,*) '---- iterative estimate of solution is ----'
      WRITE(10,606) (X(I),I=1,NCA)
 606  FORMAT(1X,6E14.6)
      OPEN(UNIT=12,FILE='SOLUTIONSP1.DAT',STATUS='NEW',
     +ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      WRITE(12) NCA
      WRITE(12) (X(I),I=1,NCA)
      CLOSE(12)
C
C FROM THE ESTIMATE OF X CALCULATE WHAT THE ESTIMATE OF THE DATA IS
C
      DO 110 I=1,NRA
 110    DX(I) = 0.0
      DO 120 IQ=1,NQ
120     DX(II(IQ)) = DX(II(IQ)) + Q(IQ)*X(JJ(IQ))
      OPEN(UNIT=9,FILE='DFIT.DAT',STATUS='NEW',FORM='UNFORMATTED',
     .ACCESS='SEQUENTIAL')
      rewind 9
      NT = NRA
      write(9) nt 
      write(9) (Dx(i),i=1,nt)
      close(9)
      close(10)
      stop
      end
c----------------------------------------------------------------------
      SUBROUTINE LSQUC(Q,II,JJ,NA,NQ,NC,NR,X,DX,D,NCYCLE,SIGMA,RES)
C----------------------------------------------------------------------
C
C       THIS VERSION MODIFIED ON 1/3/85 HANDLES A SPARSE A MATRIX WHICH
C   IS CONTAINED IN RANDOM ORDER IN THE VECTOR Q.  THE LOCATION OF EACH
C   ELEMENT IS CONTAINED IN II AND JJ.  THIS IS THE VERSION SPARS1 WHICH
C   CORRESPONDS TO APPROACH A IN PKP NOTES.  (GCB)
C
C     LEAST SQUARES SOLUTION USING RICHARDSON'S ALGORITHM
C     WITH CHEBYSHEV ACCELERATION.  THE STEP SIZE IS VARIED TO OBTAIN
C     UNIFORM CONVERGENCE OVER A PRESCRIBED RANGE OF EIGENVALUES.
C
C                             ..... ALLEN H. OLSON 6-29-85.
C----------------------------------------------------------------------
C
C                       NC                   ~
C          GIVEN :     SUM [ A(J,I) * X(J) ] = D(I)   ;  I=1, ... NR
C                      J=1
C                                               T
C          MINIMIZE :  || A*X - D || = (A*X - D) * (A*X -D)
C
C-------------------------------------------------------------------Y--
C -------
C  INPUT
C -------
C     Q(1...NQ)  ----  - VECTOR CONTAINS A MATRIX DEFINED ABOVE.
C                        NOTE  FOR THE A MATRIX THE FIRST INDEX IS THE
C                        COLUMN INDEX! I.E. THE MATRIX IS STORED IN ROW
C                        ORDER.
C     II(1...NQ) ----  - CONTAINS ROW LOCATION OF Q
C     JJ(1...NQ) ----  - CONTAINS COLUMN LOCATION OF Q 
C     RES(1...NR) ---  - STORAGE VECTOR FOR MULTIPLICATIONS
C     NA    ---------  - DIMENSION OF Q(.)
C     NQ    ---------  - DIMENSION OF Q ACTUALLY USED 
C     X(1...NC) -----  - INITIAL GUESS SOLUTION; CAN BE SET TO ZERO OR
C                        VALUES RETURNED FORM PREVIOUS CALLS TO LSQUC.
C     DX(1...NC) ----  - TEMPORARY STORAGE ARRAY
C     D(1...NR) -----  - DATA AS DEFINED ABOVE
C     NCYCLE  -------  - NUMBER OF ITERATIONS TO PERFORM.
C     SIGMA(1..NCYCLE) - ARRAY CONTAINING THE WEIGHTS FOR STEP SIZES.
C                        SEE SUBROUTINE 'CHEBYU' FOR COMPUTING THESE. 
C --------
C  OUTPUT
C --------
C     X(1..NC)   - SOLUTION VECTOR AS DEFINED ABOVE
C                  ONLY ARRAY X(..) IS OVER-WRITTEN BY LSQUC.
C
C----------------------------------------------------------------------
      REAL Q(NA),X(*),DX(*),D(*),SIGMA(*),RES(*)
      INTEGER II(NA),JJ(*)
C
      DO 3000 ICYC=1,NCYCLE
C
      DO 1000 J=1,NC
 1000 DX(J)=0.0
C
      DO 1010 I=1,NR
 1010 RES(I) = 0.0
C
      DO 100 N=1,NQ
        I = II(N)
        J = JJ(N)
 100    RES(I) = RES(I) + Q(N)*X(J)
C
      DO 101 I=1,NR
 101    RES(I) = D(I) - RES(I) 
C
      DO 200 N=1,NQ
        I = II(N)
        J = JJ(N)
 200    DX(J) = DX(J) + RES(I)*Q(N)
C
      DO 1300 J=1,NC
 1300 X(J)=X(J)+DX(J)/SIGMA(ICYC)
C
 3000 CONTINUE
C
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE CHEBYU(SIGMA,NCYCLE,SHI,SLO,WORK)
C----------------------------------------------------------------------
C
C COMPUTES THE CHEBYSHEV WEIGHTS WITH UNIFORM DISTRIBUTION.
C WEIGHTS ARE ORDERED PAIR-WISE IN SUCH A FASHION THAT AFTER AN EVEN
C NUMBER OF STEPS THEY ARE DISTRIBUTED UNIFORMLY ON THE INTERVAL [SLO,SHI].
C THIS ORDERING PROVIDES OPTIMUM NUMERICAL STABILITY OF ROUTINE LSQUC.
C
C                             ..... ALLEN H. OLSON 6-29-85.
C
C----------------------------------------------------------------------
C -------
C  INPUT
C -------
C
C     NCYCLE   ------  - MUST BE A POWER OF TWO!  NUMBER OF ITERATIONS.
C                        
C     SHI,SLO  ------  - HIGH AND LOW LIMITS DEFINING THE BAND OF EIGENVALUES
C                        TO RETAIN IN THE SOLUTION.  SHI >= LARGEST EIGENVALUE
C                        OF THE NORMAL EQUATIONS.
C     WORK(1..NCYCLE)  - WORK ARRAY FOR SORTING ARRAY SIGMA(..).
C
C -------
C  OUTPUT
C -------
C
C     SIGMA(1..NCYCLE) - WEIGHTS FOR THE STEP SIZES IN ROUTINE LSQUC.
C -------
C CALLS SUBROUTINE SPLITS.
C
C----------------------------------------------------------------------
      DIMENSION SIGMA(*),WORK(*)
C      
      PI=3.1415927
C      
C  SET UP THE CHEBYSHEV WEIGHTS IN INCREASING ORDER
C
      DO 100 I=1,NCYCLE
      SIGMA(I)=-COS( (2*I-1)*PI/2/NCYCLE )
      SIGMA(I)=( SIGMA(I)*(SHI-SLO)+(SHI+SLO) )/2
  100 CONTINUE
C
C  SORT THE WEIGHTS
C
      LEN=NCYCLE
  200 NSORT=NCYCLE/LEN
      DO 300 IS=1,NSORT
      I0=1+(IS-1)*LEN
      CALL SPLITS(SIGMA(I0),WORK,LEN)
  300 CONTINUE
      LEN=LEN/2
      IF(LEN.GT.2) GOTO 200
C
      RETURN
      END
C---------------
      SUBROUTINE SPLITS(X,T,N)
C---------  CALLED BY CHEBYU
      DIMENSION X(*),T(*)
      L=0
      DO 20 I=1,N,2
      L=L+1
   20 T(L)=X(I)
      DO 30 I=2,N,2
      L=L+1
   30 T(L)=X(I) 
C
      NB2=N/2
      NB2P1=NB2+1
      IF(NB2.GE.2) THEN
        DO 40 I=1,NB2
   40   X(I)=T(NB2-I+1)
        DO 50 I=NB2P1,N
   50   X(I)=T(I)
      ELSE
        DO 60 I=1,N
   60   X(I)=T(I)
      ENDIF 
      RETURN
      END
C----------------------------------------------------------------------
      FUNCTION ERRLIM(SIGMA,NCYCLE,SHI,SLO)
C  RETURNS LIMIT OF THE MAXIMUM THEORETICAL ERROR USING CHEBYSHEV WEIGHTS
C----------------------------------------------------------------------
      DIMENSION SIGMA(*)
      ERRLIM=1.0
      DELTA=0.25*(SHI-SLO)
      DO 10 I=1,NCYCLE
   10 ERRLIM=ERRLIM*DELTA/SIGMA(I)
      ERRLIM=2*ERRLIM
      RETURN
      END 
C----------------------------------------------------------------------
      REAL FUNCTION ERRVAL(X,SIGMA,NCYCLE)
C  COMPUTES THE THEORETICAL ERROR AT EIGENVALUE 'X'.
C----------------------------------------------------------------------
      DIMENSION SIGMA(*)
C
      ERRVAL=1.0
      DO 50 K=1,NCYCLE
   50 ERRVAL=ERRVAL*(1.0-X/SIGMA(K))
      ERRVAL=ABS(ERRVAL)
      RETURN
      END
C----------------------------------------------------------------------
      REAL FUNCTION ERRRAT(X1,X2,SIGMA,NCYCLE)
C  COMPUTES THE RATIO OF THE ERROR AT EIGENVALUE X1 TO THE ERROR AT X2.
C----------------------------------------------------------------------
      DIMENSION SIGMA(*)
C
      ERRRAT=1.0
      RAT=X1/X2
      DO 50 K=1,NCYCLE
   50 ERRRAT=ERRRAT*RAT*(1.0-SIGMA(K)/X1)/(1.0-SIGMA(K)/X2)
      ERRRAT=ABS(ERRRAT)
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE LSEIG(Q,II,JJ,NQ1,NQ,NC,NR,NCYC,NSIG,X,DX,SIGMA,W,
     .                 SMAX,ERR,SUP,RES)
C----------------------------------------------------------------------
C     LEAST-SQUARES EIGENVALUE
C
C     MODIFIED 1/3/85 TO HANDLE SPARSE RANDOMLY STORED A MATRIX HOUSED
C     IN THE Q VECTOR - "SPARS1"
C
C     ITERATIVELY ESTIMATES THE LARGEST EIGENVALUE AND EIGENVECTIOR
C     WITH ERROR BOUNDS FOR THE LEAST-SQUARES NORMAL MATRIX A'A.
C     A CHEBYSHEV CRITERION IS USED TO CALCULATE THE OPTIMUM SET OF
C     ORIGIN SHIFTS IN ORDER TO ACCELERATE CONVERGENCE.
C     BASED UPON THE RAYLEIGH QUOTIENT AND ERROR ANALYSIS PRESENTED IN
C     J. H. WILKINSON'S "THE ALGEBRAIC EIGENVALUE PROBLEM",
C     (PP 170 ...), (PP 572 ...).
C     UNDER VERY PESIMISTIC ASSUMPTIONS REGARDING THE STARTING VECTOR,
C     THE ALGORITHM WILL INITIALLY CONVERGE TO AN EIGENVALUE LESS THAN
C     THE LARGEST.  HENCE, A CORRESPONDING PESSIMISTIC ESTIMATE OF AN
C     UPPER-BOUND ON THE LARGEST EIGENVALUE IS ALSO MADE.
C
C                             ..... ALLEN H. OLSON 10-4-85.
C----------------------------------------------------------------------
C -------
C  INPUT
C -------
C     Q(1...NQ) -----  - VECTOR HOLDING MATRIX AS DEFINED IN SUBROUTINE
C                        LSQUC. DIMENSIONED AT LEAST (NC*NR)
C     II(1...NQ) ----  - KEEPS ROW INDEX OF A
C     JJ(1...NQ) ----  - KEEPS COLUMN INDEX OF A
C     RES(1...NR) ---  - UTILITY VECTOR FOR MULTIPLICATIONS
C     NA    ---------  - DIMENSION OF Q.
C     NQ    ---------  - DIMENSION OF Q THAT ACTUALLY GETS USED.
C     NC,NR  --------  - NUMBER OF COLUMNS, ROWS DEFINED IN ROUTINE LSQUC
C     NCYC  ---------  - NUMBER OF CHEBYSHEV ITERATIONS TO PERFORM. 
C                          MUST BE A POWER OF TWO.
C     NSIG  ---------  - CUMULATIVE NUMBER OF ITERATIONS PERFORMED BY
C                          PREVIOUS CALLS TO THIS ROUTINE.  MUST BE SET
C                          TO ZERO ON INITIAL CALL.  NSIG IS AUTOMATICALLY
C                          INCREMENTED BY THIS ROUTINE AND MUST NOT BE
C                          REDEFINED ON SUBSEQUENT CALLS BY CALLING PROGRAM.
C     X(1..NC)  -----  - INITAL GUESS FOR THE EIGENVECTOR.
C     DX(1...NC) ----  - TEMPORARY STORAGE ARRAY FOR X(.).
C     SIGMA(1..NSMX)-  - ARRAY FOR HOLDING THE CHEBYSHEV ORIGIN SHIFTS.
C                          EACH CALL TO LSEIG PERFORMS NCYC+1 ITERATIONS.
C                          NSMX MUST BE GREATER THAN OR EQUAL TO THE
C                          CUMULATIVE NUMBER OF ITERATIONS TO BE PERFORMED.
C     W(1..NSMX)  ---  - TEMPORARY STORAGE ARRAY FOR SIGMA(.).
C     SMAX  ---------  - INITIAL GUESS FOR THE EIGENVALUE.
C --------
C  OUTPUT
C --------
C     X(1...NC) -----  - REVISED ESTIMATE OF LARGEST EIGENVECTOR.
C     SMAX  ---------  - REIVSED ESTIMATE OF LARGEST EIGENVALUE OF A'A.
C     ERR   ---------  - ERROR BOUND FOR SMAX. WE ARE GUARANTEED THAT
C                          AT LEAST ONE EIGENVALUE IS CONTAINED IN THE
C                          INTERVAL SMAX+-ERR.  IN THE NEIGHBORHOOD OF
C                          CONVERGENCE, THIS WILL CONTAIN THE MAXIMUM.
C     SUP   ---------  - A PESSIMISTIC UPPER BOUND FOR THE LARGEST
C                          EIGENVECTOR.
C ------
C  NOTE:(1) FOR NCYC=0, AN INITIAL GUESS FOR THE EIGENVECTOR IS FORMED
C ------  BY SUMMING THE ROWS OF THE MATRIX SO THAT THE ACCUMULATED
C         VECTOR INCREASES IN LENGTH AS EACH ROW IS ADDED.  ONE
C         ITERATION OF THE POWER METHOD IS THEN PERFORMED TO ESTIMATE
C         SMAX.
C
C       (1.5) IN THIS VERSION THE INITIAL GUESS TO THE EIGENVECTOR IS
C         SIMPLY THE SUM OF THE ROWS OF THE A MATRIX --- UNLIKE THE
C         PREVIOUS VERSION, THEY ARE ALWAYS SUMMED (GCB 1/7/85).
C
C       (2) CAUTION: AN ANOMALOUS BAD GUESS FOR THE INITIAL EIGENVECTOR
C         BEING VIRTUALLY ORTHOGONAL TO THE LARGEST EIGENVECTOR
C         WILL CAUSE EARLIER ITERATIONS TO CONVERGE TO THE NEXT LARGEST
C         EIGENVECTOR.  THIS IS IMPOSSIBLE TO DETECT.  IN THIS CASE,
C         SMAX MAY BE MUCH LESS THAN THE LARGEST EIGENVALUE.
C         THE PARAMETER EPS SET BELOW IS A PESSIMISTIC ASSUMPTION ABOUT
C         THE RELATIVE SIZE OF THE COMPONENT OF THE LARGEST EIGENVECTOR
C         IN THE INITIAL ITERATION.  FROM THIS, AN UPPER-BOUND IS
C         CALCULATED FOR THE LARGEST EIGENVALUE, SUP.  SUP WILL ALWAYS
C         BE LARGER THAN SMAX AND REFLECTS THE UNCERTAINTY DUE TO AN
C         ANOMOLOUS BAD CHOICE FOR THE STARTING VECTOR.
C -------------------------
C  SAMPLE CALLING SEQUENCE
C -------------------------
C     NCYC=0
C     NSIG=0
C     CALL LSEIG(Q,II,JJ,NA,NQ,NC,NR,NCYC,NSIG,X,DX,SIGMA,W,SMAX,ERR,
C                SUP,RES)
C     NCYC=4
C     CALL LSEIG(Q,II,JJ,NA,NQ,NC,NR,NCYC,NSIG,X,DX,SIGMA,W,SMAX,ERR,
C                SUP,RES)
C     NCYC=8
C     CALL LSEIG(Q,II,JJ,NA,NQ,NC,NR,NCYC,NSIG,X,DX,SIGMA,W,SMAX,ERR,
C                SUP,RES)
C
C     THE FIRST CALL WITH NCYC=0 INITIALIZES X(.) AND SMAX.  IF THESE
C     ARE ALREADY KNOWN THEN NCYC CAN BE SET TO A NONZERO VALUE FOR
C     THE FIRST CALL.  NSIG MUST ALWAYS BE ZERO FOR THE FIRST CALL
C     HOWEVER.  THE NEXT TWO CALLS PERFORM CHEBYSHEV ITERATION TO
C     IMPROVE X(.) AND SMAX.  UPON COMPLETION, A TOTAL OF NSIG=1+5+9=15 
C     ITERATIONS HAVE ACTUALLY BEEN PERFORMED.  BY MAKING REPEATED
C     CALLS TO LSEIG, THE ERROR IN SMAX AND THE DIFFERENCE BETWEEN
C     SMAX AND SUP CAN BE MONITORED UNTIL THE DESIRED LEVEL OF
C     CERTAINTY IS OBTAINED.
C
C----------------------------------------------------------------------
C
      REAL Q(NQ1),X(*),DX(*),SIGMA(*),W(*),RES(*)
      INTEGER II(NQ1),JJ(*)
C
      EPS=1.E-6
      IF(NCYC.EQ.0) THEN
C
C -----------------  INTIALIZE X VECTOR WITH FIRST COLUMN OF A
C
        DO 53 I=1,NC
 53     X(I) = 0.0
C
        DO 54 J=1,NR
 54     RES(J) = 0.0
C
        DO 55 N=1,NQ
 55     X(JJ(N)) = Q(N)
C
        RES1=0.0
        DO 200 J=1,NC
 200    RES1=RES1+X(J)*X(J)
        RES1=1.0/SQRT(RES1)
C
        DO 300 J=1,NC
 300      X(J)=X(J)*RES1
C
       ELSE
        SLO=0.0
        CALL CHEBYU(SIGMA(NSIG+1),NCYC,SMAX,SLO,W)
      ENDIF
C
      NSIG1=NSIG+1
      NSIG =NSIG1+NCYC
      SIGMA(NSIG)=0.0
      DO 3000 ICYC=NSIG1,NSIG
C
C ------------------------   RE-ESTIMATE EIGENVECTOR
C
      DO 1000 J=1,NC
 1000 DX(J)=0.0
C
      DO 1010 I=1,NR
 1010 RES(I) = 0.0
C
      DO 1200 N=1,NQ
        I = II(N)
        J = JJ(N)
 1200   RES(I) = RES(I) + Q(N)*X(J)
C
      DO 1210 N=1,NQ
        I = II(N)
        J = JJ(N)
 1210   DX(J) = DX(J) + RES(I)*Q(N)
C
      DO 1250 J=1,NC
 1250 DX(J)=DX(J) - SIGMA(ICYC)*X(J)
C
      SMAX=0.0
      DO 1300 J=1,NC
 1300 SMAX=SMAX+ DX(J)*DX(J)
      SMAX=SQRT(SMAX)
C
      IF(ICYC.EQ.NSIG) THEN
      ERR=0.0
      DO 1350 J=1,NC
      RES1=DX(J)-SMAX*X(J)
 1350 ERR=ERR+ RES1*RES1
      ERR=SQRT(ERR)
      ENDIF
C
      DO 1400 J=1,NC
 1400 X(J)=DX(J)/SMAX
C
 3000 CONTINUE
C
      SLO=SMAX
      SUP=(1.+EPS)*SMAX/(EPS**(1.0/NSIG))
      DO 4000 ICYC=1,25
      SMP=(SUP+SLO)/2
      ERRSMP=ERRRAT(SMAX,SMP,SIGMA,NSIG)
      IF(ERRSMP.GT.EPS) THEN
        SLO=SMP
      ELSE
        SUP=SMP
      ENDIF
      RES1=(SUP-SLO)/SLO
      IF(RES1.LE.EPS) GO TO 4500
 4000 CONTINUE
 4500 CONTINUE
C
      RETURN
      END
