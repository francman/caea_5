
c================================================================ 

c SUBROUTINE:  FFT
c JIM COOLEY'S SIMPLE FFT PROGRAM--USES DECIMATION IN 
c TIME ALGORITHM
c X IS AN N=2**M POINT COMPLEX ARRAY THAT   INITIALLY 
c CONTAINS THE INPUT
c AND ON OUTPUT CONTAINS THE TRANSFORM
c THE PARAMETER INV   SPECIFIED DIRECT TRANSFORM IF 0 
c AND INVERSE IF 1
c----------------------------------------------------------------
c
      SUBROUTINE FFT(X, N, INV)
      COMPLEX X(1), U, W, T, CMPLX
c
c   X = COMPLEX ARRAY OF SIZE N--ON INPUT X CONTAINS
c       THE SEQUENCE TO BE TRANSFORMED
c       ON OUTPUT X CONTAINS THE DFT OF THE INPUT
c   N = SIZE OF FFT TO BE COMPUTED--N=2**M FOR  
c   1.LE.M.LE.15
c INV = PARAMETER TO DETERMINE WHETHER TO DO A DIRECT 
c TRANSFORM (INV=0)
c OR AN INVERSE TRANSFORM (INV=1)
c  
      M = ALOG(FLOAT(N))/ALOG(2.) + .1
      NV2 = N/2
      NM1 = N - 1
      J = 1
      DO 40 I=1,NM1
        IF (I.GE.J) GO TO 10
        T = X(J)
        X(J) = X(I)
        X(I) = T
  10    K = NV2
  20    IF (K.GE.J) GO TO 30
        J = J - K
        K = K/2
        GO TO 20
  30    J = J + K
  40  CONTINUE
      PI = 4.*ATAN(1.0)
      DO 70 L=1,M
        LE = 2**L
        LE1 = LE/2
        U = (1.0,0.0)
        W = CMPLX(COS(PI/FLOAT(LE1)),-SIN(PI/FLOAT(LE1)))
        IF (INV.NE.0) W = CONJG(W)
        DO 60 J=1,LE1
          DO 50 I=J,N,LE
            IP = I + LE1
            T = X(IP)*U
            X(IP) = X(I) - T
            X(I) = X(I) + T
  50      CONTINUE
          U = U*W
  60    CONTINUE
  70  CONTINUE
      IF (INV.EQ.0) RETURN
      DO 80 I=1,N
        X(I) = X(I)/CMPLX(FLOAT(N),0.)
  80  CONTINUE
      RETURN
      END
                                                        
