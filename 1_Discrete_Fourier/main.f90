! ****************************************************
!   FRANK MANU
!   SPRING 2020
!	REF: DR. THOMPSON
!   EECE.5200 - COMPUTER AIDED ENGINEERING ANAYLSIS
!   PROBLEM SET 5 - PART 1
! ****************************************************

!PART 1.A AND B
	IMPLICIT NONE
!==================
	INTEGER,PARAMETER:: N=16,N2=2*N
	COMPLEX,DIMENSION(0:N-1):: U,DUDX,UHALF
	REAL,   DIMENSION(0:N-1):: X
	
	COMPLEX,DIMENSION(0:N-1):: UDFT,CTMP
	REAL:: TPI,DX
	INTEGER:: I,K

	TPI = 2*ACOS(-1.0)
	DX = TPI/N

	DO I=0,N-1
	 	X(I)= I*DX
        	U(I)   =   COS(3*X(I))
		DUDX(I)= -3*SIN(3*X(I))
		UHALF(I)= COS(3* (X(I)+DX/2.0) )
	ENDDO

!FIND DFT OF U
	UDFT=U
	CALL FFT(UDFT,N,0)

! COMPUTE DERIVATIVE USING DFT	
	DO K=0,N/2
		CTMP(K) = UDFT(K) *CMPLX(0.0,K )
	ENDDO

	DO K=N/2+1,N-1
		CTMP(K) = UDFT(K)*CMPLX(0.0,N-K)*(-1)
	ENDDO

	CALL FFT(CTMP,N,1)

	WRITE(*,*) "Evaluate derivative for each x"
	WRITE(*,*) "   X              DERIVATIVE         APPROX             ERROR"

	DO I=0,N-1
		WRITE(*,*) X(I),REAL(DUDX(I)),REAL(CTMP(I)),REAL(DUDX(I)-CTMP(I))
	ENDDO

	WRITE(*,*)" "

!PART B
!-----------------------------
! COMPUTE  1/2 STEP USING DFT	
	DO K=0,N/2
		CTMP(K) = UDFT(K) *EXP(  K*CMPLX(0.0,DX/2.) )
	ENDDO

	DO K=N/2+1,N-1
		CTMP(K) = UDFT(K)*EXP((N-K)*CMPLX(0.0,DX/2.)*(-1) )
	ENDDO

	CALL FFT(CTMP,N,1)
	
	WRITE(*,*) "Evaluate derivate for each x/2 + x"
	WRITE(*,*) "   X              DERIVATIVE         APPROX             ERROR"

	DO I=0,N-1
	 WRITE(*,*) X(I)+DX/2.0,REAL(UHALF(I)),REAL(CTMP(I)),REAL(UHALF(I)-CTMP(I))
	ENDDO

	END
