PROGRAM fit
	USE analyze_timeseries
	USE analysis_functions
	USE make_model
	USE fit_library
	USE marquardt_library
 
	INTERFACE
		SUBROUTINE modelWrapper(X,A,Y,DYDA,NA)
			USE fit_library
			USE make_model
 
			REAL, DIMENSION(NA) :: A
			REAL, DIMENSION(25) :: X, Y
			REAL, DIMENSION(2) :: fRange, qyRange, zRange
			REAL, DIMENSION(25, NA) :: DYDA
			REAL :: fStep, qyStep, zStep
			INTEGER :: numRows, NA
		END SUBROUTINE modelWrapper
	END INTERFACE
 
	REAL, DIMENSION(:,:), POINTER :: model
 
	REAL, DIMENSION(:,:), POINTER :: dataFFT
	INTEGER :: numRowsFFT, numRowsModel
	REAL :: percentDiffTol, percentOutTol, sampleRate
	REAL :: alpha, AR, elongation, nu, theta, v, average, stdDev
	CHARACTER(LEN=64) :: inputFile, fftOutputFile, modelOutputFile
 
	INTEGER :: MA, MFIT, NCA, NDATA
	REAL :: CHISQ, ALAMDA
	REAL, DIMENSION(:), ALLOCATABLE :: A, SIG, X, Y
	REAL, DIMENSION(:,:), ALLOCATABLE :: AMATRIX, COVAR
	INTEGER, DIMENSION(:), ALLOCATABLE :: LISTA
 
	CALL getTimeseriesArguments(inputFile, sampleRate, percentDiffTol, percentOutTol)
	CALL getModelArguments(alpha, AR, elongation, nu, theta, v)
 
	dataFFT => analyzeTimeseries(inputFile, sampleRate, percentDiffTol, percentOutTol, numRowsFFT)
	WRITE(fftOutputFile, '(A4,A12,A5)') 'fft_', inputFile, '.data'
	CALL writeData(fftOutputFile, dataFFT, numRowsFFT)
 
	MA = 6
	MFIT = 2
	NCA = MA
	NDATA = 25
	CHISQ = 0.0
	ALAMDA = -0.1
 
	ALLOCATE(A(MA))
	ALLOCATE(AMATRIX(NCA,NCA))
	ALLOCATE(COVAR(NCA,NCA))
	ALLOCATE(LISTA(MA))
	ALLOCATE(SIG(NDATA))
	ALLOCATE(X(NDATA))
	ALLOCATE(Y(NDATA))
 
	A = (/alpha, AR, elongation, nu, theta, v/)
	LISTA = (/1,6/)
 
	X = dataFFT(1:9,1)
	Y = dataFFT(1:9,2)
 
	SIG = 0.1 * Y
 
	DO i=1,3
		CALL MRQMIN(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,COVAR,AMATRIX,NCA,CHISQ,modelWrapper,ALAMDA)
		PRINT*, "CHI, LAMBDA: ", CHISQ, ALAMDA
		PRINT*, "A: ", A
	ENDDO
 
	model => makeModel(A(1),A(2),A(3),A(4),A(5),A(6),numRowsModel)
	WRITE(modelOutputFile, '(A6,A12,A5)') 'model_', inputFile, '.data'		
	CALL writeData(modelOutputFile, model, numRowsModel)
END PROGRAM fit
 
SUBROUTINE modelWrapper(X,A,Y,DYDA,NA)
	USE marquardt_library
	USE make_model
	USE integrate_library
 
	INTERFACE
		FUNCTION integrandAlphaD(parameters, f, qy, z)
			REAL, DIMENSION(6) :: parameters
			REAL :: alpha, AR, elongation, nu, theta, v, f, qy, z, z0
			REAL :: au, PI, c, lambda, Vx, p, beta
			REAL(8) :: Fd, Fs, Ft, R, qx, qSquared
			REAL(8) :: integrandAlphaD
		END FUNCTION integrandAlphaD
 
		FUNCTION integrandARD(parameters, f, qy, z)
			REAL, DIMENSION(6) :: parameters
			REAL :: alpha, AR, elongation, nu, theta, v, f, qy, z, z0
			REAL :: au, PI, c, lambda, Vx, p, beta
			REAL(8) :: Fd, Fs, Ft, R, qx, qSquared
			REAL(8) :: integrandARD
		END FUNCTION integrandARD
 
		FUNCTION integrandThetaD(parameters, f, qy, z)
			REAL, DIMENSION(6) :: parameters
			REAL :: alpha, AR, elongation, nu, theta, v, f, qy, z, z0
			REAL :: au, PI, c, lambda, Vx, p, beta
			REAL(8) :: Fd, Fs, Ft, R, qx, qSquared
			REAL(8) :: integrandThetaD
		END FUNCTION integrandThetaD
 
		FUNCTION integrandVD(parameters, f, qy, z)
			REAL, DIMENSION(6) :: parameters
			REAL :: alpha, AR, elongation, nu, theta, v, f, qy, z, z0
			REAL :: au, PI, c, lambda, Vx, p, beta
			REAL(8) :: Fd, Fs, Ft, R, qx, qSquared
			REAL(8) :: integrandVD
		END FUNCTION integrandVD
	END INTERFACE
 
	REAL, DIMENSION(NA) :: A
	REAL, DIMENSION(25) :: X, Y
	REAL, DIMENSION(2) :: fRange, qyRange, zRange
	REAL, DIMENSION(25, NA) :: DYDA
	REAL :: fStep, qyStep, zStep
	INTEGER :: numRows, NA
 
	fRange = (/ 0.09765625,2.44140625 /)
	qyRange = (/ 1.0E-6,500.0E-6 /)
	zRange = (/ -2.0,0.0 /)
 
	fStep = 0.09765625
	qyStep = 10.0E-6
	zStep = 0.05
 
	numRows = getMaxIndex(fRange,fStep)
 
	X = getFrequencyList(fRange, fStep)
	Y = integrate(integrand, A, fRange, fStep, qyRange, qyStep, zRange, zStep)
 
	DYDA(:,1) = integrate(integrandAlphaD, A, fRange, fStep, qyRange, qyStep, zRange, zStep)
	DYDA(:,2) = integrate(integrandARD, A, fRange, fStep, qyRange, qyStep, zRange, zStep)
	DYDA(:,3) = 0
	DYDA(:,4) = 0
	DYDA(:,5) = integrate(integrandThetaD, A, fRange, fStep, qyRange, qyStep, zRange, zStep)
	DYDA(:,6) = integrate(integrandVD, A, fRange, fStep, qyRange, qyStep, zRange, zStep)
END SUBROUTINE modelWrapper
 
FUNCTION integrand(parameters, f, qy, z)
!Elongation should be input in degrees
!Theta should be input in arc seconds
!Computations done in SI units
!z should be input as multipliers of an AU (ex: z = (actual distance)/AU)
	IMPLICIT NONE
 
	REAL, DIMENSION(6) :: parameters
	REAL :: alpha, AR, elongation, nu, theta, v, f, qy, z, z0
	REAL :: au, PI, c, lambda, Vx, p, beta
	REAL(8) :: Fd, Fs, Ft, R, qx, qSquared
	REAL(8) :: integrand
 
	!REDACTED
 
	RETURN
END FUNCTION integrand
 
FUNCTION integrandAlphaD(parameters, f, qy, z)
!Elongation should be input in degrees
!Theta should be input in arc seconds
!Computations done in SI units
!z should be input as multipliers of an AU (ex: z = (actual distance)/AU)
	IMPLICIT NONE
 
	REAL, DIMENSION(6) :: parameters
	REAL :: alpha, AR, elongation, nu, theta, v, f, qy, z, z0
	REAL :: au, PI, c, lambda, Vx, p, beta
	REAL(8) :: Fd, Fs, Ft, R, qx, qSquared
	REAL(8) :: integrandAlphaD
 
	!REDACTED
 
	RETURN
END FUNCTION integrandAlphaD
 
FUNCTION integrandARD(parameters, f, qy, z)
!Elongation should be input in degrees
!Theta should be input in arc seconds
!Computations done in SI units
!z should be input as multipliers of an AU (ex: z = (actual distance)/AU)
	IMPLICIT NONE
 
	REAL, DIMENSION(6) :: parameters
	REAL :: alpha, AR, elongation, nu, theta, v, f, qy, z, z0
	REAL :: au, PI, c, lambda, Vx, p, beta
	REAL(8) :: Fd, Fs, Ft, R, qx, qSquared
	REAL(8) :: integrandARD
 
	!REDACTED
 
	RETURN
END FUNCTION integrandARD
 
FUNCTION integrandThetaD(parameters, f, qy, z)
!Elongation should be input in degrees
!Theta should be input in arc seconds
!Computations done in SI units
!z should be input as multipliers of an AU (ex: z = (actual distance)/AU)
	IMPLICIT NONE
 
	REAL, DIMENSION(6) :: parameters
	REAL :: alpha, AR, elongation, nu, theta, v, f, qy, z, z0
	REAL :: au, PI, c, lambda, Vx, p, beta
	REAL(8) :: Fd, Fs, Ft, R, qx, qSquared
	REAL(8) :: integrandThetaD
 
	!REDACTED
 
	RETURN
END FUNCTION integrandThetaD
 
FUNCTION integrandVD(parameters, f, qy, z)
!Elongation should be input in degrees
!Theta should be input in arc seconds
!Computations done in SI units
!z should be input as multipliers of an AU (ex: z = (actual distance)/AU)
	IMPLICIT NONE
 
	REAL, DIMENSION(6) :: parameters
	REAL :: alpha, AR, elongation, nu, theta, v, f, qy, z, z0
	REAL :: au, PI, c, lambda, Vx, p, beta
	REAL(8) :: Fd, dvFd, Fs, dvFs, Ft, dvFt, R, qx, qSquared
	REAL(8) :: integrandVD
 
	!REDACTED
 
	RETURN
END FUNCTION integrandVD
