PROGRAM main

USE smoother
USE hipfort
USE hipfort_check
USE ISO_C_BINDING

IMPLICIT NONE

  INTEGER, PARAMETER :: nW = 2
  INTEGER :: nX, nY, nIter
  REAL(prec), POINTER :: f(:,:)
  REAL(prec), POINTER :: smoothF(:,:)
  REAL(prec), POINTER :: weights(:,:)
  TYPE(c_ptr) :: f_dev
  TYPE(c_ptr) :: smoothF_dev
  TYPE(c_ptr) :: weights_dev
  REAL(prec) :: dx, dy, x, y
  INTEGER :: i, j, iter

    CALL GetCLIConf( nX, nY, nIter )

    ALLOCATE( f(1:nX,1:nY), smoothF(1:nX,1:nY), weights(-nW:nW,-nW:nW) )

    ! Allocate device memory
    CALL hipCheck( hipMalloc( f_dev, SIZEOF(f) ) )
    CALL hipCheck( hipMalloc( smoothF_dev, SIZEOF(smoothF) ) )
    CALL hipCheck( hipMalloc( weights_dev, SIZEOF(weights) ) )

    ! Create gaussian weights
    DO j = -nW, nW
      y = REAL(j,prec)
      DO i = -nW, nW
        x = REAL(i,prec)
        weights(i,j) = exp( -(x**2 + y**2)/(0.5_prec) )
      ENDDO
    ENDDO

    dx = 2.0_prec/REAL(nX, prec)
    dy = 2.0_prec/REAL(nY, prec)

    ! Initialize f and smoothF
    DO j = 1, nY
      y = -1.0_prec + dy*REAL(j-1,prec)
      DO i = 1, nX
        x = -1.0_prec + dx*REAL(i-1,prec)
        f(i,j) = tanh(x/0.001_prec)*tanh(y/0.001_prec)
        smoothF(i,j) = f(i,j)
      ENDDO
    ENDDO

    ! Write the initial condition to file
    OPEN(UNIT=2, FILE='function.txt', STATUS='REPLACE', ACTION='WRITE')
    DO j = 1, nY
      DO i = 1, nX
        WRITE(2,'(E9.4)') f(i,j)
      ENDDO
    ENDDO
    CLOSE(UNIT=2)

    ! Copy weights to weights_dev
    CALL hipCheck(hipMemcpy(weights_dev, c_loc(weights), SIZEOF(weights), hipMemcpyHostToDevice))

    DO iter = 1, nIter

      ! Copy f to f_dev
      CALL hipCheck(hipMemcpy(f_dev, c_loc(f), SIZEOF(f), hipMemcpyHostToDevice))

      smoothF = ApplySmoother( f, weights, nW, nX, nY )

      ! Copy smoothF_dev to smoothF
      CALL hipCheck(hipMemcpy(c_loc(smoothF), smoothF_dev, SIZEOF(smoothF), hipMemcpyDeviceToHost))

      CALL ResetF( f, smoothF, nW, nX, nY )

    ENDDO

    ! Write the initial condition to file
    OPEN(UNIT=2, FILE='smooth-function.txt', STATUS='REPLACE', ACTION='WRITE')
    DO j = 1, nY
      DO i = 1, nX
        WRITE(2,'(E9.4)') f(i,j)
      ENDDO
    ENDDO
    CLOSE(UNIT=2)

    ! Deallocate GPU memory
    CALL hipCheck( hipFree(f_dev) )
    CALL hipCheck( hipFree(smoothF_dev) )
    CALL hipCheck( hipFree(weights_dev) )

    DEALLOCATE( f, smoothF, weights )


END PROGRAM main

