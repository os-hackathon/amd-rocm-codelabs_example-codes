MODULE smoother

USE ISO_FORTRAN_ENV

IMPLICIT NONE

#ifdef DOUBLE_PRECISION
  INTEGER, PARAMETER :: prec = real64
#else
  INTEGER, PARAMETER :: prec = real32
#endif


CONTAINS

FUNCTION str2int( aString ) RESULT( aNumber )
  IMPLICIT NONE
  CHARACTER(*) :: aString
  INTEGER :: aNumber
  ! Local
  INTEGER :: stat

    READ(aString,*,iostat=stat) aNumber
    IF( stat /= 0 )THEN
      PRINT*, 'Error : Invalid input : '//TRIM(aString)
      STOP
    ENDIF

END FUNCTION str2int

SUBROUTINE GetCLIConf( nX, nY, nIter )
  IMPLICIT NONE
  INTEGER, INTENT(out) :: nX, nY, nIter
  ! Local
  INTEGER :: nArg, argID
  CHARACTER(6) :: arg


    nArg = command_argument_count( )
    IF( nArg /= 3 )THEN
      PRINT*, 'Error : Incorrect number of arguments. nArg = ', nArg
      CALL CLIHelp()
      STOP
    ENDIF

    DO argID = 1, nArg

      CALL get_command_argument( argID, arg )

      IF( argID == 1 )THEN
        nX = str2int( arg )
      ELSE IF( argID == 2 )THEN
        nY = str2int( arg )
      ELSE IF( argID == 3 )THEN
        nIter = str2int( arg )
      ENDIF

    ENDDO

    PRINT*, 'Info : nX = ', nX
    PRINT*, 'Info : nY = ', nY
    PRINT*, 'Info : nIter = ', nIter


END SUBROUTINE GetCLIConf

SUBROUTINE CLIHelp( )
  IMPLICIT NONE

      PRINT*, ' ------------------------------------------------------------------------------- '
      PRINT*, ''
      PRINT*, 'smoother '      
      PRINT*, ' An example 2d smoother application.'      
      PRINT*, ''
      PRINT*, ' Usage :'
      PRINT*, '    smoother nX nY nIter'      
      PRINT*, ''
      PRINT*, ' Input :'
      PRINT*, ''
      PRINT*, '    nX  - The number of grid cells in the x-direction'
      PRINT*, ''
      PRINT*, '    nY  - The number of grid cells in the y-direction'
      PRINT*, ''
      PRINT*, '    nIter  - The number of iterations to apply the smoothing operator'
      PRINT*, ''
      PRINT*, ' ------------------------------------------------------------------------------- '
END SUBROUTINE CLIHelp

FUNCTION ApplySmoother( f, weights, nW, nX, nY ) RESULT( smoothF )
  IMPLICIT NONE
  REAL(prec) :: f(1:nX,1:nY)
  REAL(preC) :: weights(-nW:nW,-nW:nW)
  INTEGER :: nW, nX, nY
  REAL(prec) :: smoothF(1:nX,1:nY)
  ! Local
  INTEGER :: i, j, ii, jj


    DO j = 1+nW, nY-nW
      DO i = 1+nW, nX-nW

        ! Take the weighted sum of f to compute the smoothF field
        smoothF(i,j) = 0.0_prec
        DO jj = -nW, nW
          DO ii = -nW, nW

            smoothF(i,j) = smoothF(i,j) + f(i+ii,j+jj)*weights(ii,jj)

          ENDDO
        ENDDO

      ENDDO
    ENDDO

END FUNCTION ApplySmoother

SUBROUTINE ResetF( f, smoothF, nW, nX, nY )
  IMPLICIT NONE
  REAL(prec), INTENT(inout) :: f(1:nX,1:nY)
  REAL(prec), INTENT(in) :: smoothF(1:nX,1:nY)
  INTEGER, INTENT(in) :: nW, nX, nY
  ! Local
  INTEGER :: i, j

    DO j = 1+nW, nY-nW 
      DO i = 1+nW, nX-nW 
        f(i,j) = smoothF(i,j)
      ENDDO
    ENDDO

END SUBROUTINE ResetF

END MODULE smoother
