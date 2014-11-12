! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! Copyright (C) 2014 by Matthew Clay <mpclay@gmail.com>
!
!> @file FDUpwind.F90
!> @author Matthew Clay
!> @brief Solve the Burgers' equation using a non-conservative FD scheme.
!!
!! In this program, we use a non-conservative FD scheme to solve the Burgers'
!! equation. The purpose of this investigation is to see what happens when a
!! non-conservative scheme is used to solve a conservation-law problem.
PROGRAM FDUpwind_p

   ! Required modules.
   USE Parameters_m,ONLY: RWP, IWP, PI

   IMPLICIT NONE

   ! User inputs for the simulation.
   !
   !> Length of the domain.
   REAL(KIND=RWP),PARAMETER :: L = 1.0_RWP
   !> Number of grid points.
   INTEGER(KIND=IWP),PARAMETER :: n = 1024_IWP
   !> CFL for time integration.
   REAL(KIND=RWP),PARAMETER :: cfl = 0.1_RWP
   !> End time for the simulation.
   REAL(KIND=RWP),PARAMETER :: tend = 1.0_RWP
   !> Time frequency with which to write data files.
   REAL(KIND=RWP),PARAMETER :: writePeriod = 0.01_RWP
   !> Time frequency with which to print information to the user.
   REAL(KIND=RWP),PARAMETER :: printPeriod = 0.01_RWP

   ! Variables required for simulation.
   !
   !> The computational grid.
   REAL(KIND=RWP),DIMENSION(1:n) :: x
   !> The local grid spacing.
   REAL(KIND=RWP) :: dx
   !> Solution array. We include 1 left ghost layer for the spatial derivative.
   REAL(KIND=RWP),DIMENSION(0:n) :: u
   !> Array for the RHS of the ODE system.
   REAL(KIND=RWP),DIMENSION(1:n) :: Lu
   !> Time step.
   REAL(KIND=RWP) :: dt
   !> Current simulation time.
   REAL(KIND=RWP) :: t
   !> Current simulation step.
   INTEGER(KIND=IWP) :: nadv
   !> Time after which a data file will be written.
   REAL(KIND=RWP) :: writeTime = writePeriod
   !> Time after which information will be printed to the user.
   REAL(KIND=RWP) :: printTime = printPeriod

   ! Extraneous variables.
   !
   !> Looping index.
   INTEGER(KIND=IWP) :: i
   !> Variables to control when to stop looping.
   LOGICAL :: loopBool, exitThisStep
   !> Logic to determine when to write data and print info to the user.
   LOGICAL :: printBool, writeBool

   ! Zero out necessary variables.
   t = 0.0_RWP
   nadv = 0_IWP

   ! Fill in the computational grid.
   dx = L/REAL(n-1, RWP)
   DO i = 0, n-1
      x(i+1) = REAL(i, RWP)*dx
   END DO
   x(n) = L

   ! Form the initial conditions.
   DO i = 1, n
      u(i) = 1.0_RWP + SIN(2.0_RWP*PI*x(i)/L)
   END DO
   !
   ! Write out the initial conditions.
   CALL WriteData(n, x, u, nadv)

   ! Calculate the initial time step.
   CALL CalculateTimeStep(n, u, dx, cfl, dt)

   ! Print some information to the user.
   WRITE(*,100) '-------------------------------------------------------------'
   WRITE(*,100) 'FDUpwind: A non-conservative scheme to solve Burgers equation'
   WRITE(*,100) '-------------------------------------------------------------'
   WRITE(*,100) ''
   WRITE(*,120) 'Domain length:', L
   WRITE(*,130) 'Num. of points:', n
   WRITE(*,120) 'CFL number:', cfl
   WRITE(*,120) 'End time:', tend
   WRITE(*,100) ''
   100 FORMAT (A)
   120 FORMAT (A,T16,ES15.8)
   130 FORMAT (A,T17,I10.10)

   ! Enter the main time stepping loop, which uses the forward Euler scheme.
   loopBool = .TRUE.
   exitThisStep = .FALSE.
   tloop: DO WHILE (loopBool)
      ! Fill the boundary conditions.
      CALL BoundaryConditions(n, u)
      !
      ! Calculate the RHS with the data at the start of the time step.
      CALL RHS(n, u, dx, Lu)
      !
      ! Determine u at the next step.
      u(1:n) = u(1:n) + dt*Lu(1:n)

      ! Increment the time and step counters.
      t = t + dt
      nadv = nadv + 1_IWP

      ! Calculate the next time step.
      CALL CalculateTimeStep(n, u, dx, cfl, dt)

      ! Check whether or not we will write out the data to file.
      IF (writeBool) THEN
         CALL WriteData(n, x, u, nadv)
         writeTime = writeTime + writePeriod
         writeBool = .FALSE.
      ELSE
         IF (t + dt >= writeTime) THEN
            writeBool = .TRUE.
         END IF
      END IF

      ! Check whether or not we will print information to the user.
      IF (printBool .OR. exitThisStep) THEN
         printTime = printTime + printPeriod
         printBool = .FALSE.
         WRITE(*,20) 'Simulation step number: ', nadv, &
                     '; Simulation time: ', t, &
                     '; Max: ', MAXVAL(u(1:n)), &
                     '; Min: ', MINVAL(u(1:n))
         20 FORMAT (A,I10.10,A,ES15.8,A,ES15.8,A,ES15.8)
      ELSE
         IF (t + dt >= printTime) THEN
            printBool = .TRUE.
         END IF
      END IF

      ! Evaluate the exit conditions.
      IF (exitThisStep) THEN
         EXIT tloop
      END IF
      !
      ! Adjust the time step so the desired end time is exactly reached.
      IF (t + dt > tend) THEN
         dt = tend - t
         exitThisStep = .TRUE.
      END IF
   END DO tloop

   ! Write out the final simulation step and time to the user.
   WRITE(*,100) ''
   WRITE(*,130) 'Final step:', nadv
   WRITE(*,120) 'Final time:', t

   ! Write out the final solution.
   CALL WriteData(n, x, u, nadv)

CONTAINS

   !> Subroutine to calculate the RHS of the ODE system.
   !!
   !> @param[in] n Number of grid points.
   !> @param[in] u Current solution array.
   !> @param[in] dx Local grid size.
   !> @param[out] Lu RHS of the ODE system.
   SUBROUTINE RHS(n, u, dx, Lu)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n
      REAL(KIND=RWP),DIMENSION(0:n),INTENT(IN) :: u
      REAL(KIND=RWP),INTENT(IN) :: dx
      REAL(KIND=RWP),DIMENSION(1:n),INTENT(OUT) :: Lu
      ! Local variables.
      ! Looping index.
      INTEGER(KIND=IWP) :: i

      ! Form the RHS using an upwinded FD approximation to the derivative.
      DO i = 1, n
         Lu(i) = -u(i)*(u(i) - u(i-1))/dx
      END DO
   END SUBROUTINE RHS

   !> Subroutine to apply periodic boundary conditions.
   !!
   !! In the finite difference framework, the first and last points are the
   !! same, so the point copied to the left ghost layer is one left of the
   !! right-most point in the domain.
   !!
   !! Again, because an upwinded derivative is used in this code, only one
   !! ghost layer is copied.
   !!
   !> @param[in] n Number of finite difference points.
   !> @param[in,out] u Solution array.
   SUBROUTINE BoundaryConditions(n, u)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n
      REAL(KIND=RWP),DIMENSION(0:n),INTENT(INOUT) :: u

      ! Set periodic boundary conditions for finite difference.
      u(0) = u(n-1)
   END SUBROUTINE BoundaryConditions

   !> Subroutine to calculate the time step.
   !!
   !! We use the CFL criterion to evaluate a stable time step.
   !!
   !> @param[in] n Number of finite difference points.
   !> @param[in] u Current solution array.
   !> @param[in] dx Local grid size.
   !> @param[in] cfl The desired CFL number.
   !> @param[out] dt Time step satisfying the CFL criterion.
   SUBROUTINE CalculateTimeStep(n, u, dx, cfl, dt)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n
      REAL(KIND=RWP),DIMENSION(0:n),INTENT(IN) :: u
      REAL(KIND=RWP),INTENT(IN) :: dx, cfl
      REAL(KIND=RWP),INTENT(OUT) :: dt

      ! Use the minimum dt for stability.
      dt = MINVAL(cfl*dx/ABS(u(1:n)))
   END SUBROUTINE CalculateTimeStep

   !> Routine to write the solution to file.
   !!
   !> @param[in] n Number of points in the domain.
   !> @param[in] x Locations of the points.
   !> @param[in] u Solution to be written.
   !> @param[in] nadv Current simulation step.
   SUBROUTINE WriteData(n, x, u, nadv)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n, nadv
      REAL(KIND=RWP),DIMENSION(1:n),INTENT(IN) :: x
      REAL(KIND=RWP),DIMENSION(0:n),INTENT(IN) :: u
      ! Local variables.
      ! Output file name.
      CHARACTER(LEN=256) :: fname
      ! Looping index.
      INTEGER(KIND=IWP) :: i

      ! Form the output file name.
      WRITE(fname,10) nadv, '.dat'
      10 FORMAT (I12.12,A)

      ! Loop over the cells and record the data to file.
      OPEN(UNIT=100,FILE=fname,FORM="FORMATTED",ACTION="WRITE",STATUS="REPLACE")
      DO i = 1, n
         WRITE(100,20) x(i), u(i)
         20 FORMAT (ES15.8,2X,ES15.8)
      END DO
      CLOSE(UNIT=100)
   END SUBROUTINE WriteData

END PROGRAM FDUpwind_p

