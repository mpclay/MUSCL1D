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
!> @file MUSCL.F90
!> @author Matthew Clay
!> @brief Solve the Burgers' equation using the MUSCL scheme.
PROGRAM MUSCL_p

   ! Required modules.
   USE Parameters_m,ONLY: RWP, IWP, PI
   USE Flux_m,ONLY: MINMOD, MONOTONIZED_CENTRAL, SUPERBEE, SetLimiter, Flux

   IMPLICIT NONE

   ! User inputs for the simulation.
   !
   !> Length of the domain.
   REAL(KIND=RWP),PARAMETER :: L = 1.0_RWP
   !> Number of finite volume cells.
   INTEGER(KIND=IWP),PARAMETER :: n = 4096_IWP
   !> CFL for time integration.
   REAL(KIND=RWP),PARAMETER :: cfl = 0.1_RWP
   !> End time for the simulation.
   REAL(KIND=RWP),PARAMETER :: tend = 1.0_RWP
   !> Time frequency with which to write data files.
   REAL(KIND=RWP),PARAMETER :: writePeriod = 0.05_RWP
   !> Time frequency with which to print information to the user.
   REAL(KIND=RWP),PARAMETER :: printPeriod = 0.05_RWP
   !> Which limiter to use.
   INTEGER(KIND=IWP),PARAMETER :: limiter = MINMOD

   ! Variables required for simulation.
   !
   !> The computational grid. Note that for a given cell 'i', x(i) is the x
   !! location of the right cell face.
   REAL(KIND=RWP),DIMENSION(0:n) :: x
   !> Locations of the cell centers.
   REAL(KIND=RWP),DIMENSION(1:n) :: xc
   !> The local grid spacing.
   REAL(KIND=RWP) :: dx
   !> Solution array at the start and end of the time step. We include 1 ghost
   !! layer for the second order scheme.
   REAL(KIND=RWP),DIMENSION(0:n+1) :: u0
   !> Solution array for the first intermediate step during time integration.
   REAL(KIND=RWP),DIMENSION(0:n+1) :: u1
   !> Solution array for the second intermediate step during time integration.
   REAL(KIND=RWP),DIMENSION(1:n) :: u2
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
   !> Total mass in the domain.
   REAL(KIND=RWP) :: mass

   ! Zero out necessary variables.
   t = 0.0_RWP
   nadv = 0_IWP

   ! Set up the flux module.
   CALL SetLimiter(limiter)

   ! Fill in the computational grid.
   dx = L/REAL(n, RWP)
   DO i = 0, n-1
      x(i) = REAL(i, RWP)*dx
      xc(i+1) = dx/2.0_RWP + REAL(i, RWP)*dx
   END DO
   x(n) = L

   ! Form the initial conditions.
   DO i = 1, n
      u0(i) = 1_RWP + L/(2.0_RWP*PI*dx)*(COS(2.0_RWP*PI*x(i-1)/L) - &
                                         COS(2.0_RWP*PI*x(i)/L))
   END DO
   !
   ! Write out the initial conditions.
   CALL WriteData(n, xc, u0, nadv)
   !
   ! Calculate the total mass of the initial conditions.
   CALL TotalMass(n, dx, u0, mass)

   ! Calculate the initial time step.
   CALL CalculateTimeStep(n, u0, dx, cfl, dt)

   ! Print some information to the user.
   WRITE(*,100) '---------------------------------------------------------'
   WRITE(*,100) 'MUSCL1D: Using the MUSCL scheme to solve Burgers equation'
   WRITE(*,100) '---------------------------------------------------------'
   WRITE(*,100) ''
   WRITE(*,120) 'Domain length:', L
   WRITE(*,130) 'Num. of cells:', n
   WRITE(*,120) 'CFL number:', cfl
   WRITE(*,120) 'End time:', tend
   WRITE(*,120) 'Initial mass:', mass
   WRITE(*,100) ''
   100 FORMAT (A)
   110 FORMAT (A,T17,A)
   120 FORMAT (A,T16,ES15.8)
   130 FORMAT (A,T17,I10.10)

   ! Enter the main time stepping loop, which uses the improved Euler scheme.
   loopBool = .TRUE.
   exitThisStep = .FALSE.
   tloop: DO WHILE (loopBool)
      ! Fill the boundary conditions.
      CALL BoundaryConditions(n, u0)
      !
      ! Calculate the RHS with the data at the start of the time step.
      CALL Flux(n, u0, dx, Lu)
      !
      ! Determine the intermediate Euler stage.
      u1(1:n) = u0(1:n) + dt*Lu(1:n)

      ! Fill the boundary conditions.
      CALL BoundaryConditions(n, u1)
      !
      ! Calculate the RHS using the intermediate stage.
      CALL Flux(n, u1, dx, Lu)
      !
      ! Determine the next intermediate solution.
      u2(1:n) = u1(1:n) + dt*Lu(1:n)

      ! Average the two intermediate stages to get the next solution.
      u0(1:n) = 0.5_RWP*(u1(1:n) + u2(1:n))

      ! Increment the time and step counters.
      t = t + dt
      nadv = nadv + 1_IWP

      ! Calculate the next time step.
      CALL CalculateTimeStep(n, u0, dx, cfl, dt)

      ! Check whether or not we will write out the data to file.
      IF (writeBool) THEN
         CALL WriteData(n, xc, u0, nadv)
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
         CALL TotalMass(n, dx, u0, mass)
         WRITE(*,20) 'Simulation step number: ', nadv, &
                     '; Simulation time: ', t, &
                     '; Max: ', MAXVAL(u0(1:n)), &
                     '; Min: ', MINVAL(u0(1:n)), &
                     '; Total: ', mass
         20 FORMAT (A,I10.10,A,ES15.8,A,ES15.8,A,ES15.8,A,ES15.8)
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

   ! Write out the final simulation step, time, and mass to the user.
   CALL TotalMass(n, dx, u0, mass)
   WRITE(*,100) ''
   WRITE(*,130) 'Final step:', nadv
   WRITE(*,120) 'Final time:', t
   WRITE(*,120) 'Final mass:', mass

   ! Write out the final solution.
   CALL WriteData(n, xc, u0, nadv)

CONTAINS

   !> Subroutine to apply periodic boundary conditions.
   !!
   !! In the finite volume framework, periodic boundary conditions are enforced
   !! by gluing the leftmost and rightmost faces together.
   !!
   !> @param[in] n Number of finite volume cells.
   !> @param[in,out] u Solution array.
   SUBROUTINE BoundaryConditions(n, u)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n
      REAL(KIND=RWP),DIMENSION(0:n+1),INTENT(INOUT) :: u

      ! Set periodic boundary conditions for finite volume.
      u(0) = u(n)
      u(n+1) = u(1)
   END SUBROUTINE BoundaryConditions

   !> Subroutine to calculate the time step.
   !!
   !! We use the CFL criterion to evaluate a stable time step.
   !!
   !> @param[in] n Number of finite volume cells.
   !> @param[in] u Current solution array.
   !> @param[in] dx Local grid size.
   !> @param[in] cfl The desired CFL number.
   !> @param[out] dt Time step satisfying the CFL criterion.
   SUBROUTINE CalculateTimeStep(n, u, dx, cfl, dt)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n
      REAL(KIND=RWP),DIMENSION(0:n+1),INTENT(IN) :: u
      REAL(KIND=RWP),INTENT(IN) :: dx, cfl
      REAL(KIND=RWP),INTENT(OUT) :: dt

      ! Use the minimum dt for stability.
      dt = MINVAL(cfl*dx/ABS(u(1:n)))
   END SUBROUTINE CalculateTimeStep

   !> Subroutine to calculate the total mass in the domain.
   !!
   !> @param[in] n Number of FV cells.
   !> @param[in] dx Local grid spacing.
   !> @param[in] u Solution array.
   !> @param[out] mass Total mass
   SUBROUTINE TotalMass(n, dx, u, mass)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n
      REAL(KIND=RWP),INTENT(IN) :: dx
      REAL(KIND=RWP),DIMENSION(0:n+1),INTENT(IN) :: u
      REAL(KIND=RWP),INTENT(OUT) :: mass
      ! Local variables.
      INTEGER(KIND=IWP) :: i

      ! Increment the mass to get the total mass.
      mass = 0.0_RWP
      DO i = 1, n
         mass = mass + dx*u(i)
      END DO
   END SUBROUTINE TotalMass

   !> Routine to write the solution to file.
   !!
   !> @param[in] n Number of cells in the domain.
   !> @param[in] xc Locations of cell centers.
   !> @param[in] us Solution to be written.
   !> @param[in] nadv Current simulation step.
   SUBROUTINE WriteData(n, xc, us, nadv)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n, nadv
      REAL(KIND=RWP),DIMENSION(1:n),INTENT(IN) :: xc
      REAL(KIND=RWP),DIMENSION(0:n+1),INTENT(IN) :: us
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
         WRITE(100,20) xc(i), us(i)
         20 FORMAT (ES15.8,2X,ES15.8)
      END DO
      CLOSE(UNIT=100)
   END SUBROUTINE WriteData

END PROGRAM MUSCL_p

