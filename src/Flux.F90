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
!> @file Flux.F90
!> @author Matthew Clay
!> @brief Calculate the RHS of the finite volume scheme.
MODULE Flux_m

   ! Required modules.
   USE Parameters_m,ONLY: IWP, RWP

   IMPLICIT NONE

   ! Available limiters.
   !
   !> Minmod limiter.
   INTEGER(KIND=IWP),PARAMETER,PUBLIC :: MINMOD = 0_IWP
   !> Monotonized central limiter.
   INTEGER(KIND=IWP),PARAMETER,PUBLIC :: MONOTONIZED_CENTRAL = 1_IWP
   !> Superbee limiter.
   INTEGER(KIND=IWP),PARAMETER,PUBLIC :: SUPERBEE = 2_IWP

   !> Which limiter is used during the simulation.
   INTEGER(KIND=IWP),PRIVATE :: limiter

   ! Module procedures.
   PUBLIC :: SetLimiter, Flux
   PRIVATE :: MinModLimiter, MonotonizedCentralLimiter, SuperbeeLimiter

CONTAINS

   !> Routine to set the limiter to use for calculations.
   !!
   !> @param[in] limiter_ The desired limiter.
   SUBROUTINE SetLimiter(limiter_)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: limiter_

      ! Check to make sure it is a valid option.
      SELECT CASE (limiter_)
         CASE (MINMOD)
            limiter = MINMOD
         CASE (MONOTONIZED_CENTRAL)
            limiter = MONOTONIZED_CENTRAL
         CASE (SUPERBEE)
            limiter = SUPERBEE
         CASE DEFAULT
            WRITE(*,100) 'Invalid limiter. Halting.'
            STOP
            100 FORMAT (A)
      END SELECT
   END SUBROUTINE SetLimiter

   !> Routine to calculate the RHS of the finite volume scheme.
   !!
   !! To be clear: indexing for cells starts at 1 and goes to n, while indexing
   !! for faces starts from 0 and goes to n.
   !!
   !! CELLS:    1     2     3          n-1    n
   !!        |-----|-----|-----| ... |-----|-----|
   !! FACES: 0     1     2     3    n-2   n-1    n
   !!
   !> @param[in] n Number of finite volume cells.
   !> @param[in] u Current state vector.
   !> @param[in] dx Local grid size.
   !> @param[in,out] Lu RHS of the system.
   SUBROUTINE Flux(n, u, dx, Lu)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n
      REAL(KIND=RWP),DIMENSION(0:n+1),INTENT(IN) :: u
      REAL(KIND=RWP),INTENT(IN) :: dx
      REAL(KIND=RWP),DIMENSION(1:n),INTENT(OUT) :: Lu
      ! Local variables.
      ! Primitives extrapolated from the left for a given cell face.
      REAL(KIND=RWP),DIMENSION(0:n) :: uL
      ! Primitives extrapolated from the right for a given cell face.
      REAL(KIND=RWP),DIMENSION(0:n) :: uR
      ! Flux calculated with left-extrapolated values.
      REAL(KIND=RWP),DIMENSION(0:n) :: fL
      ! Flux calculated with right-extrapolated values.
      REAL(KIND=RWP),DIMENSION(0:n) :: fR
      ! Upwinded flux.
      REAL(KIND=RWP),DIMENSION(0:n) :: fU
      ! The divided differences for each cell.
      REAL(KIND=RWP),DIMENSION(1:n) :: r
      ! Value of the slope-limiter function for each cell.
      REAL(KIND=RWP),DIMENSION(1:n) :: phi
      ! Maximum wavespeed in the domain, which is used for the LF flux.
      REAL(KIND=RWP) :: alpha
      ! Looping index for cells and faces.
      INTEGER(KIND=IWP) :: i

      ! Calculate the divided differences for each cell.
      DO i = 1, n
         r(i) = (u(i) - u(i-1))/(u(i+1) - u(i))
      END DO

      ! Calculate the limiter function for each cell.
      SELECT CASE (limiter)
         CASE (MINMOD)
            CALL MinModLimiter(n, r, phi)
         CASE (MONOTONIZED_CENTRAL)
            CALL MonotonizedCentralLimiter(n, r, phi)
         CASE (SUPERBEE)
            CALL SuperbeeLimiter(n, r, phi)
      END SELECT

      ! Loop over each cell and form the left and right extrapolated values.
      !
      ! CELL:          i
      !          |-----------|
      ! VALUE: uR(i-1)     uL(i)
      !
      DO i = 1, n
         uR(i-1) = u(i) - 0.5_RWP*phi(i)*(u(i+1) - u(i))
         uL(i) = u(i) + 0.5_RWP*phi(i)*(u(i+1) - u(i))
      END DO

      ! Copy the left-extrapolated values at the right of the domain to the left
      ! face at the beginning of the domain.
      uL(0) = uL(n)
      !
      ! Copy the right-extrapolated values at the left of the domain to the
      ! right face at the end of the domain.
      uR(n) = uR(0)

      ! At each face, calculate fluxes based on the left and right states.
      DO i = 0, n
         fL(i) = 0.5_RWP*uL(i)**2
         fR(i) = 0.5_RWP*uR(i)**2
      END DO

      ! Calculate the maximum wavespeed in the domain for the LF flux.
      alpha = MAXVAL(ABS(u(1:n)))
      !
      ! Combine the separate fluxes on each face with a LF flux function.
      DO i = 0, n
         fU(i) = 0.5_RWP*(fL(i) + fR(i) - alpha*(uR(i) - uL(i)))
      END DO

      ! For each cell, calculate its increment based on the flux difference.
      ! Remember, for each cell, a flux array indexed at the cell index is for
      ! the right face, and the flux array indexed at one less than the cell
      ! index is for the left face.
      DO i = 1, n
         Lu(i) = -(fU(i) - fU(i-1))/dx
      END DO
   END SUBROUTINE Flux

   !> Minmod slope limiter.
   !!
   !> @param[in] n Number of cells in the domain.
   !> @param[in] r Array of divided differences.
   !> @param[out] phi Slope limiting function.
   SUBROUTINE MinModLimiter(n, r, phi)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n
      REAL(KIND=RWP),DIMENSION(1:n),INTENT(IN) :: r
      REAL(KIND=RWP),DIMENSION(1:n),INTENT(OUT) :: phi
      ! Local variables.
      ! Looping index.
      INTEGER(KIND=IWP) :: i

      ! Loop over each cell to calculate the limiter.
      DO i = 1, n
         phi(i) = MAX(0.0_RWP, MIN(1.0_RWP, r(i)))
      END DO
   END SUBROUTINE MinModLimiter

   !> Monotonized-central slope limiter.
   !!
   !> @param[in] n Number of cells in the domain.
   !> @param[in] r Array of divided differences.
   !> @param[out] phi Slope limiting function.
   SUBROUTINE MonotonizedCentralLimiter(n, r, phi)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n
      REAL(KIND=RWP),DIMENSION(1:n),INTENT(IN) :: r
      REAL(KIND=RWP),DIMENSION(1:n),INTENT(OUT) :: phi
      ! Local variables.
      ! Looping index.
      INTEGER(KIND=IWP) :: i

      ! Loop over each cell to calculate the limiter.
      DO i = 1, n
         phi(i) = MAX(0.0_RWP, MIN(2.0_RWP*r(i), &
                                   0.5_RWP*(r(i) + 1.0_RWP), &
                                   2.0_RWP))
      END DO
   END SUBROUTINE MonotonizedCentralLimiter

   !> Superbee slope limiter.
   !!
   !> @param[in] n Number of cells in the domain.
   !> @param[in] r Array of divided differences.
   !> @param[out] phi Slope limiting function.
   SUBROUTINE SuperbeeLimiter(n, r, phi)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n
      REAL(KIND=RWP),DIMENSION(1:n),INTENT(IN) :: r
      REAL(KIND=RWP),DIMENSION(1:n),INTENT(OUT) :: phi
      ! Local variables.
      ! Looping index.
      INTEGER(KIND=IWP) :: i

      ! Loop over each cell to calculate the limiter.
      DO i = 1, n
         phi(i) = MAX(0.0_RWP, MIN(2.0_RWP*r(i), 1.0_RWP), MIN(r(i), 2.0_RWP))
      END DO
   END SUBROUTINE SuperbeeLimiter

END MODULE Flux_m

