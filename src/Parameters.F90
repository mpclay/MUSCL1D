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
!> @file Parameters.F90
!> @author Matthew Clay
!> @brief Working parameters for the fortran enviornment.
MODULE Parameters_m

   ! Required modules.
   USE ISO_FORTRAN_ENV,ONLY: REAL64, INT32

   !> Working precision for reals.
   INTEGER,PARAMETER,PUBLIC :: IWP = INT32
   !> Working precision for integers.
   INTEGER,PARAMETER,PUBLIC :: RWP = REAL64
   !> Pi.
   REAL(KIND=RWP),PARAMETER,PUBLIC :: PI = ACOS(-1.0_RWP)

END MODULE Parameters_m

