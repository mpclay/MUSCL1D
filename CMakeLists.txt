# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# Copyright (C) 2014 by Matthew Clay <mpclay@gmail.com>

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.9)

PROJECT(MUSCL)
ENABLE_LANGUAGE(Fortran)

# Main source directory.
ADD_SUBDIRECTORY(src)

# Build the main executable.
OPTION(BUILD_MUSCL "Build the MUSCL FV program." OFF)
IF(BUILD_MUSCL)
   ADD_EXECUTABLE(muscl ${CORE_SOURCE_FILES} ${MUSCL_SOURCE_FILES})
   INSTALL(TARGETS muscl DESTINATION ${CMAKE_INSTALL_PREFIX})
ENDIF(BUILD_MUSCL)

# Build the FD upwind code.
OPTION(BUILD_FDUPWIND "Build the upwind FD program." OFF)
IF(BUILD_FDUPWIND)
   ADD_EXECUTABLE(fdupwind ${CORE_SOURCE_FILES} ${FDUPWIND_SOURCE_FILES})
   INSTALL(TARGETS fdupwind DESTINATION ${CMAKE_INSTALL_PREFIX})
ENDIF(BUILD_FDUPWIND)

