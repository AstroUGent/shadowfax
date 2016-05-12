################################################################################
# This file is part of Shadowfax
# Copyright (C) 2015 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# Shadowfax is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Shadowfax is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Shadowfax. If not, see <http://www.gnu.org/licenses/>.
################################################################################

find_package(Git)
execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags
                OUTPUT_VARIABLE GIT_BUILD_STRING
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND date "+%d/%m/%y"
                OUTPUT_VARIABLE DATE_STRING
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND date "+%H:%M:%S"
                OUTPUT_VARIABLE TIME_STRING
                OUTPUT_STRIP_TRAILING_WHITESPACE)
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(RELEASE_BUILD True)
endif()

configure_file(${PROJECT_SOURCE_DIR}/src/ShadowfaxVersion.hpp.in
               ${PROJECT_BINARY_DIR}/inc/ShadowfaxVersion.hpp)
