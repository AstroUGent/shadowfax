# -*- coding: utf-8 -*-

################################################################################
# This file is part of Shadowfax
# Copyright (C) 2016 Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
#                    Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

"""
Created on Fri Feb 26 15:54:25 2016

@author: Yorick Van Den Bossche (yorick.vandenbossche@ugent.be)
"""
import os
import shutil

files = os.listdir(".")
for file in files:
    if file.startswith("RadLoss_-99.00"):
        file2 = file.split("_")
        file2[2] = "-0.55"
        newfile = ""
        for part in file2:
            newfile = newfile + part + "_"
        newfile = newfile.strip("_")
        shutil.copy(file, newfile)
        
        file2[2] = "0.47"
        newfile = ""
        for part in file2:
            newfile = newfile + part + "_"
        newfile = newfile.strip("_")
        shutil.copy(file, newfile)
