#! /usr/bin/python

################################################################################
# This file is part of Shadowfax
# Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# This script reads a parameter file and initial condition file used for a
# Gadget2 simulation, and converts it into an initial condition file and
# parameter file that can be used with Shadowfax.

# First: check if we have h5py and numpy. Without them, we cannot work.
try:
  import h5py
except ImportError:
  print "It seems you do not have h5py. However, this script cannot operate " \
        "without it. Please install h5py and try again."
  exit()
try:
  import numpy as np
except ImportError:
  print "It seems you do not have numpy. However, this script cannot operate " \
        "without it. Please install numpy and try again."
  exit()

# Now make sure we have the necessary command line arguments: the name of the
# initial condition file and the name of the parameter file.
import argparse
parser = argparse.ArgumentParser(
           description = "Convert Gadget2 initial conditions to Shadowfax " \
                         "initial conditions.")
parser.add_argument("icfile", help = "Gadget2 initial condition file (using " \
                                     "snapshot format 3 (HDF5)")
parser.add_argument("paramfile", help = "Gadget2 parameter file")
parser.add_argument("--no_questions_asked", action = "store_true",
                    help = "Do not prompt for user and use default values " \
                           "for everything")

args = parser.parse_args()

icfile_name = args.icfile
paramfile_name = args.paramfile
if args.no_questions_asked:
  user_input = False
else:
  user_input = True

# check if the files exist and can be opened for reading
import os.path
if not os.path.isfile(icfile_name):
  print "Error! Initial condition file \"{icfile_name}\" not found!".format(
          icfile_name = icfile_name)
  exit()
if not os.path.isfile(paramfile_name):
  print "Error! Parameter file \"{paramfile_name}\" not found!".format(
          paramfile_name = paramfile_name)
  exit()

try:
  icfile = h5py.File(icfile_name, "r")
except Exception:
  print "Error opening initial condition file \"{icfile_name}\"!".format(
          icfile_name = icfile_name)
  exit()

try:
  paramfile = open(paramfile_name, "r")
except Exception:
  print "Error opening parameter file \"{paramfile_name}\"!".format(
          paramfile_name = paramfile_name)
  exit()

# first parse the parameter file and store its contents in a dictionary
params = {}
paramlines = paramfile.readlines()
for line in paramlines:
  # skip empty lines
  line = line.strip()
  if line:
    # skip comment lines
    if not line[0] == '%':
      paramrow = line.split()
      params[paramrow[0]] = paramrow[1]

# now open the initial condition file and collect general data
# we do not want to store all particle data in memory just yet, as the file
# might be quite big
header = {}
for item in icfile["/Header"].attrs.items():
  header[item[0]] = item[1]

groups = {}
for i in range(6):
  parttype = "PartType{i}".format(i = i)
  if parttype in icfile.keys():
    groups[parttype] = True
  else:
    groups[parttype] = False

# if we have gas, check if all primitive variables are available
if groups["PartType0"]:
  has_density = ("Density" in icfile["/PartType0"].keys())

# we currently do not support more than one DM particle type
if groups["PartType2"] or groups["PartType3"] or groups["PartType4"] or \
   groups["PartType5"]:
  print "Particles of a type other than 0 or 1 were found. We do currently " \
        "not support these."
  exit()

# if we have DM, check if masses are in a dataset
# if not, we read them from the header
if groups["PartType1"]:
  has_mass = ("Masses" in icfile["/PartType1"].keys())

# check if we have a box size
if header["BoxSize"]:
  boxsize = header["BoxSize"]
elif float(params["BoxSize"]):
  boxsize = float(params["BoxSize"])
else:
  print "No BoxSize found. The script will analyze all particle positions " \
        "and propose a BoxSize based on this."
  boxsize = 0.

# we need to ask some questions. However, if user_input is False, we do not
# ask them and just use the default values.

# default values:
shadowfax_paramname = paramfile_name.replace(".param", "_shadowfax.ini")
shadowfax_icname = icfile_name.replace(".hdf5", "_shadowfax.hdf5")
snapshot_prefix = "snapshot_"
riemannsolvertype = "Exact"
adiabatic_index = 5./3.
gravity_on = "true"
internal_units = "galactic"
output_units = "galactic"

# questions:
if user_input:
  input = raw_input("Enter a name for the Shadowfax parameter file (default: " \
                    "\"{shadowfax_paramname}\"): ".format(
                      shadowfax_paramname = shadowfax_paramname))
  if input:
    shadowfax_paramname = input
    if not shadowfax_paramname.endswith(".ini"):
      shadowfax_paramname += ".ini"
  print "Shadowfax parameter file name will be " \
        "\"{shadowfax_paramname}\"".format(
          shadowfax_paramname = shadowfax_paramname)

  input = raw_input("Enter a name for the Shadowfax initial condition file " \
                    "(default: \"{shadowfax_icname}\"): ".format(
                      shadowfax_icname = shadowfax_icname))
  if input:
    shadowfax_icname = input
    if not shadowfax_icname.endswith(".hdf5"):
      shadowfax_icname += ".hdf5"
  print "Shadowfax initial condition file name will be " \
        "\"{shadowfax_icname}\"".format(
          shadowfax_icname = shadowfax_icname)

  input = raw_input("Enter a prefix for the snapshot files that will be " \
                    "written out during the simulation (default: " \
                    "\"{snapshot_prefix}\"): ".format(
                      snapshot_prefix = snapshot_prefix))
  if input:
    snapshot_prefix = input
  print "Snapshot prefix will be \"{snapshot_prefix}\"".format(
          snapshot_prefix = snapshot_prefix)

  if groups["PartType0"]:
    input = raw_input("Enter the type of Riemann solver you want to use " \
                      "(possible values: Exact/TRRS/HLLC, default: " \
                      "\"{riemannsolvertype}\": ".format(
                        riemannsolvertype = riemannsolvertype))
    if input:
      riemannsolvertype = input
      if not riemannsolvertype in ["Exact", "HLLC", "TRRS"]:
        print "Error! Unknow Riemann solver type:", riemannsolvertype
        exit()
    print "Using Riemann solver type \"{riemannsolvertype}\"".format(
            riemannsolvertype = riemannsolvertype)

  if groups["PartType0"]:
    input = raw_input("Enter a value for the adiabatic index (default: " \
                      "{adiabatic_index}): ".format(
                        adiabatic_index = adiabatic_index))
    if input:
      adiabatic_index = float(input)
      if adiabatic_index <= 1.:
        print "Error! Invalid adiabatic index:", adiabatic_index
        exit()
    print "Using adiabatic index", adiabatic_index

  input = raw_input("Should the simulation use gravity? (true/false, " \
                    "default: {gravity_on}) ".format(gravity_on = gravity_on))
  if input:
    gravity_on = input
    if not gravity_on in ["true", "false"]:
      print "Error! Invalid choice:", gravity_on
      exit()
  if gravity_on == "true":
    print "Running simulation with gravity"
  else:
    print "Running simulation without gravity"

  input = raw_input("Enter the internal units that should be used (CGS/SI/" \
                    "galactic, default: {internal_units}): ".format(
                      internal_units = internal_units))
  if input:
    internal_units = input
    if not internal_units in ["CGS", "SI", "galactic"]:
      print "Error! Invalid choice of internal units:", internal_units
      exit()
  print "Using internal units:", internal_units

  input = raw_input("Enter the output units that should be used (CGS/SI/" \
                    "galactic, default: {output_units}): ".format(
                      output_units = output_units))
  if input:
    output_units = input
    if not output_units in ["CGS", "SI", "galactic"]:
      print "Error! Invalid choice of output units:", output_units
      exit()
  print "Using output units:", output_units

# the hard part: find out how to convert from gadget units to the requested
# units
gadget_unit_length_in_cgs = float(params["UnitLength_in_cm"])
gadget_unit_mass_in_cgs = float(params["UnitMass_in_g"])
gadget_unit_velocity_in_cgs = float(params["UnitVelocity_in_cm_per_s"])
gadget_unit_time_in_cgs = gadget_unit_length_in_cgs / \
                          gadget_unit_velocity_in_cgs
if internal_units == "CGS":
  shadowfax_unit_length_in_cgs = 1.
  shadowfax_unit_mass_in_cgs = 1.
  shadowfax_unit_time_in_cgs = 1.
elif internal_units == "SI":
  shadowfax_unit_length_in_cgs = 100.
  shadowfax_unit_mass_in_cgs = 1000.
  shadowfax_unit_time_in_cgs = 1.
elif internal_units == "galactic":
  shadowfax_unit_length_in_cgs = 3.08567758e21
  shadowfax_unit_mass_in_cgs = 1.9891e33
  shadowfax_unit_time_in_cgs = 3.154e16
else:
  print "Unknown unit type:", internal_units
  exit()
# multiply length values with this constant to convert them
unit_length_gadget_to_shadowfax = gadget_unit_length_in_cgs / \
                                  shadowfax_unit_length_in_cgs
# multiply time values with this constant to convert them
unit_time_gadget_to_shadowfax = gadget_unit_time_in_cgs / \
                                shadowfax_unit_time_in_cgs

# write the Shadowfax parameter file
shadowfax_paramfile = open(shadowfax_paramname, 'w')
shadowfax_paramfile.write("; Parameter file automatically generated by " \
                          "convert_gadget2_ic.py\n")

# write Time block
shadowfax_paramfile.write("[Time]\n")
maxtime = float(params["TimeMax"])*unit_time_gadget_to_shadowfax
shadowfax_paramfile.write("MaxTime = {maxtime}\n".format(
                            maxtime = maxtime))
shadowfax_paramfile.write("GlobalTimestep = false\n")
maxtimestep = float(params["MaxSizeTimestep"])*unit_time_gadget_to_shadowfax
shadowfax_paramfile.write("MaxTimeStep = {maxtimestep}\n".format(
                            maxtimestep = maxtimestep))
mintimestep = float(params["MinSizeTimestep"])*unit_time_gadget_to_shadowfax
shadowfax_paramfile.write("MinTimeStep = {mintimestep}\n".format(
                            mintimestep = mintimestep))
shadowfax_paramfile.write("TreeTime = false\n\n")

# write Snapshots block
shadowfax_paramfile.write("[Snapshots]\n")
shadowfax_paramfile.write("BaseName = {prefix}\n".format(
                            prefix = snapshot_prefix))
snaptime = float(params["TimeBetSnapshot"])*unit_time_gadget_to_shadowfax
shadowfax_paramfile.write("SnapTime = {snaptime}\n".format(
                            snaptime = snaptime))
# we do need to calculate the index of the first snapshot...
firstsnap_interval = float(params["TimeOfFirstSnapshot"]) - \
                     float(params["TimeBegin"])
firstsnap = int(firstsnap_interval / float(params["TimeBetSnapshot"]) )
shadowfax_paramfile.write("FirstSnap = {firstsnap}\n".format(
                            firstsnap = firstsnap))
shadowfax_paramfile.write("OutputDir = {outputdir}\n".format(
                            outputdir = params["OutputDir"]))
shadowfax_paramfile.write("Type = Gadget\n")
shadowfax_paramfile.write("PerNodeOutput = false\n\n")

# write IC block
shadowfax_paramfile.write("[IC]\n")
shadowfax_paramfile.write("FileName = {filename}\n".format(
                            filename = shadowfax_icname))
shadowfax_paramfile.write("Type = Gadget\n\n")

# write RiemannSolver block
shadowfax_paramfile.write("[RiemannSolver]\n")
shadowfax_paramfile.write("Type = {type}\n".format(
                            type = riemannsolvertype))
shadowfax_paramfile.write("Tolerance = 1.e-8\n")
shadowfax_paramfile.write("CutOff = -5.\n")
shadowfax_paramfile.write("CFL = 0.4\n\n")

# write Hydro block
shadowfax_paramfile.write("[Hydro]\n")
shadowfax_paramfile.write("Gamma = {gamma}\n\n".format(
                            gamma = adiabatic_index))

# write the Gravity block
shadowfax_paramfile.write("[Gravity]\n")
shadowfax_paramfile.write("Gravity = {gravity}\n".format(
                            gravity = gravity_on))
if groups["PartType0"]:
  hgas = float(params["SofteningGas"])
else:
  hgas = 0.
if groups["PartType1"]:
  hdm = float(params["SofteningHalo"])
else:
  hdm = hgas
if not hgas == hdm:
  hmin = min(hgas, hdm)
  print "Warning: different softening lengths found for gas ({hgas}) and " \
        "DM ({hdm}). Using the smallest value ({hmin})".format(
          hgas = hgas, hdm = hdm, hmin = hmin)
else:
  hmin = hdm
shadowfax_paramfile.write("Softening = {softening}\n".format(
                            softening = hmin*unit_length_gadget_to_shadowfax))
eta = float(params["ErrTolIntAccuracy"])/2.8
shadowfax_paramfile.write("Eta = {eta}\n".format(
                            eta = eta))
shadowfax_paramfile.write("Alpha = {alpha}\n\n".format(
                            alpha = float(params["ErrTolForceAcc"])))

# write Voronoi block
shadowfax_paramfile.write("[Voronoi]\n")
shadowfax_paramfile.write("Tolerance = 1.e-9\n\n")

# write Tree block
shadowfax_paramfile.write("[Tree]\n")
shadowfax_paramfile.write("EwaldSize = 64\n")
shadowfax_paramfile.write("EwaldAlpha = 2.\n\n")

# write Memory block
shadowfax_paramfile.write("[Memory]\n")
shadowfax_paramfile.write("MaximumSize = {size} MB\n\n".format(
                            size = int(params["BufferSize"])))

# write Code block
shadowfax_paramfile.write("[Code]\n")
shadowfax_paramfile.write("RestartTime = {restarttime}\n\n".format(
                          restarttime = float(params["CpuTimeBetRestartFile"])))

# write Units block
shadowfax_paramfile.write("[Units]\n")
shadowfax_paramfile.write("InternalUnits = {internal_units}\n".format(
                            internal_units = internal_units))
shadowfax_paramfile.write("OutputUnits = {output_units}\n\n".format(
                            output_units = output_units))

# write Physics block
shadowfax_paramfile.write("[Physics]\n")
if float(params["GravityConstantInternal"]) == 1.:
  real_physics = "false"
else:
  real_physics = "true"
shadowfax_paramfile.write("RealPhysics = {real_physics}\n".format(
                            real_physics = real_physics))
shadowfax_paramfile.write("Cooling = false\n")
# this value is ignored without cooling
shadowfax_paramfile.write("MeanMolWeight = 1.219512195\n")
shadowfax_paramfile.write("StarFormation = false\n")
shadowfax_paramfile.write("StellarFeedback = false\n")

# do not write Cooling, StarFormation, StellarFeedback and Cosmology blocks
# we use default values for these

# find out the dimensions of the box (if necessary)
if groups["PartType0"]:
  coords = np.array(icfile["/PartType0/Coordinates"])
  gmin = np.amin(coords, axis = 0)
  gmax = np.amax(coords, axis = 0)
if groups["PartType1"]:
  coords = np.array(icfile["/PartType1/Coordinates"])
  dmin = np.amin(coords, axis = 0)
  dmax = np.amax(coords, axis = 0)
if groups["PartType0"] and groups["PartType1"]:
  xmin = np.minimum(gmin, dmin)
  xmax = np.maximum(gmax, dmax)
elif groups["PartType0"]:
  xmin = gmin
  xmax = gmax
else:
  xmin = dmin
  xmax = dmax
if boxsize == 0.:
  boxsize = np.ceil(xmax - xmin)
print "Box size will be", boxsize
boxcenter = 0.5*(xmin+xmax)
print "Box center is currently", boxcenter
translation = 0.5*boxsize - boxcenter
print "All particles will be translated by", translation

# write the Shadowfax IC file
shadowfax_icfile = h5py.File(shadowfax_icname, 'w')

# write Header
grp = shadowfax_icfile.create_group("/Header")
header["BoxSize"] = boxsize
for key, value in header.iteritems():
  grp.attrs[key] = value

# write Runtime parameters
grp = shadowfax_icfile.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = int(params["PeriodicBoundariesOn"])

# write Units
# since we use the actual units used by Gadget, we do not have to convert any
# values in the initial condition file
# the conversions will be done by Shadowfax when it reads in the IC
grp = shadowfax_icfile.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = gadget_unit_length_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = gadget_unit_mass_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = gadget_unit_time_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

# write Particle group(s)
if groups["PartType0"]:
  grp = shadowfax_icfile.create_group("/PartType0")
  grp.create_dataset("Coordinates",
                     data = np.array(icfile["/PartType0/Coordinates"]) + \
                            translation,
                     dtype = 'd')
  grp.create_dataset("Velocities",
                     data = np.array(icfile["/PartType0/Velocities"]),
                     dtype = 'd')
  grp.create_dataset("Masses",
                     data = np.array(icfile["/PartType0/Masses"]),
                     dtype = 'd')
  if has_density:
    grp.create_dataset("Density",
                       data = np.array(icfile["/PartType0/Density"]),
                       dtype = 'd')
  grp.create_dataset("SmoothingLength",
                     data = np.array(icfile["/PartType0/SmoothingLength"]),
                     dtype = 'd')
  grp.create_dataset("InternalEnergy",
                     data = np.array(icfile["/PartType0/InternalEnergy"]),
                     dtype = 'd')
  grp.create_dataset("ParticleIDs",
                     data = np.array(icfile["/PartType0/ParticleIDs"]),
                     dtype = 'L')

if groups["PartType1"]:
  grp = shadowfax_icfile.create_group("PartType1")
  grp.create_dataset("Coordinates",
                     data = np.array(icfile["/PartType1/Coordinates"]) + \
                            translation,
                     dtype = 'd')
  grp.create_dataset("Velocities",
                     data = np.array(icfile["/PartType1/Velocities"]),
                     dtype = 'd')
  if has_mass:
    masses = np.array(icfile["/PartType1/Masses"])
  else:
    masses = np.zeros(header["NumPart_ThisFile"][1]) + header["MassTable"][1]
  grp.create_dataset("Masses", data = masses, dtype = 'd')
  grp.create_dataset("ParticleIDs",
                     data = np.array(icfile["/PartType0/ParticleIDs"]),
                     dtype = 'L')

# print some final instructions
print ""
print "Initial conditions successfully generated"
print "Run the simulation using the following command:"
if groups["PartType0"] and not has_density:
  print "\t../shadowfax3d --read_mass --params", shadowfax_paramname
else:
  print "\t../shadowfax3d --params", shadowfax_paramname
