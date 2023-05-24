import numpy as np
import os
import shutil

#########################################################################################################
# Task 1 - Find out which calculation type								#
#########################################################################################################

cwd = os.getcwd()
if '/B3LYP/' in cwd:
 print "Preparing input for a B3LYP calculation"
 Func="B3LYP"
if '/mB3LYP/' in cwd:
 print "Preparing input for a mB3LYP calculation"
 Func="mB3LYP"
if '/M062X/' in cwd:
 print "Preparing input for a M062X calculation"
 Func="M062X"
if '/M06-2X/' in cwd:
 print "Preparing input for a M062X calculation"
 Func="M062X"

cwd = os.getcwd()
if '/Water/' in cwd:
 print "Preparing input for a Water calculation"
 Env="Wat"
if '/Vacuum/' in cwd:
 print "Preparing input for a Vacuum calculation"
 Env="Vac"

for file in os.listdir("./../../../../../"):
 if file.endswith('.out'):
  Nam=file
g=open("./../../../../../"+str(Nam), 'r')
loc=g.readlines()
g.close()
for line in range(len(loc)):
 if "Program Version 3.0.3 - RELEASE   -" in loc[line]:
  PathMOExec="/net/hpc-home/llk/Programs/Orca_3.0.3/orca"
  Ver="3_0_3"
  print "Using Orca 3.0.3 for MOs"
 elif "Program Version 4.0.0.2 - RELEASE -" in loc[line]:
  PathMOExec="/net/hpc-home/llk/Programs/Orca_4_0_0_2/orca"
  Ver="4_0_0_2"
  print "Using Orca 4.0.0.2 for MOs"

if Env=="Vac":
 for file in os.listdir("./../../Localize/"):
  if file.endswith('.out'):
   Nam=file
 g=open("./../../Localize/"+str(Nam), 'r')
 loc=g.readlines()
 g.close()
 for line in range(len(loc)):
  if "Program Version 3.0.3 - RELEASE   -" in loc[line]:
   PathLMOExec="/net/hpc-home/llk/Programs/Orca_3.0.3/orca"
   Ver="3_0_3"
   print "Using Orca 3.0.3 for LMOs"
  elif "Program Version 4.0.0.2 - RELEASE -" in loc[line]:
   PathLMOExec="/net/hpc-home/llk/Programs/Orca_4_0_0_2/orca"
   Ver="4_0_0_2"
   print "Using Orca 4.0.0.2 for LMOs"

# Make the MO and LMO folders
os.makedirs("./../../MOs/")
if Env=="Vac":
 os.makedirs("./../../Localize/LMOs/")

# Copy the necessary .gbw and .loc files
for file in os.listdir("./../../../../../"):
 if file.endswith('.gbw'):
  shutil.copyfile("./../../../../../"+file, "./../../MOs/MOs.gbw")
if Env=="Vac":
 for file in os.listdir("./../../Localize/"):
  if file.endswith('.loc'):
   shutil.copyfile("./../../Localize/"+file, "./../../Localize/LMOs/LMOs.loc")

# Construct the necessary input-file-name from information obtained so far
MOInp="./Templates/MOs_"+str(Env)+"/MO_"+str(Func)+"_"+str(Ver)+".inp"
LMOInp="./Templates/LMOs_"+str(Env)+"/MO_"+str(Func)+"_"+str(Ver)+".inp"

# Copy necessary files
shutil.copyfile(MOInp, "./../../MOs/MOs.inp")
if Env=="Vac":
 shutil.copyfile(LMOInp, "./../../Localize/LMOs/LMOs.inp")

# open the xyz file of the structure optimization and graft the content onto the inputs
for file in os.listdir("./../../../../../"):
 if file.endswith('.xyz'):
  Nam=file

g=open("./../../../../../"+str(Nam),'r')
loc=g.readlines()
g.close()

h=open("./../../MOs/MOs.inp",'a')
h.write('* xyz 0 1\n')
for line in range(2,len(loc)):
 h.write(loc[line])
h.write('*')
h.close()

if Env=="Vac":
 l=open("./../../Localize/LMOs/LMOs.inp",'a')
 l.write('* xyz 0 1\n')
 for line in range(2,len(loc)):
  l.write(loc[line])
 l.write('*')
 l.close()

#########################################################################################################
# Task 2 - Generate submission script that starts the correct versions					#
#########################################################################################################

m=open("Submit.sh",'w')

m.write('''#########################################################
#							#
#	Script starting the Orbital calculations	#
#							#
#########################################################''')

m.write('''
cd ./../../MOs/
'''+str(PathMOExec)+''' MOs.inp > MOs.out &
cd -
''')
if Env=="Vac":
 m.write('''#cd ./../../Localize/LMOs/
#'''+str(PathLMOExec)+''' LMOs.inp > LMOs.out &
 #cd -''')
m.close()
