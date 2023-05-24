import numpy as np

#    DEFAULTS    #
DOCarbon=True
DOOxygen=False
DONitrogen=False

f=open('MOs.out','r')
loc=f.readlines()
f.close()
NumC=0
NumN=0
NumO=0

#########################################################################################################
# Task 1 - Read in the first file
#########################################################################################################

line=0
while not "* xyz " in loc[line]:
 line=line+1

# The next line is then the xyz line - hence skip the next line and then read until "> *"
line=line+1

while not "> *" in loc[line]:
# print(loc[line].split()[2])
 Read=((loc[line].split(">")[1]).split()[0])
 line=line+1
 if Read=="C":
  NumC=NumC+1
 if Read=="N":
  NumN=NumN+1
 if Read=="O":
  NumO=NumO+1

### Allocate Energy and core MO Vectors (CMO)
CMOTotal=NumC+NumN+NumO 
CMORead=0
# Initialize temporary sorting Matrices
TempEn=np.zeros((CMOTotal))
TempVec=np.zeros((CMOTotal,CMOTotal))
TempLabel=[]

# The next lines are uninteresting - hence read until the MO compositions come up
while not "MOLECULAR ORBITALS" in loc[line]:
 line=line+1

# The next line is the MOLECULAR ORBITALS line - skip next two lines which are empty
line=line+2

#### Read all the CMOs and associated 1s vectors until all CMOs were included ####
while CMORead < CMOTotal:
 # save rewind line
 rewind=line
 
 # Read how many elements can be split up = number of MOs read in simultaneously
 CMOsimul=len(loc[line].split())
 
 # Read the first (next) CMO informations, if it is still a CMO
 for i in range(CMOsimul):
  basis=0
  # Remember - there are n atoms, but we read orbitals from 0 to n-1, hence .lt. instead of .le.!!!
  if CMORead < CMOTotal:
 # Advance a line and then read Energy of the Orbital i
   line=line+1
   TempEn[CMORead]=loc[line].split()[i]
 
 # Advance three more lines until the compositions start
   line=line+3
   
 # Read the lines until next delimiter "--------"
   while not "--------" in loc[line]:
    if "C   1s" in loc[line]:
     TempVec[CMORead, basis]=loc[line].split()[i+2]
     TempLabel.append(loc[line].split()[0])
#     print (loc[line].split()[0]), basis
     basis=basis+1
    if "N   1s" in loc[line]:
     TempVec[CMORead, basis]=loc[line].split()[i+2]
     TempLabel.append(loc[line].split()[0])
#     print (loc[line].split()[0]), basis
     basis=basis+1
    if "O   1s" in loc[line]:
     TempVec[CMORead, basis]=loc[line].split()[i+2]
     TempLabel.append(loc[line].split()[0])
#     print (loc[line].split()[0]), basis
     basis=basis+1
    line=line+1
   CMORead=CMORead+1
 
 # Rewind, as long as we have not reached the last one to store
   if i+1 < CMOsimul:
    line=rewind
 
 # Now, we are in the delimiter line of the next block. Rewind 3 lines to get back to the orbital numbers and start again
 line=line-3

loc=[]

#########################################################################################################
# Task 2 - Analyze the initial localization and store energies and MO indices
#########################################################################################################

#GL abbreviates GreenListed
GLCMO=[]
GLEn=[]
GLLabel=[]

for MO in range(CMOTotal):
 for Bas in range(CMOTotal):
  if abs(TempVec[MO,Bas]) > 0.95:
#   print "Orbital number", MO, "with energy", TempEn[MO], "is strongly localized at", TempLabel[Bas]
   GLLabel.append(TempLabel[Bas])
   GLEn.append(TempEn[MO])
   GLCMO.append(MO)

#########################################################################################################
# Task 3 - Calculate Median of remaining O, C and N Edges
#########################################################################################################

# initialize Edge Medians
OEM=0.0
NEM=0.0
CEM=0.0
RLO=0
RLN=0
RLC=0

# start looking for left-over energies and which atom they are assigned to
for MO in range(CMOTotal):
 for Bas in range(CMOTotal):
  if not MO in GLCMO:
   if max(abs(TempVec[MO,:])) == abs(TempVec[MO,Bas]):
#    print "Orbital", MO, "is mostly localized on", TempLabel[Bas]
    if "C" in TempLabel[Bas]:
     CEM=CEM+TempEn[MO]
     RLC=RLC+1
    if "N" in TempLabel[Bas]:
     NEM=NEM+TempEn[MO]
     RLN=RLN+1
    if "O" in TempLabel[Bas]:
     OEM=OEM+TempEn[MO]
     RLO=RLO+1

# Build the medians
if RLC > 0:
 CEM=CEM/RLC*27.2114
if RLN > 0:
 NEM=NEM/RLN*27.2114
if RLO > 0:
 OEM=OEM/RLO*27.2114

#########################################################################################################
# Task 4 - Read every strongly localized orbital with respective Atom Label
#########################################################################################################

g=open('Loc.out','r')
loc=g.readlines()
g.close()
line=0
LocMOs=0
LMONum=[]
LMOLabel=[]
LMORead=0

while not "Rather strongly localized orbitals" in loc[line]:
 line=line+1

# Jump to the MO's line and read the number of total localized MOs. Head down that many lines
line=line+1
FirstMO=line
LocMOs=int((loc[line].split(":")[0]).split()[1])
line=line+LocMOs

# Now read upwards and store the MO numbers as well as the atom labels
while line > FirstMO-1:
 
 LMONum.append((loc[line].split(":")[0]).split()[1])
 LMOLabel.append((loc[line].split(":")[1]).split()[0])
 line=line-1

TotalMOs=len(LMONum)

#########################################################################################################
# Task 5 - Read the original energies of these N orbitals, since these are the "wrongly assigned ones"  #
#########################################################################################################

# Rewind file and go to the ORBITAL ENERGIES
line=0
while not "ORBITAL ENERGIES" in loc[line]:
 line=line+1

# Jump to the first orbital - these are always in energetic order
line=line+4
OldLabel=[]
OldEn=np.zeros((TotalMOs))

for i in range(TotalMOs):
 OldLabel.append(loc[line].split()[0])
 OldEn[i]=(float(loc[line].split()[3]))
# print "Old Orbital number ", OldLabel[i]," has had energy of ", OldEn[i]
 line=line+1

#########################################################################################################
# Task 6 - Regroup, Reassign and Calculate Shifts atomwise                                              #
#########################################################################################################

CCoreAtomLabel=[]
OCoreAtomLabel=[]
NCoreAtomLabel=[]
CCoreMOLabel=[]
OCoreMOLabel=[]
NCoreMOLabel=[]
CCoreEn=np.zeros((NumC))
OCoreEn=np.zeros((NumO))
NCoreEn=np.zeros((NumN))
CCoreShift=np.zeros((NumC))
OCoreShift=np.zeros((NumO))
NCoreShift=np.zeros((NumN))
CCnt=0
OCnt=0
NCnt=0
NewEn=0

for i in range(TotalMOs):
 print "Probing MO Number", LMONum[i], " that is localized at ", LMOLabel[i]
 for j in range(CMOTotal):
  if TempLabel[j] == LMOLabel[i]:
   if "C" in LMOLabel[i]:
    # if not already assigned, put inside the group
    if LMOLabel[i] not in CCoreAtomLabel:
     CCoreAtomLabel.append(LMOLabel[i])
     CCoreMOLabel.append(LMONum[i])
     CCoreEn[CCnt]=OldEn[i]
     print "The Carbon Atom ", CCoreAtomLabel[CCnt], " was assigned the LMO Number", CCoreMOLabel[CCnt],". It had an original energy of ", OldEn[i],".\n"

     ###proceed with calculating the Shift
     # check if the MO was on a strongly localized Atom before
     if CCoreAtomLabel[CCnt] in GLLabel:
      # check which one it was and assign the respective energy as TempVar NewEn
      for k in range(len(GLLabel)):
       if CCoreAtomLabel[CCnt] in GLLabel[k]:
        NewEn=float(GLEn[k])*27.2114
        CCoreShift[CCnt]=NewEn-CCoreEn[CCnt]
     else:
      NewEn=CEM
      CCoreShift[CCnt]=NewEn-CCoreEn[CCnt]
     
     #continue counting
     CCnt=CCnt+1

   if "O" in LMOLabel[i]:
    # if not already assigned, put inside the group
    if LMOLabel[i] not in OCoreAtomLabel:
     OCoreAtomLabel.append(LMOLabel[i])
     OCoreMOLabel.append(LMONum[i])
     OCoreEn[OCnt]=OldEn[i]
     print "The Oxygen Atom ", OCoreAtomLabel[OCnt], " was assigned the LMO Number", OCoreMOLabel[OCnt],". It had an original energy of ", OldEn[i],".\n"

     ###proceed with calculating the Shift
     # check if the MO was on a strongly localized Atom before
     if OCoreAtomLabel[OCnt] in GLLabel:
      # check which one it was and assign the respective energy as TempVar NewEn
      for k in range(len(GLLabel)):
       if OCoreAtomLabel[OCnt] in GLLabel[k]:
        NewEn=float(GLEn[k])*27.2114
        OCoreShift[OCnt]=NewEn-OCoreEn[OCnt]
     else:
      NewEn=OEM
      OCoreShift[OCnt]=NewEn-OCoreEn[OCnt]
     
     #continue counting
     OCnt=OCnt+1

   if "N" in LMOLabel[i]:
    # if not already assigned, put inside the group
    if LMOLabel[i] not in NCoreAtomLabel:
     NCoreAtomLabel.append(LMOLabel[i])
     NCoreMOLabel.append(LMONum[i])
     NCoreEn[NCnt]=OldEn[i]
     print "The Nitrogen Atom ", NCoreAtomLabel[NCnt], " was assigned the LMO Number", NCoreMOLabel[NCnt],". It had an original energy of ", OldEn[i],".\n"

     ###proceed with calculating the Shift
     # check if the MO was on a strongly localized Atom before
     if NCoreAtomLabel[NCnt] in GLLabel:
      # check which one it was and assign the respective energy as TempVar NewEn
      for k in range(len(GLLabel)):
       if NCoreAtomLabel[NCnt] in GLLabel[k]:
        NewEn=float(GLEn[k])*27.2114
        NCoreShift[NCnt]=NewEn-NCoreEn[NCnt]
     else:
      NewEn=NEM
      NCoreShift[NCnt]=NewEn-NCoreEn[NCnt]
     
     #continue counting
     NCnt=NCnt+1

#########################################################################################################
# Task 8 - Write a script that generate the "false spectra" around the former edge energies
#########################################################################################################

g=open('XASSpec.sh','w')
g.write('''#########################################################
#							#
#	script for generating the "false" XAS spectra	#
#							#
#########################################################''')

g.write('''

XASSPEC(){
Nam="Job.out"
suffold=".absq.dat"
D1="-x0"
D2="-x1"
R1=$D1$1
R2=$D2$2
a=$Nam$suffold
orca_mapspc $Nam ABSQ -eV $R1 $R2 -n10000 -w1.0
mv $a XAS_eV_wrong.dat
}
''')
if DOCarbon==True:
 for i in range(NumC):
 # Calculate the Window of the wrong spectrum
  WinM=-np.round(CCoreEn[i], decimals=0)
  WinL=int(WinM-50)
  WinH=int(WinM+50)
  if WinL < int(0):
   WinL=0
   WinH=100

  g.write('''
cd ../../Single_Carbon/C_Atom_'''+str(i)+'''
XASSPEC '''+str(WinL)+''' '''+str(WinH)+'''
cd -
  ''')

if DOOxygen==True:
 for i in range(NumO):
 # Calculate the Window of the wrong spectrum
  WinM=-np.round(OCoreEn[i], decimals=0)
  WinL=int(WinM-50)
  WinH=int(WinM+50)
  if WinL < int(0):
   WinL=0
   WinH=100

  g.write('''
cd ../../Single_Oxygen/O_Atom_'''+str(i)+'''
XASSPEC '''+str(WinL)+''' '''+str(WinH)+'''
cd -
  ''')

if DONitrogen==True:
 for i in range(NumN):
 # Calculate the Window of the wrong spectrum
  WinM=-np.round(NCoreEn[i], decimals=0)
  WinL=int(WinM-50)
  WinH=int(WinM+50)
  if WinL < int(0):
   WinL=0
   WinH=100
 
  g.write('''
cd ../../Single_Nitrogen/N_Atom_'''+str(i)+'''
XASSPEC '''+str(WinL)+''' '''+str(WinH)+'''
cd -
  ''')

g.close()

#########################################################################################################
# Task 9 - Write the shifts to an awk script
#########################################################################################################

k=open('Shift.sh','w')
k.write('''#########################################################
#							#
#	awk-shift script for shifting XAS spectra	#
#							#
#########################################################''')

if DOCarbon==True:
 for i in range(NumC):
  k.write('''
cd ../../Single_Carbon/C_Atom_'''+str(i))
  if CCoreShift[i] > 0.0:
   Target=np.round(CCoreShift[i], decimals=2)
   k.write('''
awk '{print $1-'''+str(abs(Target))+''', $2/'''+str(abs(CCoreEn[i]))+'''}' XAS_eV_wrong.dat > XAS_shift.dat
cd -
   ''')
  else:
   Target=np.round(CCoreShift[i], decimals=2)
   k.write('''
awk '{print $1+'''+str(abs(Target))+''', $2/'''+str(abs(CCoreEn[i]))+'''}' XAS_eV_wrong.dat > XAS_shift.dat
cd -
   ''')

if DOOxygen==True:
 for i in range(NumO):
  k.write('''
cd ../../Single_Oxygen/O_Atom_'''+str(i))
  if OCoreShift[i] > 0.0:
   Target=np.round(OCoreShift[i], decimals=2)
   k.write('''
awk '{print $1-'''+str(abs(Target))+''', $2/'''+str(abs(OCoreEn[i]))+'''}' XAS_eV_wrong.dat > XAS_shift.dat
cd -
   ''')
  else:
   Target=np.round(OCoreShift[i], decimals=2)
   k.write('''
awk '{print $1+'''+str(abs(Target))+''', $2/'''+str(abs(OCoreEn[i]))+'''}' XAS_eV_wrong.dat > XAS_shift.dat
cd -
   ''')

if DONitrogen==True:
 for i in range(NumN):
  k.write('''
cd ../../Single_Nitrogen/N_Atom_'''+str(i))
  if NCoreShift[i] > 0.0:
   Target=np.round(NCoreShift[i], decimals=2)
   k.write('''
awk '{print $1-'''+str(abs(Target))+''', $2/'''+str(abs(NCoreEn[i]))+'''}' XAS_eV_wrong.dat > XAS_shift.dat
cd -
   ''')
  else:
   Target=np.round(NCoreShift[i], decimals=2)
   k.write('''
awk '{print $1+'''+str(abs(Target))+''', $2/'''+str(abs(NCoreEn[i]))+'''}' XAS_eV_wrong.dat > XAS_shift.dat
cd -
   ''')

if DOCarbon==True and DOOxygen==False and DONitrogen==False:
 k.write("""
sed -i 's/xxxaxxx/"""+str(NumC)+"""/g' SumSpec.py
 """)
if DOCarbon==True and DOOxygen==True and DONitrogen==False:
 k.write("""
sed -i 's/xxxaxxx/"""+str(NumC+NumO)+"""/g' SumSpec.py
 """)
if DOCarbon==True and DOOxygen==False and DONitrogen==True:
 k.write("""
sed -i 's/xxxaxxx/"""+str(NumC+NumN)+"""/g' SumSpec.py
 """)
if DOCarbon==True and DOOxygen==True and DONitrogen==True:
 k.write("""
sed -i 's/xxxaxxx/"""+str(NumC+NumO+NumN)+"""/g' SumSpec.py
 """)

k.close()
