#################################################################################
#										#
#	This Script collects the original orbitals and localized orbitals	#
#	and then runs the ShiftXas.py Script, which then generates two		#
#	additional .sh-Scripts.							#
#										#
#	The XASSpec.sh Script then generates all the XAS_eV.dat files		#
#	that have the "wrong" position. Afterwards, the Shift.sh Script		#
#	is called to shift these XAS_eV.dat towards the correct position	#
#	with respect to the original energies.					#
#										#
#	If the Orbitals have been localized well enough in the beginning	#
#	(i.e. at least a coefficient of 0.98), the energy of the orbital	#
#	is kept. Otherwise, the orbital will be shifted towards the mean	#
#	value of the delocalized orbitals.					#
#										#
#################################################################################

cp ../../Localize/*.out ./Loc.out
cp ../../Localize/*.out ./MOs.out
cp ../../MOs/*.out ./MOs.out
python2.7 ShiftXas.py
