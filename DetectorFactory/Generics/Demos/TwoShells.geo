#
# Ascii file that described a collection of perfect shells for HermeShellDetectorConstruction
# You must write from the smallest radius to the largest !!
# Two shells are mandatory one named Inner and a second one named Outer
# In case you want to work with only one active shell the second is composed of air and set inative
# (see Example 1)
#  
# Full description of one shell per line with the following conventions:
#   (lengths are in centimeter and angles in degres)
#   (is_active is equal to 1 if it is an active shell 0 otherwise. active means a detector-like shell)
#
# name_of-the-shell  material rMin rMax starting_phi delta_phi starting_theta delta_theta is_active
#
# Example 1 - one unique shell (the second one is then composed of air and inactive)
#	Shell:0	NaI 	10.	15.	0.	360.	0.	180.	1
#	Shell:1	Air 	25.	40.	0.	360.	0.	180.	0
#
# Example 2 - two active shells inner of NaI and outer BGO
#	Shell:0	NaI 	10.	15.	0.	360.	0.	180.	1
#	Shell:1	BGO 	25.	40.	0.	360.	0.	180.	1
#
# Example 3 - two active shells, the target chamber and a dead layer between the inner and outer
#	Target	Aluminium	8.	8.5	0.	360.	0.	180.	0
#	Shell:0	NaI 	10.	15.	0.	360.	0.	180.	1
#	Abs1	Lead	16.0	16.05	0.	360.	0.	180.	0	
#	Shell:1	BGO 	25.	40.	0.	360.	0.	180.	1
#
Shell:0     LaBr3 	10.0	15.	0.	360. 0.	180.	1
Shell:1     NaI 	15.1	35.	0.	360. 0.	180.	1
#
#
#