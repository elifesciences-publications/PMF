
# Following scripts have been used for generating the data represented in Figures 2C and S10. Note that the example scripts below is used for a subset of our umberalla sampling windows. An identical set of scripts have been used for the rest of the windows but on different *.psf files.


# Calculating local thickness around the channel protein

cd ~/PMF-final/0A

mol delete all
resetpsf
set j {17 18 19 20 21 22}
set output [open thcknss_cm_locChan.txt w]

foreach p $j {


	
mol load psf pos0A.psf dcd PMF-${p}A.dcd 

set nf [molinfo top get numframes]
for {set i 100} {$i<= $nf} {incr i } {
set upp [atomselect top "name C1 and z<-12 and within 8 of protein" frame $i]
set down [atomselect top "name C1 and z>10 and within 8 of protein" frame $i]


set zmax [lindex [measure center $upp] 2]
set zmin [lindex [measure center $down] 2]
set thickness [expr $zmax - $zmin]


puts $output "$i $thickness "

}
puts $output " \n \n "
}
close $output 

mol delete all
resetpsf



# To calculate the  hydrophobic length of the protein



cd ~/PMF-final/0A

mol delete all
resetpsf
set j {17 18 19 20 21 22}

set output [open 0A-thcknss-Protein_all.txt w]
foreach p $j {

mol delete all
resetpsf
	
mol load psf pos0A.psf dcd PMF-${p}A.dcd 

set nf [molinfo top get numframes]
for {set i 100} {$i<= $nf} {incr i } {
set upp [atomselect top "protein and resid 16" frame $i]
set down [atomselect top "protein and resid 48" frame $i]


set zmax [lindex [measure center $upp] 2]
set zmin [lindex [measure center $down] 2]
set thickness [expr $zmin - $zmax]


puts $output "$i $thickness "

}
puts $output " \n \n "

}
close $output 

mol delete all
resetpsf



# To calculate the local curvature of the bilayer


cd ~/PMF-final/0A

mol delete all
resetpsf
set j {17 18 19 20 21 22}	
foreach p $j {

mol load psf pos0A.psf dcd PMF-${p}A.dcd 

set output [open R-${p}-final.dat w]
set nf [molinfo top get numframes]
for {set i 100} {$i<= 200} {incr i } {
set uppR [atomselect top "name P and z<-12 and (x>65 or x<-65)" frame $i]
set uppL [atomselect top "name P and z<-12 and within 8 of protein" frame $i]

set zmax [lindex [measure center $uppR] 2]
set zmin [lindex [measure center $uppL] 2]
set height [expr $zmax - $zmin]

set xmin [lindex [measure minmax $uppR] 0 0]
set xmax [lindex [measure minmax $uppL] 1 0]

set length [expr 2*($xmax - $xmin)]


set radius [expr $height/2 + $length*$length/(8*$height)]

puts $output "$i $length $height $radius"
}
close $output
mol delete all
resetpsf
puts $output " \n \n "
}





# To calculate the global lipid bilayer thickness


cd ~/PMF-final/0A

mol delete all
resetpsf
set j {17 18 19 20 21 22}	
set output [open 0a_global_mem-thcknss.txt w]

foreach p $j {

	
mol load psf pos0A.psf dcd PMF-${p}A.dcd 

set nf [molinfo top get numframes]
for {set i 100} {$i<= $nf} {incr i } {
set upp [atomselect top "name P and z<-12" frame $i]
set down [atomselect top "name P and z>10" frame $i]


set zmax [lindex [measure center $upp] 2]
set zmin [lindex [measure center $down] 2]
set thickness [expr $zmax - $zmin]


puts $output "$i $thickness "

}
puts $output " \n \n "

}
close $output 



----------------- Analysing the PMF trajectories for PMF calculations 

awk '($1>=0000000 && $1 <= 25000000 && $1%500==0) {print $1, $2}' ./PMF-33A.colvars.traj > ./PMF-33A-colv.txt

---------------------------- Using WHAM for free energy calculations

./wham/wham/wham P=0 16 32 80 0.001 298 0 metafile1 PMF-MscS-3-2019.txt > PMF-MscS-3-2019.log



