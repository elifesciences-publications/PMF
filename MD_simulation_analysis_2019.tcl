
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




---------------------------------------




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




# To calculate the global area of the lipid bilayer

cd ~/PMF-final/0A

mol delete all
resetpsf
set j {17 18 19 20 21 22}
set output [open 0a_global_mem-area.txt w]

foreach p $j {


mol load psf pos0A.psf dcd PMF-${p}A.dcd

set nf [molinfo top get numframes]
for {set i 180} {$i<= $nf} {incr i } {
set upp [atomselect top "name P and z<-12" frame $i]
set down [atomselect top "name P and z>10" frame $i]




set xmin [lindex [measure minmax $upp] 0 0]
set xmax [lindex [measure minmax $upp] 1 0]

set Xlength [expr $xmax - $xmin]

set ymin [lindex [measure minmax $upp] 0 1]
set ymax [lindex [measure minmax $upp] 1 1]

set Ylength [expr $ymax - $ymin]

set area [expr $Xlength * $Ylength]
puts $output "$i $area "

}
puts $output " \n \n "

}
close $output




# To calculate curvature across the bilayer as a function of distance from the protein center.

cd ~/..

mol delete all

resetpsf


mol delete all
resetpsf
set s {17 18 19 20 21 22} 
foreach p $s {

mol load psf pos0A.psf dcd lessFramesPMF-${p}A.dcd 



set output [open curvature-${p}A.dat w]

set nf [molinfo top get numframes]





# define the appropriate bin size for integration



set no_bin 15

	


for {set i 1} {$i<= $no_bin} {incr i} {
set AveHeight($i) 0
}



for {set i 1} {$i< $nf} {incr i } {

set LBond 30
set UBond 35


set ProtCm [atomselect top "protein and resid 88 to 105" frame $i]
set XProt [lindex [measure center $ProtCm] 0]
set YProt [lindex [measure center $ProtCm] 1]


		for {set j 1} {$j<= $no_bin} {incr j } {

set LBondp [expr $LBond+1]
set LBondm [expr $LBond-1]

set UBondp [expr $UBond+1]
set UBondm [expr $UBond-1]


			set LowR [atomselect top "name P and (z<-12) and (sqrt( (x-$XProt)^2 + (y-$YProt)^2 ) - $LBondp <0) and (sqrt((x-$XProt)^2 + (y-$YProt)^2) - $LBondm >0)" frame $i]

			set HighR [atomselect top "name P and (z<-12) and (sqrt( (x-$XProt)^2 + (y-$YProt)^2 ) - $UBondp <0) and (sqrt((x-$XProt)^2 + (y-$YProt)^2) - $UBondm >0)" frame $i]



			set zmin [lindex [measure center $LowR] 2]

			set zmax [lindex [measure center $HighR] 2]

			



		

			

			

			set height [expr $zmax - $zmin]

		    set delta_r [expr 30/$no_bin]



			set mean_Radius_bin [expr $LBond + ($delta_r/2)]

			set Area_bin [expr 2*3.14*$mean_Radius_bin*$delta_r]

			

			
			puts $output " Frame: $i  Bin: $j  LBond: $LBond UBond: $UBond Height: $height Delta_R: $delta_r  Mean Radius: $mean_Radius_bin  Mean area: $Area_bin"



			set AveHeight($j) [expr ($AveHeight($j) + $height)]
            set AveLength($j) [expr ($AveLength($j) + $mean_Radius_bin)]
			 
			set LBond [expr $UBond]

			set UBond [expr $UBond+ $delta_r]


			

		}



puts $output \n



#set radius [expr $height/2 + $length*$length/(8*$height)]



}

 puts $output \n\n\n\n

puts $output " average over all  $i frames: "

for {set n 1} {$n<= $no_bin} {incr n } {

puts $output " height($n): [expr $AveHeight($n) / ($nf-1)]"

}




close $output

}
close $output


----------------- Analysing the PMF trajectories for PMF calculations 

awk '($1>=0000000 && $1 <= 25000000 && $1%500==0) {print $1, $2}' ./PMF-33A.colvars.traj > ./PMF-33A-colv.txt

---------------------------- Using WHAM for free energy calculations

./wham/wham/wham P=0 16 32 80 0.001 298 0 metafile1 PMF-MscS-3-2019.txt > PMF-MscS-3-2019.log



