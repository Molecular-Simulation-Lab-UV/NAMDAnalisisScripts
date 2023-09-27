#Input Parameters
##Ideally run with a filtered trajectory, as this code load the full trajectory in memory (modificar linea 3 & 4)###
set input_pdb  eq5_prot_mem.pdb
set input_dcd  eq5_prot_mem_stride.dcd
set RMSFselection "protein and chain D and name CA and resid 93 to 230"
set RMSDslection "protein and chain D and name CA and resid 93 to 230"
set FITselection "protein and chain D and name CA and resid 93 to 230"
set CenterSelection "resname LIG and name C8"
set frames_begin 0
#Axix alignment parameters
set selVecA "resname LIG and name C13"
set selVecB "resname LIG and name C8"
set selLig "resname LIG"
set extraR 0
#Set the flag for RMSD matrix calculation (heavy calculation)
set CalcMAtrix 1
set CentWriteTraj 0

 ## END INPUT ##

#PROCEDURES

# Calculates rmsf
proc rmsf { {mol top} {outfile rmsf.dat} {first_at 1} selection } {
    set sel [atomselect top "$selection"]
    set out [open $outfile w]
    set rmsf_list [measure rmsf $sel]
    set name_list [$sel get name]
    set i $first_at
    puts $out [format {%10s %10s %10s} "#atom Name" "atom Index" "rmsf"]
    foreach rmsf $rmsf_list name $name_list {
	puts $out [format {%10s %10d %10.2f} $name $i $rmsf]
	incr i
    }         
    close $out 
}

# stuctural fit
proc fit { {mol top} selection refFrame  } {
    set ref [atomselect $mol "$selection" frame $refFrame]
    set sel [atomselect top "$selection"]
    set all [atomselect top all]
    set nf [molinfo $mol get numframes]
    for {set frame 0} {$frame < $nf} {incr frame} {
	$sel frame $frame
	$all frame $frame 
	$all move [measure fit $sel $ref]
    }
} 

# Computes a RMSD matrix
proc RMSD_MAT {frames_begin numFrame NumframesTitle  sel1 sel2} {

    set outfile1 [open rmsd_matx.dat  w]
    set outfile2 [open rmsd_hist.dat w]
    puts $outfile1  "TITLE"
    puts $outfile1 "        rmsd-matrix for $NumframesTitle =  $numFrame structures "
    puts $outfile1 "END"
    puts $outfile1 "RMSDMAT"
    puts $outfile1 "# number of frames  skip  stride"
    puts $outfile1 "$numFrame      0       1"
    puts $outfile1 "#precision"
    puts $outfile1 "10000"

    for { set f1 $frames_begin } { $f1 < [expr $numFrame -1] } { incr f1 } {
	$sel1 frame $f1
	for { set f2 [expr $f1 +1] } { $f2 < $numFrame } { incr f2 } {
	    $sel2 frame $f2
            set rmsd [measure rmsd $sel1 $sel2]
            puts $outfile2 "$rmsd"
            set rmsd [expr int($rmsd*1000)]
            set ref [string range [format " %7i" $f1] end-7 end]
            set str [string range [format "%7i" $f2] end-7 end]
            set rms_f [string range [format "%7i" $rmsd] end-7 end]
            puts $outfile1 "$ref $str $rms_f"
        }
    }
    puts $outfile1 "END"
    close $outfile1
    close $outfile2
}

#center
proc center {selText1 selText2 } {
    #Molecules
    set sel1 [atomselect top "$selText1"  frame 0]
    #Center
    set sel2 [atomselect top "$selText2"  frame 0]
    $sel1 moveby [vecinvert [measure center $sel2]]
}
# requires slight modifications based on molecules
proc Align {selVecA selVecB  selLig {extraR 1} } {
    
    set A   [atomselect top "$selVecA" frame 0]
    set B   [atomselect top "$selVecB" frame 0]
     
    set sel [atomselect top "$selLig"  frame 0]

    set vecA "[$A get x] [$A get y] [$A get z]"
    set vecB "[$B get x] [$B get y] [$B get z]"
  
    set DistAB [vecsub $vecA $vecB]
    #Aling vector DistAB in the x axis
    set M [transvecinv $DistAB]
    #move sel along vector aligned in x 
    $sel move $M
    #Now align in z
    set M [transaxis y -90] 
    $sel move $M        
    #Extra rotations 
    if {$extraR} {
	set M [transabout {0 0 1} $extraR]
	$sel move $M
     
     
	#set M [transabout {0 1 0} 40]
	#$sel move $M
    }
}

#MAIN#
proc main {} {
    global input_pdb input_dcd RMSFselection FITselection RMSDslection frame frames_begin CalcMAtrix CentWriteTraj
    global selVecA selVecB  selLig CenterSelection extraR
    #loads psf and dcd
    mol load pdb $input_pdb
    mol addfile $input_dcd waitfor all
 
    
    #center $selLig $CenterSelection
    #Align $selVecA $selVecB  $selLig $extraR
    # structural fit agasint initial frame
    fit  top $FITselection 0
    #Write centered trajectory
    #if {$CentWriteTraj} {
#	set all [atomselect top all]
#	animate write dcd $input_dcd  $all
 #   }

    #Computes RMSF
    rmsf top rmsf.dat 1  $RMSFselection
    
    #sets selections for RMSD matrix calculation
    if {$CalcMAtrix} {
	set sel1 [atomselect top "$RMSDslection"]
	set sel2 [atomselect top "$RMSDslection"]
	set numFrame [molinfo top get numframes]
	set NumframesTitle [expr $numFrame -1]
	# Computes rmsd matrix
	RMSD_MAT $frames_begin $numFrame $NumframesTitle  $sel1 $sel2
    }
}


###END MAIN##

main 
exit
