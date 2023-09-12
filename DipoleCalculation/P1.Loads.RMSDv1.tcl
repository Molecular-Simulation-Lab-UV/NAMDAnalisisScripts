#!/usr/bin/tclsh

###### Change here ###########################################################
set refSel  "chain A and resid 78 162 189 204"; #reference need to be centered FAPIP REFERENCE
set Wrapped 0;# is the trajectory wrapped,  0 no 1 yes
set wrapSel "chain A"              
set outName "LOADS.RMSDA.dat"
#set WatSel "name OH2"
#set WatSel  "chain A and resid 78 162 189 204 and name CB";
set WatSel "chain A and resid 78 162 189 204 and (name NE2 or name NH2 or name CG2)"
source /opt/vmd-1.9.3/bin/bigdcd.tcl
set firstDCD 10
set lastDCD 17

#### END #####################################################################

############# binning and pore definitions ################################
for {set i $firstDCD} {$i <= $lastDCD} {incr i} {
	set dcd($i) /home/nespinoza/NAS-mnt/jgarate/Aquaporins/FaPIP/FaPIP2.1.one.template.2b5f/MD/MD1Corrected/eq$i.dcd
}

set minZ   -20; #Min Z value for pore
set maxZ    20; #Max Z value for pore
set binNumZ 80; #Numbers of bins along axial Dim
set BinSizeZ  [expr 1.000*($maxZ-$minZ)/$binNumZ]
set rad      8; #Radius of pore
set rad2  [expr $rad*$rad];
for { set i 0}  {$i < $binNumZ} {incr i} {
    set P1Axis($i)  0
    set LoadsCounter($i) 0
}
set steps 0; # counter for averages
################################################################################
proc BinWatersDipZ { &arrName1 &arrName2 } {
    global minZ maxZ rad rad2;#cylinder parameters
    global BinSizeZ binNumZ;  #Binning parameters
    global WatSel;            #water selection 
    upvar 1 ${&arrName1}  LoadsCounter;# Total Counter of observations binned in z
    upvar 1 ${&arrName2}  P1Axis;
    set WatLoadCounter 0
    #Get indexes of  all waterr within cylinder
    set indexes [ [atomselect top "$WatSel and z> $minZ and z < $maxZ and (x^2 + y^2)< $rad2"] get index]
    foreach index $indexes {
	# Do selections and collect vectors
	set indexH1 [expr $index +1]
	set indexH2 [expr $index +2]
	set sel [atomselect top "index $index $indexH1 $indexH2"]
	set coord [$sel get {x y z}]
	set x [lindex [lindex $coord 0] 0]
	set y [lindex [lindex $coord 0] 1]
	set z [lindex [lindex $coord 0] 2]
	# Binning
	set stepZ [expr int(( ($z-$minZ)/($BinSizeZ) ))]
	# Accumulate Axial Radial histogram
	if { $stepZ < $binNumZ && $stepZ >= 0 } {
	    #Accumulate P1 values, binned along z axis
	    set vectorDip [vecnorm [measure dipole $sel -masscenter]]
	    set P1Axis($stepZ) [expr $P1Axis($stepZ) +[lindex $vectorDip 2]]
	    incr LoadsCounter($stepZ); # and accumulate for averages
	}
	$sel delete
    }
    return [llength $indexes];
}

#Ubin Average dipole along axial dim
proc UnbinDipZ { &arrName1 &arrName2 } {
    upvar 1 ${&arrName1} P1Axis
    upvar 1 ${&arrName2} LoadsCounter
    global binNumZ BinSizeZ minZ
    set out [open "P1_axial2A.dat" w]
    puts $out "#AvgP1 along z axis"
    puts $out "#z          <P1>"
    for { set i 0}  {$i < $binNumZ} {incr i} { 
	set z [format {%8.2f} [expr $i*$BinSizeZ + $minZ + $BinSizeZ*0.5]]
	if { $LoadsCounter($i) > 0} {
	    set AverageP1 [format {%8.2f} [expr 1.00*$P1Axis($i)/$LoadsCounter($i)] ]
	} else {
	    set AverageP1 [format {%8.2f} 0 ]
	}
	puts -nonewline $out $z
	puts -nonewline $out $AverageP1
	puts $out ""
    }
    close $out
}
#Ubin Average Load along axial dim
proc UnbinLoadZ {&arrName } {
    upvar 1 ${&arrName} LoadsCounter
    global minZ BinSizeZ binNumZ steps
    set out [open "Loads_axial2A.dat" w]
    puts $out "#AvgLoad along z axis"
    puts $out "#z           <Load>"
    for { set i 0}  {$i < $binNumZ} {incr i} { 
	set z [format {%8.2f} [expr $i*$BinSizeZ + $minZ + $BinSizeZ*0.5]]
	set AverageLoad [format {%8.2f} [expr 1.00*$LoadsCounter($i)/$steps] ]
	puts -nonewline $out $z
	puts -nonewline $out $AverageLoad
	puts $out ""
    }
    close $out
}

#RMSD agains ref strcuture
proc rmsd_bigdcd { {mol top} selection } {
    set ref [atomselect $mol "$selection" frame 0]
    set sel [atomselect $mol "$selection"]
    set all [atomselect $mol all]
    $all move [measure fit $sel $ref]
    set rmsd [measure rmsd $sel $ref]
    return $rmsd
 }
#Reads each frame at a time and collects data
proc ts {frame} {
#    global outName
    global wrapSel refSel Wrapped
    global LoadsCounter P1Axis steps
    if {$Wrapped == 0} {
	pbc wrap -center com -centersel "$wrapSel" -compound chain -now
    }
    set RMSD  [rmsd_bigdcd top $refSel];# Always before any calculation
    set LOADS [BinWatersDipZ LoadsCounter P1Axis]
    set out   [open "LOADS.RMSD2A.dat" a+]
    puts $out [format "%10d %10.3f %10.3f" $frame $LOADS $RMSD]
    close $out
    incr steps
}
# Reads arguments and call ts via bigdcd
proc main {&arrName} {
    global argc argv
	upvar 1 ${&arrName} dcd
    global outName argc argv
    global firstDCD lastDCD
    package require pbctools
#     if { $argc != 3 } {
#	 puts "The P1.loads.tcl script requires 3 variables  to be inputed."
#	 puts "For example, vmd -dispdev text -e path/to/program -args path/to/inputpsf  path/to/inputdcd path/to/refdcd"
#	 puts "Please try again."
#	 exit        
#     }

    set input_psf     [lindex $argv 0]
#    set input_dcd     [lindex $argv 1]
    set ref_dcd       [lindex $argv 1]

    
    mol load psf $input_psf
    mol addfile  $ref_dcd
    set out2 [open "LOADS.RMSD2A.dat" w]
    puts $out2 "#Frame      LOADS      RMSDsel"
    close $out2;
	for {set i $firstDCD} {$i <= $lastDCD} { incr i } {
#    bigdcd ts $input_dcd
    	bigdcd ts $dcd($i)
		bigdcd_wait
	
	}
}
###########################################
main dcd
UnbinDipZ P1Axis LoadsCounter
UnbinLoadZ LoadsCounter
exit
