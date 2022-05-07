#!/usr/bin/tclsh

set refSel "chain D and backbone"
#set wrapSel "protein"

proc main {} {

global argv argc
global refSel wrapSel

set input_psf		[lindex $argv 0]
set ref_pdb		[lindex $argv 1]
set input_dcd		[lindex $argv 2]
set outName		[lindex $argv 3]

mol new $input_psf
mol addfile $ref_pdb
mol addfile $input_dcd waitfor all

set ref [atomselect top "$refSel" frame 0]
set number [molinfo top get numframes]
for {set i 1} {$i <= $number} {incr i} {
set sel [atomselect top "$refSel" frame $i]
set all [atomselect top all frame $i]
$all move [measure fit $sel $ref]
}

animate write dcd $outName beg 1

}

main
exit

