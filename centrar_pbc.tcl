# ================================================================
# SCRIPT FINAL REPARADO - TETRÁMEROS + LOOPS (Z-axis)
# ================================================================

mol delete all
mol new BvPIP2.1.hmr.psf
mol addfile BvPIP2.1.hmr.pdb

puts "Cargando trayectoria DCD..."
mol addfile 24DynCm.dcd waitfor all

package require pbctools

puts "Leyendo XST..."
pbc readxst 24DynCm_sync.xst

set n_dcd [molinfo top get numframes]
set prot_sel "protein and segname PROA PROB PROC PROD"

# 1. UNWRAP (Continuidad temporal)
puts "Unwrapping temporal..."
pbc unwrap -sel $prot_sel -all

# 2. JOIN ESPACIAL (Unión de monómeros)
puts "Join de segmentos..."
pbc join seg -sel $prot_sel -all

# --- LÍNEA AÑADIDA PARA REPARAR LOOPS EN EJE Z ---
puts "Reparando loops y residuos cortados..."
pbc join residue -sel $prot_sel -all
# ------------------------------------------------

# 3. WRAP FINAL (Centrado en la proteína)
puts "Wrapping final..."
pbc wrap \
    -centersel $prot_sel \
    -center com \
    -compound seg \
    -all

# GUARDADO
set all [atomselect top all]
set outname "24DynCm_FINAL_TETRAMER_OK_PRUEBA.dcd"

puts "Guardando DCD final..."
animate write dcd $outname sel $all beg 0 end [expr {$n_dcd - 1}] waitfor all

puts "========================================================="
puts "PROCESO COMPLETADO - LOOPS CORREGIDOS"
puts "========================================================="
