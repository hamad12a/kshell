#!/bin/sh
# export OMP_STACKSIZE=1g
export GFORTRAN_UNBUFFERED_PRECONNECTED=y
ulimit -s unlimited
# export OMP_NUM_THREADS=1
# export KMP_NUM_THREADS=20
export OMP_DYNAMIC='true'

rm Cr48.ptn 
echo "0" | ../bin/gen_partition.py kb3.snt Cr48.ptn 4 4 1

cat > cr48.input <<EOF
&input
  fn_int = "kb3.snt"
  fn_ptn = "Cr48.ptn"
!  fn_ptn = "Cr48t4.ptn"
  mtot = 0
  n_eigen = 1 ! 5
  n_restart_vec = 10
  max_lanc_vec = 200
  maxiter = 1000
  hw_type = 1
  is_double_j = .false.
!  is_double_j = .true. 
!  fn_ptn_init = "Cr48t4.ptn"
!  fn_load_wave = "init.wav"
  fn_save_wave = "cr48j0.wav"
! is_load_snapshot = .true.
  eff_charge=1.5, 0.5
! time_limit_hour = 0.01
  mode_lv_hdd = 0 ! 0
! is_calc_tbme = .true.
! ss_e_range = -33.5d0, -30.d0
! ss_e_range = -33.d0, -32.8d0
!n_block =  2
&end
EOF

# cp ../bin/kshell.exe . 
./kshell.exe cr48.input 
# ../src_jump/kshell cr48.input 

# rm tmp_snapshot_Cr48.ptn*

