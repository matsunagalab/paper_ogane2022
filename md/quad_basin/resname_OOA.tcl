# changes GXXX resnames in KBGO pdb file to real resnames
# input: aa pdb file (1 chain)
# input: kbgo pdb file
# output: kbgo pdb file with GXX resnames switched by the real resnames
###############################################################################################################
set PDB_NAME_AA AA_OOA.pdb
set PDB_NAME_KBGO GO_OOA.pdb
#1wdn_A.pdb
#GO_1wdna.pdb
###############################################################################################################
# i/o
#####

# displays a variable with its name
# input: var_name
# output: var_name, var value
proc proc_var_disp {var_name} {
	set call [dict get [info frame -1] cmd]
	puts "[lindex $call 1]\t$var_name"
}

# txt files
###########

# reads a txt file into a list
# input: file with 1-col entries
# output: list with each line as an element
proc proc_file_read_to_list {file_name} {
	set file_in [open $file_name r]
	set list_out {}
	while {[gets $file_in line]>=0} {
		lappend list_out $line
	}
	close $file_in
	return $list_out
}

# read col-n from a file
# input: file_name
# input: col_n (from 0)
# output: list with all elements of col_n
proc proc_file_read_col_n {file_name col_n} {
	set file_in [open $file_name r]
	set list_col {}
	while {[gets $file_in line]>=0} {
		lappend list_col [lindex $line $col_n]
	}
	close $file_in
	return $list_col
}

# inserts a string into a txt file name
# input: file name XXX.YYY
# input: str_add
# output: XXX_str_add.YYY
proc proc_file_name_insert_str {file_name str_add} {
	set list_file_name [split $file_name .]
	set file_name_add [lindex $list_file_name 0]_$str_add.[lindex $list_file_name 1]
	return $file_name_add
}

# return file name without a path
# input: file name with path (../../XXX/../filename.txt)
# output: name of last file in path (filename.txt)
proc proc_file_name_no_path {str_path} {
	set list_str_path [split $str_path /]
	return [lindex $list_str_path end]
}

# returns a string without its ending
# input: XXX.YYY
# output: XXX
proc proc_str_no_end {str_full} {
	set list_str [split $str_full .]
	set str_no_end [lindex $list_str 0]
	return $str_no_end
}

# creating output files

# creating output file name from script file name and a given input file name
# input: input_name.xxx
# ouptut: SCRIPT_NAME_input_name.out
proc proc_file_out_name_script_param1 {input_name} {
	global SCRIPT_NAME
	set input_name_noend [proc_str_no_end $input_name]
	set file_out_name "${SCRIPT_NAME}_${input_name_noend}.out"
	return $file_out_name
}

# creating output file name from script file name and 2 given input file names
# input: input_name_1.xxx, input_name_2.xxx
# ouptut: SCRIPT_NAME_input_name.out
proc proc_file_out_name_script_param2 {input_name_1 input_name_2} {
	global SCRIPT_NAME
	set input_name_1_noend [proc_str_no_end $input_name_1]
	set input_name_2_noend [proc_str_no_end $input_name_2]
	set file_out_name "${SCRIPT_NAME}_${input_name_1_noend}_${input_name_2_noend}.out"
	return $file_out_name
}

# manipulating numbers
######################

# rounds a number to the desired decimal digits
# input: num_tar (eg 1.2345676)
# input: num_digit (eg 2)
# output: number rounded to num_digit decimal digits (1.23)
proc proc_num_round_dec {num_tar num_digit} {
	set num_pow_10 [expr 10**$num_digit]
	set num_round [expr double(round($num_pow_10*$num_tar))/$num_pow_10]
	return $num_round
}

# manipulating pdb/psf/dcd
##########################

# open a pdb file
# input: pdb_name
# output: molecule
proc proc_pdb_open {pdb_name} {
	set mol_pdb [mol new $pdb_name]
	return $mol_pdb
}

# open a dcd traj
# input: pdb, psf, dcd
# output: a molecule element
proc proc_dcd_open {pdb_name psf_name dcd_name} {
	set mol_dcd [mol new $psf_name]
	mol addfile $pdb_name $mol_dcd
	mol addfile $dcd_name waitfor all $mol_dcd
	return $mol_dcd
}

# lists
#######

# list edit

# display a list, line per element
# input: list_tar
# output (screen): elements of list, line per element
proc proc_puts_list_el_per_line {list_tar} {
	foreach el $list_tar {
		puts $el
	}
}

# writes a list to a file, line per element
proc proc_write_list_el_per_line {list_tar file_out} {
	foreach el $list_tar {
		puts $file_out $el
	}
}


# remove all empty elements from list
# input: list_tar
# output: list_tar without any empty elements
proc proc_list_clean_empty {list_tar} {
	set list_clean {}
	foreach el $list_tar {
		if {$el!={} && $el!=""} {
			lappend list_clean $el
		}
	}
	return $list_clean
}

# inverts a 2D list dimension
# input: list_list_dim1_dim2
# output: list_list_dim2_dim1
proc proc_list_invert_dimension {list_list_dim1_dim2} {
	set list_list_dim2_dim1 {}
	set list_tmp [lindex $list_list_dim1_dim2 0]
	set num_dim2 [llength $list_tmp]
	for {set ind 0} {$ind<$num_dim2} {incr ind} {
		set list_dim2 {}
		foreach list_dim1 $list_list_dim1_dim2 {
			lappend list_dim2 [lindex $list_dim1 $ind]
		}
		lappend list_list_dim2_dim1 $list_dim2
	}
	return $list_list_dim2_dim1
}

# creating lists
 
# create a list of all integers between 2 integers, not including 
# input: int_begin, int_end
# output: list of integers between int_begin to int_end, not including the end integers
proc proc_get_list_int_between_beg_end_no_inc {int_begin int_end} {
	set list_out {}
	for {set int_now [expr $int_begin+1]} {$int_now<$int_end} {incr int_now} {
		lappend list_out $int_now
	}
	return $list_out
} 

# create a list of all integers between 2 integers, including ends
# input: int_begin, int_end
# output: list of integers between int_begin to int_end, including the end integers
proc proc_get_list_int_between_beg_end_yes_inc {int_begin int_end} {
	set list_out {}
	for {set int_now $int_begin} {$int_now<=$int_end} {incr int_now} {
		lappend list_out $int_now
	}
	return $list_out
} 

# calculations on lists

# finds max val in list
# input: list_tar 
# output: max val in list_tar 
proc proc_list_get_max {list_tar} {
	set list_sort [lsort -real $list_tar]
	set val_max [lindex $list_sort end]
	return $val_max
} 

# finds min val in list
# input: list_tar 
# output: min val in list_tar 
proc proc_list_get_min {list_tar} {
	set list_sort [lsort -real $list_tar]
	set val_min [lindex $list_sort 0]
	return $val_min
} 

# GO model related
##################

# turns a KBGO resname into a resid
# input: KBGO resname (e.g. G22)
# output: only resid (22)
proc proc_go_resname_to_resid {resname_go} {
	set resid_go [string trimleft $resname_go "G"]
	return $resid_go
}

###############################################################################################################
# open pdb
set mol_kbgo [mol new $PDB_NAME_KBGO]
set mol_aa [mol new $PDB_NAME_AA]

# change resname
set sel_kbgo_ca [atomselect $mol_kbgo "name CA"]
set sel_aa_ca [atomselect $mol_aa "name CA"]
puts [$sel_kbgo_ca num]
puts [$sel_aa_ca num]
set i 1
foreach resname_aa [$sel_aa_ca get resname] resid_kbgo [$sel_kbgo_ca get resid] {
	set sel_kbgo_sing_res [atomselect $mol_kbgo "resid $resid_kbgo"]
  puts "$resid_kbgo\t$resname_aa"
  $sel_kbgo_sing_res set resname $resname_aa
  $sel_kbgo_sing_res set resid $i
  $sel_kbgo_sing_res delete
  incr i 1  
}
$sel_kbgo_ca delete
$sel_aa_ca delete

# write pdb
set pdb_new_name [proc_file_name_insert_str $PDB_NAME_KBGO resname]
set sel_kbgo_all [atomselect $mol_kbgo "all"]
$sel_kbgo_all writepdb $pdb_new_name
$sel_kbgo_all delete

mol delete $mol_kbgo
mol delete $mol_aa
exit

