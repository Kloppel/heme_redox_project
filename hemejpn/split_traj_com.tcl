# load trajectory into vmd

#set allexes [glob -directory "" -- "database/pdb/from_rcsb/*.pdb"]
set allexes [glob -directory "" -- "*.pdb"]
foreach f $allexes {
    puts "$f"

    set name $f
    echo $name 123123
    #set name 1f1f.pdb
    mol load  pdb $name
    set num_frames [molinfo top get numframes]

    for {set i 0} {$i<=$num_frames-1} {incr i} {

        set listName {"HEM" "HEA" "HEC" "HED"}
        set number 0
        for {set j 0} {$j<=3 && $number==0} {incr j} {
            set heme [lindex $listName  $j]
            set all [atomselect top "(resname $heme)" frame $i]
            set number  [$all num ]
            echo $heme $number
            echo j $j
        }
        echo 1
        for {set i 0} {$i<=0 && $number==0} {incr i} {
            echo $name failed
        }
        echo 2
        #set heme "HEC"
        #set all [atomselect top "(resname $heme)" frame $i]
        $all writepdb prepared/$name
        echo 3
    }
    mol delete all
}
quit
# vmd -dispdev text -e split_traj_com.tcl 