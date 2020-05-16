mol delrep top top
for {set j 0} {$j<10000} {incr j} {
mol representation QuickSurf 0.7 0.3 0.5 3
set temp [atomselect top "resid = $j"]
if {[$temp num] != 0} {
if {$j > 5000} {
mol color colorid 0
} else {
  if {$j > 3000} {
  mol color colorid 7
  } else {
    if {$j > 2000} {
    mol color colorid 27	} else {
	mol color colorid 1
	}
    }
  }
mol selection "resid $j"
mol addrep top
}
}
