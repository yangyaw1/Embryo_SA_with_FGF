mol delrep top top
for {set j 0} {$j<10000} {incr j} {
mol representation QuickSurf 0.7 0.3 0.5 3
set temp [atomselect top "resid = $j"]
if {[$temp num] != 0} {
if { ($j >= 5600 && $j < 6000) || ($j >= 2600 && $j < 3000) || ($j >= 1600 && $j < 2000)} {
mol color colorid 3
} else {
  if {($j >= 5400 && $j < 5600) || ($j >= 2400 && $j < 2600) || ($j >= 1400 && $j < 1600) } {
  mol color colorid 15
  } else {
    if {($j >= 5200 && $j < 5400) || ($j >= 2200 && $j < 2400) || ($j >= 1200 && $j < 1400) } {
    mol color colorid 10
    } else {
    if {($j >= 5000 && $j < 5200) || ($j >= 2000 && $j < 2200) || ($j >= 1000 && $j < 1200) } {
    mol color colorid 1} else {
	mol color colorid 7
	}
    }
    }
  }
mol selection "resid $j"
mol addrep top
}
}
