program test_spider

use precalc_spider
use mersenne_twister
implicit none

integer, parameter :: MAXFILENAMELEN = 500
character(MAXFILENAMELEN) filename
type(spider), dimension(:), allocatable :: spiders
integer numberOfSpiders

filename = 'spiders'
call load_precalc_spiders(filename,spiders,numberOfSpiders)
call print_precalc_spiders(spiders,numberOfSpiders)

end program
