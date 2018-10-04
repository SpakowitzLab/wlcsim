module precalc_spider
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: IEEE_ARITHMETIC
    use mersenne_twister

    implicit none

    public

    !!!     hardcoded params. will need to change if certain parts of code change
    integer, parameter :: maxNumberOfSections = 40
    integer, parameter :: maxNumberOfLegs = 2*maxNumberOfSections


    type spider
        integer nSections
        integer nLegs
        integer sections(2,maxNumberOfSections)
        integer legs(3,maxNumberOfLegs)
        integer moved_sections(2,maxNumberOfSections)
    end type


contains

    subroutine load_precalc_spiders(filename,spiders,numberOfSpiders)
    implicit none
    integer, parameter :: MAXFILENAMELEN = 500
    character(MAXFILENAMELEN) filename
    type(spider), intent(out), dimension(:), allocatable :: spiders

    integer, intent(out) :: numberOfSpiders
    integer spider_id
    integer numberOfSections
    integer ii
    integer numberOfLegs


    open (UNIT = 5, FILE = filename, STATUS = 'OLD')
    read(5,*) numberOfSpiders
    allocate( spiders(numberOfSpiders+1)  )  ! +1 for the extra variable-spider
    Do spider_id = 1,numberOfSpiders
        !read number of sections
        read(5,*) numberOfSections
        spiders(spider_id)%nSections = numberOfSections
        if (numberOfSections > maxNumberOfSections) then
            print*, "Not enough sections"
            stop 1
        endif

        !read sections
        do ii = 1,numberOfSections
            read(5,*) spiders(spider_id)%sections(:,ii)
        enddo

        !read number of legs
        read(5,*) numberOfLegs
        if (numberOfLegs > maxNumberOfLegs) then
            print*, "Not enough legs"
        endif
        spiders(spider_id)%nLegs = numberOfLegs

        !read legs
        do ii = 1,numberOfLegs
            read(5,*) spiders(spider_id)%legs(:,ii)
        enddo

        !read moved sections
        do ii = 1,numberOfSections
            read(5,*) spiders(spider_id)%moved_sections(:,ii)
        enddo

    enddo
    CLOSE(5)

    end subroutine
    function get_highestNumberOfLegs(spiders,numberOfSpiders)
    implicit none
    integer get_highestNumberOfLegs
    integer, intent(in) :: numberOfSpiders
    type(spider), intent(in), dimension(numberOfSpiders) :: spiders
    integer ii
    get_highestNumberOfLegs=0
    do ii=1,numberOfSpiders
        if (spiders(ii)%nLegs>get_highestNumberOfLegs) then
            get_highestNumberOfLegs=spiders(ii)%nLegs
        endif
    enddo
    return
    end function
    subroutine print_precalc_spiders(spiders,numberOfSpiders)
    implicit none
    integer, intent(in) :: numberOfSpiders
    type(spider), intent(in) :: spiders(numberOfSpiders)

    integer spider_id
    integer numberOfSections
    integer ii
    integer numberOfLegs
    Do spider_id = 1,numberOfSpiders
        !read number of sections
        numberOfSections =  spiders(spider_id)%nSections

        print*, "spider", spider_id
        print*, "sections"
        !read sections
        do ii = 1,numberOfSections
            print*, spiders(spider_id)%sections(:,ii)
        enddo

        !read number of legs
        numberOfLegs = spiders(spider_id)%nLegs

        print*, "legs"
        !read legs
        do ii = 1,numberOfLegs
            print*, spiders(spider_id)%legs(:,ii)
        enddo

        !read moved sections
        print*, "moveable sections"
        do ii = 1,numberOfSections
            print*, spiders(spider_id)%moved_sections(:,ii)
        enddo

    enddo
    end subroutine

end module
