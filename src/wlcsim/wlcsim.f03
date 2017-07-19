program main
    ! Loads in parameters from the input file.
    ! Calculates parameters used in simulation from input parameters.
    ! For each save point requested, uses either Bruno, Quinn, or Brad's
    ! simulators to step forwards (in time or Monte Carlo steps) and writes the
    ! output and (optionally) the simulation state needed to restart the
    ! simulation.
    !

    ! using the FLAP library (https://github.com/szaghi/FLAP) to leverage F03
    ! command-line interface features for taking input and output locations as
    ! command line arguments
    use flap, only: command_line_interface

    ! structs that will hold simulation params and state
    use params, only: wlcsim_params, wlcsim_data, &
        MAXFILENAMELEN, save_simulation_state, get_input_from_file, &
        initialize_wlcsim_data, save_parameters, setup_runtime_floats, &
        printDescription, printWindowStats

    implicit none

    integer :: i

    ! CLI argument variables
    type(command_line_interface) :: cli
    character(MAXFILENAMELEN)    :: infile
    character(MAXFILENAMELEN)    :: outfile
    integer                      :: err

    ! Simulation state
    type(wlcsim_params) :: wlc_p
    type(wlcsim_data)   :: wlc_d

    call cli%init(progname='wlcsim.exe', &
                  description='WLCSIM: A Mesoscopic Polymer Simulator')
    call cli%add(switch='--input', &
                 switch_ab='-i',   &
                 help='Input filename (absolute or relative to CWD)', &
                 required=.false., &
                 act='store',      &
                 def='input/input',&
                 error=err)
    call stop_if_err(err, 'Internal argument parsing error.')
    call cli%add(switch='--output', &
                 switch_ab='-o',   &
                 help='Output filename base (absolute or relative to CWD)', &
                 required=.false., &
                 act='store',      &
                 def='data/',      &
                 error=err)
    call stop_if_err(err, 'Internal argument parsing error.')
    call cli%get(switch='-i', val=infile, error=err)
    call stop_if_err(err, 'Unable to parse input file name.')
    call cli%get(switch='-o', val=outfile, error=err)
    call stop_if_err(err, 'Unable to parse output file base.')

    call setup_runtime_floats()

    call get_input_from_file(infile, wlc_d, wlc_p)

    call initialize_wlcsim_data(wlc_d, wlc_p)

    call save_parameters(wlc_p, outfile)

    i = 0
    call save_simulation_state(i, wlc_d, wlc_p, outfile, 'NEW')

    select case (trim(adjustL(wlc_p%codeName)))
    case ('quinn', 'parallel temper continuous parameters')
        do i=1,wlc_p%numSavePoints
            call wlcsim_quinn(i, wlc_d, wlc_p)
            call save_simulation_state(i, wlc_d, wlc_p, outfile, 'NEW')
        enddo
    case ('brad', 'parallel temper discrete parameters', 'twist')
        do i=1,wlc_p%numSavePoints
           call wlcsim_brad(wlc_d,wlc_p)
           call save_simulation_state(i, wlc_d, wlc_p, outfile, 'NEW')
           print *, 'i is', i
           print *, '*******************'
        enddo
    case ('bruno', 'brownian dynamics')
        do i=1,wlc_p%numSavePoints
            call wlcsim_bruno(i, wlc_d, wlc_p)
            call save_simulation_state(i, wlc_d, wlc_p, outfile, 'NEW')
        enddo
    case ('bruno_mc', 'simple monte carlo')
        do i=1,wlc_p%numSavePoints
            call wlcsim_bruno_mc(i, wlc_d, wlc_p)
            call save_simulation_state(i, wlc_d, wlc_p, outfile, 'NEW')
        enddo
    case ('bruno_loop', 'get looping events')
        do i=1,wlc_p%numSavePoints
            call wlcsim_bruno_looping_events(i, wlc_d, wlc_p, outfile)
            call save_simulation_state(i, wlc_d, wlc_p, outfile, 'NEW')
        enddo
    case default
        call stop_if_err(1, 'Invalid simulation code name specified.')
    end select

end program main

