program main
    ! Loads in parameters from the input file.
    ! Calculates simulation parameters from input parameters.
    ! For each save point requested, uses either Bruno, Quinn, or Brad's
    ! simulators to step forwards (in time or Monte Carlo steps) and writes the
    ! output and (optionally) the simulation state needed to restart the
    ! simulation.
    !

    ! using the FLAP library (https://github.com/szaghi/FLAP) to leverage F03
    ! command-line interface features for taking input and output locations as
    ! command line arguments
    use flap
    ! structs that will hold simulation state
    use params

    implicit none

    ! CLI argument variables
    type(command_line_interface) :: cli
    character(1024)              :: infile  ! allow <=1024byte filenames
    character(1024)              :: outfile ! allow <=1024byte filenames
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
                 def='input/input'
                 error=err)
    stop_if_err(err, 'Internal argument parsing error.')
    call cli%add(switch='--output', &
                 switch_ab='-o',   &
                 help='Output filename base (absolute or relative to CWD)', &
                 required=.false., &
                 act='store',      &
                 def='data/'
                 error=err)
    stop_if_err(err, 'Internal argument parsing error.')
    call cli%get(switch='-i', val=infile, error=err)
    stop_if_err(err, 'Unable to parse input file name.')
    call cli%get(switch='-i', val=outfile, error=err)
    stop_if_err(err, 'Unable to parse output file base.')

    get_input_from_file(wlc_p=wlc_p, infile=infile)
    get_params_from_input(wlc_p)
    initialize_wlcsim_data(wlc_d=wlc_d, wlc_p=wlc_p)

    select case (wlc_p%simType)
    case ('quinn', 'parallel temper continuous parameters')
        do i=1,wlc_p%num_save_points
            wlcsim_quinn(i, wlc_d, wlc_p)
            save_simulation_state(wlc_d, wlc_p)
        enddo
    case ('brad', 'parallel temper discrete parameters', 'twist')
        do i=1,wlc_p%num_save_points
            wlcsim_brad(i, wlc_d, wlc_p)
            save_simulation_state(wlc_d, wlc_p)
        enddo
    case ('bruno', 'brownian dynamics')
        do i=1,wlc_p%num_save_points
            wlcsim_brad(i, wlc_d, wlc_p)
            save_simulation_state(wlc_d, wlc_p)
        enddo
    case default
        stop_if_err(1, 'Invalid simulation code specified.')
    end select
end program main

