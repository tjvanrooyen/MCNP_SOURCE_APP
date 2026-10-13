Program mcnpsource
! Program to read source term definition from library and generate MCNP SDEF input
!-------------------------------------------------------

USE FileHandle         ! Generic File handling routines
USE Kinds

IMPLICIT NONE

!Revision Record:
!----------------------------------
!  1.)   Creation of basic structure              - 12 / 11 / 13       - Rian Prinsloo


TYPE ISO_RAD

   character(len=255)                            :: iso_name                             ! isotope name
   character(len=255)                            :: hl                                   ! halflife
   integer                                       :: num_entries                          ! number of source entries
   real(kind=doubel),allocatable,dimension(:)    :: em_yield                             ! emission yield
   real(kind=doubel),allocatable,dimension(:)    :: em_line                              ! emission line energy
   integer          ,allocatable,dimension(:)    :: em_type                              ! emission type
   character(len=255),allocatable,dimension(:)   :: em_part                              ! emission particle designation

END TYPE ISO_RAD

TYPE ISO_BET

   character(len=255)                            :: iso_name                             ! isotope name
   integer                                       :: num_entries                          ! number of source entries
   real(kind=doubel),allocatable,dimension(:)    :: em_yield                             ! emission yield
   real(kind=doubel),allocatable,dimension(:)    :: em_yield_bep                         ! emission yield
   real(kind=doubel),allocatable,dimension(:)    :: em_line                              ! emission line energy

END TYPE ISO_BET


TYPE ATOMIC_DATA

   character(len=255)                            :: Element_name                             ! isotope name
   integer                                       :: z                          ! number of source entries

END TYPE ATOMIC_DATA


TYPE ISO_MIX

   character(len=255)                            :: iso_name                             ! isotope name in the mix
   real(kind = doubel)                           :: activity                             ! activity in Bq
   integer                                       :: biso, riso                           ! beta and rad isotopic index
   real(kind = doubel)                           :: weight_factor_ben                    ! weight factor in mixture (ben)
   real(kind = doubel)                           :: weight_factor_bep                    ! weight factor in mixture (bep)

END TYPE ISO_MIX


TYPE RADFILE

   character(len=255)                            :: filename                             ! filename for radiation emmision file
   integer                                       :: numIsotopes                          ! Number of isotopes in the file (not known upfront)
   TYPE (ISO_RAD), dimension(5000)               :: iso_rad_data

END TYPE RADFILE


TYPE BETFILE

   integer                                       :: numIsotopes                          ! Number of isotopes in the file (not known upfront)
   character(len=255)                            :: filename                             ! filename for radiation emmision file
   TYPE (ISO_BET), dimension(5000)               :: iso_bet_data

END TYPE BETFILE



character(len=255)                               :: inputfilename,tempstring,mode        ! inputfile
character(len=255)                               :: outputfilename                       ! outputfile base name
character(len=255)                               :: zlistfile                            ! zlist data file
integer                                          :: num_elements, num_mix                ! number of elements in z data file, number of iso in mix
TYPE (ATOMIC_DATA),dimension(:),allocatable      :: element_data
real(kind = doubel)                              :: e_cutoff                             ! energy cut-off for ph and el.
integer                                          :: unit_main_dump,error                 ! main debug file
character(len=10)                                :: rundate,runtime                      ! date and time 
real                                             :: time,time2
TYPE (RADFILE)                                   :: radfiledata
TYPE (BETFILE)                                   :: betfiledata
TYPE (ISO_MIX), allocatable, dimension(:)        :: isomix_data
integer                                          :: i,j,k
logical                                          :: found
integer                                          :: num_meshes


   unit_main_dump = openFile('DUMP_main.txt','REPLACE','WRITE')
   CALL DATE_AND_TIME (rundate,runtime)                                            ! Get date and time of program run
   runtime = runtime(1:2)//':'//runtime(3:4)//':'//runtime(5:6)                    ! Process runtime
   call GETARG(1,inputfilename)
   write (unit_main_dump,*) 'Commencing processing '//TRIM(runtime)//' '//rundate(1:8)
   write (*,*) 'Commencing processing '//TRIM(runtime)//' '//rundate(1:8)
   call readInput(trim(adjustl(inputfilename)),error)
   
   if (error == 0) then


      write(unit_main_dump,*) 'Input file processed.'
      
      CALL CPU_TIME(time)
      call readRadFile(trim(adjustl(radfiledata%filename)),error)
      call readBetFile(trim(adjustl(betfiledata%filename)),error)
      call readElements(trim(adjustl(zlistfile)),error)
      if (trim(adjustl(mode)) == 'MIX') call buildMixIso()
      call sortRadData()
      if (trim(adjustl(mode)) == 'ALL') then
         write(*,*) 'Dumping data for all isotopes in RAD file...'
         do i=1, radfiledata%numIsotopes
            call writeMCNPSource(trim(adjustl(outputfilename)),i)
         end do
      else if (trim(adjustl(mode)) =='MIX') then
         write(*,*) 'Constructing isotopic mix '
         call writeMCNPSource(trim(adjustl(outputfilename)),radfiledata%numIsotopes + 1)

      else

         write(*,*) 'Searching RAD file for isotope: ',trim(adjustl(mode))
         found = .false.
         do i=1, radfiledata%numIsotopes
            if (trim(adjustl(mode)) == trim(adjustl(radfiledata%iso_rad_data(i)%iso_name))) then
               write(*,*) 'Writing output only for isotope: ',trim(adjustl(radfiledata%iso_rad_data(i)%iso_name))
               found = .true.
               call writeMCNPSource(trim(adjustl(outputfilename)),i)
            end if
         end do
         if (.not.(found)) then
            write(*,*) trim(adjustl(mode)),' not found in RAD file.'
            stop
         end if

      end if
      call CPU_TIME(time2)
      
   else
      write(unit_main_dump,*) 'ERROR reading input file.'
      write(*,*) 'ERROR reading input file.'
   end if

   CALL DATE_AND_TIME (rundate,runtime)                                            ! Get date and time of program run
   runtime = runtime(1:2)//':'//runtime(3:4)//':'//runtime(5:6)                    ! Process runtime
   write (unit_main_dump,*) 'Processing finished at '//TRIM(runtime)//' '//rundate(1:8)
   write (*,*) 'Calculation finished at '//TRIM(runtime)//' '//rundate(1:8)
   write(*,*) 'Calculational time: ',time2-time
   write(unit_main_dump,*) 'Calculational time: ',time2-time

CONTAINS


!----------------------------------------------------------------------
  Subroutine readInput(   &
       inputfilename,error)
    
!Subroutine to read the Generic input file
!----------------------------------------------------------------------
    character(len=*),INTENT(IN)               :: inputfilename
    integer, intent(OUT)                      :: error
    integer                                   :: unit_input, i
    
    unit_input = openFile(trim(adjustl(inputfilename)),"OLD","READ")
    error = 0
    
    read(unit_input,*,IOSTAT=error) radfiledata%filename
    read(unit_input,*,IOSTAT=error) betfiledata%filename
    read(unit_input,*,IOSTAT=error) zlistfile
    read(unit_input,*,IOSTAT=error) e_cutoff
    read(unit_input,*,IOSTAT=error) mode
    read(unit_input,*,IOSTAT=error) outputfilename
    
    write(*,*) 'RAD file to process: ',trim(adjustl(radfiledata%filename))
    write(*,*) 'BET file to process: ',trim(adjustl(betfiledata%filename))
    write(*,*) 'Atomic number per isotope from ',trim(adjustl(zlistfile))
    write(*,*) 'Output to be written to: ',trim(adjustl(outputfilename)),'_isoname.txt'

    if (trim(adjustl(mode)) == 'MIX') then
       
      read(unit_input,*,IOSTAT=error) num_mix, num_meshes 
      write(*,*) 'Reading ',num_mix,' isotopes for single source definition.'
      allocate(isomix_data(num_mix))
      do i=1, num_mix
         read(unit_input,*,IOSTAT=error) isomix_data(i)%iso_name,isomix_data(i)%activity
         write(*,*) i,trim(adjustl(isomix_data(i)%iso_name)),isomix_data(i)%activity
      end do

    end if

    call closeFile(unit_input)
    
  End Subroutine readInput
!----------------------------------------------------------------------
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  Subroutine readElements(   &
       inputfilename,error)
    
!Subroutine to read the Generic input file
!----------------------------------------------------------------------
    character(len=*),INTENT(IN)               :: inputfilename
    integer, intent(OUT)                      :: error
    integer                                   :: unit_input,i
    
    unit_input = openFile(trim(adjustl(inputfilename)),"OLD","READ")
    error = 0
    
    read(unit_input,*,IOSTAT=error) num_elements
    allocate(element_data(num_elements))
    do i=1, num_elements
        read(unit_input,*, IOSTAT=error) element_data(i)%element_name,element_data(i)%z
        if (error > 0) then
           write(*,*) 'Error reading element data file on entry: ',i
        end if
       
    end do


    call closeFile(unit_input)
    write(*,*) 'Read element data file.'

  End Subroutine readElements
!----------------------------------------------------------------------

!----------------------------------------------------------------------
  Subroutine buildMixIso(   &
       )
    
!Subroutine to construct mixed effective isotope from provided list (in num_iso + 1 pos)
!----------------------------------------------------------------------
    integer                                   :: unit_input,i,j,k, total, startlist,endlist, counter
    logical                                   :: found
    real(kind = doubel)                       :: norm_factor, energy_max, energy_min, energy_delta, yield, yieldt, tempyield_ben, tempyield_bep
    real(kind = doubel)                       :: yieldt_ben, yieldt_bep, norm_factor_ben, norm_factor_bep
    
    
    radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%num_entries = 0

! Normalisation factor
    norm_factor = 0.0_doubel
    do i=1, num_mix
       norm_factor = norm_factor + isomix_data(i)%activity
    end do


! determine total size
    do i=1, num_mix
       found = .false.

       
       do j=1, radfiledata%numIsotopes
!          write(900,*) 'Comparing: ',trim(adjustl(radfiledata%iso_rad_data(j)%iso_name)),' to ',trim(adjustl(isomix_data(i)%iso_name))
          if (trim(adjustl(radfiledata%iso_rad_data(j)%iso_name)) == trim(adjustl(isomix_data(i)%iso_name))) then
             write(*,*) 'Found match for mix isotope: ',trim(adjustl(radfiledata%iso_rad_data(j)%iso_name))
             isomix_data(i)%riso = j
! increment number of entries
             radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%num_entries =  &   ! increment effective isotope num_entries
                  radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%num_entries + radfiledata%iso_rad_data(j)%num_entries
             found = .true.
             exit

          end if
          
          
       end do
       if (.not.found) then
          write(*,*) 'Error: Isotope in mix ',trim(adjustl(isomix_data(i)%iso_name)),' not found in rad file.'
          stop
       end if
       
       
    end do
    
! rad file data - concatenate selected isotopes
    allocate(radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%em_yield(radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%num_entries))
    allocate(radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%em_line(radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%num_entries))
    allocate(radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%em_type(radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%num_entries))
    allocate(radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%em_part(radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%num_entries))    
    radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%iso_name = 'MIX'
    counter = 0
    do i=1, num_mix

       
       do j=1, radfiledata%numIsotopes
!          write(900,*) 'Comparing: ',trim(adjustl(radfiledata%iso_rad_data(j)%iso_name)),' to ',trim(adjustl(isomix_data(i)%iso_name))
          if (trim(adjustl(radfiledata%iso_rad_data(j)%iso_name)) == trim(adjustl(isomix_data(i)%iso_name))) then
! append elements  
           
             startlist = counter + 1
             endlist = counter + radfiledata%iso_rad_data(j)%num_entries
             
             radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%em_yield(startlist:endlist) = &
                  radfiledata%iso_rad_data(j)%em_yield(1:radfiledata%iso_rad_data(j)%num_entries) * isomix_data(i)%activity / norm_factor
             radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%em_line(startlist:endlist) =  &
                  radfiledata%iso_rad_data(j)%em_line(1:radfiledata%iso_rad_data(j)%num_entries)
             radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%em_type(startlist:endlist) =  &
                  radfiledata%iso_rad_data(j)%em_type(1:radfiledata%iso_rad_data(j)%num_entries)
             radfiledata%iso_rad_data(radfiledata%numIsotopes + 1)%em_part(startlist:endlist) =  &
                  radfiledata%iso_rad_data(j)%em_part(1:radfiledata%iso_rad_data(j)%num_entries)
             
! increment number of entries
             counter =  &   ! increment effective isotope num_entries
                  counter + radfiledata%iso_rad_data(j)%num_entries
             exit
          end if
          
          
       end do
       
    end do
    
    
    write(*,*) 'Built mixed rad equivalent isotope'

! build beta equivalent isotope

! find equivalent beta file iso index 
    energy_max = 0.0_doubel
    energy_min = 1E20_doubel
    do i=1, num_mix
       found = .false.
       
       isomix_data(i)%biso = 0
       
       do j=1, betfiledata%numIsotopes
          if (trim(adjustl(betfiledata%iso_bet_data(j)%iso_name)) == trim(adjustl(isomix_data(i)%iso_name))) then
             isomix_data(i)%biso = j
             found = .true.
             if (maxval(betfiledata%iso_bet_data(j)%em_line(:)) > energy_max) energy_max = maxval(betfiledata%iso_bet_data(j)%em_line(:))
             if (minval(betfiledata%iso_bet_data(j)%em_line(:)) < energy_min) energy_min = minval(betfiledata%iso_bet_data(j)%em_line(:))
             exit
          end if
       end do
       if (.not.(found)) write(*,*) 'Warning: no beta data for requested isotope: ',trim(adjustl(isomix_data(i)%iso_name))
    end do
    write(*,'(A43,F9.4,A4,F9.4,A4)') 'Global energy limits for beta data set as: ',energy_min,' to ',energy_max,' MeV'
    betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%num_entries = num_meshes
    allocate(betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_yield(betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%num_entries))
    allocate(betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_yield_bep(betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%num_entries))
    allocate(betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_line(betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%num_entries))

    write(*,*) 'Populating mesh with ',num_meshes,' equidistant points.'
    energy_delta = (energy_max - energy_min) / (num_meshes - 1)
    
    
    
    betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_line(1)=energy_min
    do i=2, num_meshes-1
       betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_line(i) =  betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_line(i-1) + energy_delta
    end do
    betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_line(num_meshes)=energy_max
    
!interpolate on new energy grid
    norm_factor = 0.0_doubel
    norm_factor_ben = 0.0_doubel
    norm_factor_bep = 0.0_doubel
    
    do j=1, num_mix
       write(901,*) trim(adjustl(isomix_data(j)%iso_name)), isomix_data(j)%riso
       tempyield_ben = 0.0_doubel
       tempyield_bep = 0.0_doubel
       call totalYield(isomix_data(j)%riso,"B-",tempyield_ben)
       call totalYield(isomix_data(j)%riso,"B+",tempyield_bep)
       write(901,*) j,tempyield_ben,tempyield_bep
       norm_factor = norm_factor + isomix_data(j)%activity
       norm_factor_ben = norm_factor_ben + isomix_data(j)%activity * tempyield_ben
       norm_factor_bep = norm_factor_bep + isomix_data(j)%activity * tempyield_bep
    end do
!    write(901,*) 'Totals bes: ',norm_factor_ben,norm_factor_bep
    do j=1, num_mix
       tempyield_ben = 0.0_doubel
       tempyield_bep = 0.0_doubel
       call totalYield(isomix_data(j)%riso,"B-",tempyield_ben)
       call totalYield(isomix_data(j)%riso,"B+",tempyield_bep)
       isomix_data(j)%weight_factor_ben = 1.0_doubel
       isomix_data(j)%weight_factor_bep = 1.0_doubel
       if (norm_factor /= 0.0_doubel) then
          if (norm_factor_ben /= 0.0_doubel) &
               isomix_data(j)%weight_factor_ben = (isomix_data(j)%activity * tempyield_ben) / (norm_factor_ben + norm_factor_bep)
          if (norm_factor_bep /= 0.0_doubel) &
               isomix_data(j)%weight_factor_bep = (isomix_data(j)%activity * tempyield_bep) / (norm_factor_ben + norm_factor_bep)
       end if
!       write(901,*) 'Isotope: ',trim(adjustl(isomix_data(j)%iso_name))
!       write(901,*) 'BEN: ',isomix_data(j)%weight_factor_ben
!       write(901,*) 'BEP: ',isomix_data(j)%weight_factor_bep
    end do
    
    
    
    do i=1,num_meshes
       yieldt_ben = 0.0_doubel
       yieldt_bep = 0.0_doubel

       do j=1,num_mix
          if (isomix_data(j)%biso > 0) then
!             write(901,*) 'Isotope: ',isomix_data(j)%biso,isomix_data(j)%riso, 
             call interpolate(i,isomix_data(j)%biso,yield)
             if (yield > 0.0_doubel) then
                yieldt_ben = yieldt_ben + yield * isomix_data(j)%weight_factor_ben
                yieldt_bep = yieldt_bep + yield * isomix_data(j)%weight_factor_bep
  
!                write(901,*) i,betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_line(i),j,yieldt_ben,yieldt_bep,j
                
             end if
          end if

       end do

       betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_yield(i) = yieldt_ben
       betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_yield_bep(i) = yieldt_bep
    end do
       
  End Subroutine buildMixIso

!----------------------------------------------------------------------
  Subroutine interpolate(mesh,isotope,value   &
       )
    
!Subroutine to interpolate on energy grid (linear)
!----------------------------------------------------------------------
    integer, INTENT(IN)                                 :: mesh
    integer, INTENT(IN)                                 :: isotope
    integer                                             :: i,j,left,right
    real(kind = doubel), INTENT(OUT)                    :: value
    real(kind = doubel)                                 :: m,c
    
   
! find left and rightend point

    if ((betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_line(mesh)  < betfiledata%iso_bet_data(isotope)%em_line(1))) then
       left = 0
       right = 0
    end if
    if ((betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_line(mesh)  > betfiledata%iso_bet_data(isotope)%em_line(betfiledata%iso_bet_data(isotope)%num_entries))) then
       left = 0
       right = 0
    end if
    if ((betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_line(mesh)  == betfiledata%iso_bet_data(isotope)%em_line(betfiledata%iso_bet_data(isotope)%num_entries))) then
       left = betfiledata%iso_bet_data(isotope)%num_entries
       right = betfiledata%iso_bet_data(isotope)%num_entries
    end if
    do k=1,betfiledata%iso_bet_data(isotope)%num_entries-1
       
       if (betfiledata%iso_bet_data(isotope)%em_line(k) == betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_line(mesh)) then
          left = k
          right = k
       end if
       
       if ((betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_line(mesh) > betfiledata%iso_bet_data(isotope)%em_line(k)).and.&
            (betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_line(mesh) < betfiledata%iso_bet_data(isotope)%em_line(k+1))) then
          left = k
          right = k+1
       end if
       
    end do
    write(999,*) 'For isotope: ',trim(adjustl(betfiledata%iso_bet_data(isotope)%iso_name)),' at ',betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_line(mesh)
    if ((left == 0).and.(right == 0)) then
       value = 0.0_doubel
       write(999,*) 'Out of range: ',value
       
    else if (left == right) then
       value = betfiledata%iso_bet_data(isotope)%em_yield(left)
    else
       m = (betfiledata%iso_bet_data(isotope)%em_yield(right) - betfiledata%iso_bet_data(isotope)%em_yield(left)) / &
            (betfiledata%iso_bet_data(isotope)%em_line(right) - betfiledata%iso_bet_data(isotope)%em_line(left))
       c = betfiledata%iso_bet_data(isotope)%em_yield(right) - m * betfiledata%iso_bet_data(isotope)%em_line(right)
       value = m * betfiledata%iso_bet_data(betfiledata%numIsotopes + 1)%em_line(mesh) + c
!       write(999,*) 'Interpolate',value,betfiledata%iso_bet_data(isotope)%em_yield(left),betfiledata%iso_bet_data(isotope)%em_yield(right)
    end if
 
 
    
    
  End Subroutine interpolate
  
!----------------------------------------------------------------------



  Subroutine readRadFile(   &
       inputfilename,error)
    
!Subroutine to read the RAD data file
!----------------------------------------------------------------------
! Revision record:
!--------------------------------------------------------------------------------------
! No.   Date          Programmer     Description of changes
!--------------------------------------------------------------------------------------
! 0     2013/11/11    R H Prinsloo    Original
!
    character(len=*),INTENT(IN)               :: inputfilename
    integer, intent(OUT)                      :: error
    integer                                   :: unit_rad
    logical                                   :: eof_flag
    integer                                   :: count,i    
    character(len=255)                        :: dummy
    
    unit_rad = openFile(trim(adjustl(inputfilename)),"OLD","READ")
    eof_flag = .false.
    count = 0
    do while (.not.(eof_flag))
       count = count + 1
       read(unit_rad,*) dummy
       if (trim(adjustl(dummy)) == 'END') then
          eof_flag = .true.
       else
          backspace(unit = unit_rad)
          read(unit_rad,*, IOSTAT=error) radfiledata%iso_rad_data(count)%iso_name, radfiledata%iso_rad_data(count)%hl, radfiledata%iso_rad_data(count)%num_entries
          if (error > 0) then
             write(*,*) 'Error reading isotope: ',count
             stop
          else
             allocate(radfiledata%iso_rad_data(count)%em_yield(radfiledata%iso_rad_data(count)%num_entries))
             allocate(radfiledata%iso_rad_data(count)%em_line(radfiledata%iso_rad_data(count)%num_entries))
             allocate(radfiledata%iso_rad_data(count)%em_type(radfiledata%iso_rad_data(count)%num_entries))
             allocate(radfiledata%iso_rad_data(count)%em_part(radfiledata%iso_rad_data(count)%num_entries))
             
             
             do i=1, radfiledata%iso_rad_data(count)%num_entries
                read(unit_rad,*, IOSTAT=error) radfiledata%iso_rad_data(count)%em_type(i), radfiledata%iso_rad_data(count)%em_yield(i), radfiledata%iso_rad_data(count)%em_line(i), radfiledata%iso_rad_data(count)%em_part(i)
                if (error > 0) then
                   write(*,*) 'Error reading RAD isotope: ',count,' at entry: ',i,' (',trim(adjustl(radfiledata%iso_rad_data(count)%iso_name)),')'
                   stop
                end if
             end do
          end if
          
          
       end if
       
       
    end do
    radfiledata%numIsotopes = count
    write(*,*) 'Successfully processed RAD file. Read isotopes: ',radfiledata%numIsotopes
    
    
  End Subroutine readRadFile
!----------------------------------------------------------------------
!----------------------------------------------------------------------


!----------------------------------------------------------------------
  Subroutine readBetFile(   &
       inputfilename,error)
    
!Subroutine to read the RAD data file
!----------------------------------------------------------------------
! Revision record:
!--------------------------------------------------------------------------------------
! No.   Date          Programmer     Description of changes
!--------------------------------------------------------------------------------------
! 0     2013/11/11    R H Prinsloo    Original
!
    character(len=*),INTENT(IN)               :: inputfilename
    integer, intent(OUT)                      :: error
    integer                                   :: unit_bet
    logical                                   :: eof_flag
    integer                                   :: count,i    
    character(len=255)                        :: dummy
    real(kind = doubel)                       :: prev
    
    unit_bet = openFile(trim(adjustl(inputfilename)),"OLD","READ")
    eof_flag = .false.
    count = 0
    do while (.not.(eof_flag))
       count = count + 1
       read(unit_bet,*) dummy
       if (trim(adjustl(dummy)) == 'END') then
          eof_flag = .true.
       else
          backspace(unit = unit_bet)
          read(unit_bet,*, IOSTAT=error) betfiledata%iso_bet_data(count)%iso_name, betfiledata%iso_bet_data(count)%num_entries
          
          if (error > 0) then
             write(*,*) 'Error reading isotope: ',count
             stop
          else
             allocate(betfiledata%iso_bet_data(count)%em_line(betfiledata%iso_bet_data(count)%num_entries))
             allocate(betfiledata%iso_bet_data(count)%em_yield(betfiledata%iso_bet_data(count)%num_entries))
            
             prev = -1.0_doubel
             do i=1, betfiledata%iso_bet_data(count)%num_entries
                read(unit_bet,*, IOSTAT=error) betfiledata%iso_bet_data(count)%em_line(i), betfiledata%iso_bet_data(count)%em_yield(i)
                if (error > 0) then
                   write(*,*) 'Error reading BET isotope: ',count,' at entry: ',i,' (',trim(adjustl(betfiledata%iso_bet_data(count)%iso_name)),')'
                   stop
                else
                   
                   if (betfiledata%iso_bet_data(count)%em_line(i) <= prev) then
                      write(*,*) 'Error: Energy not monotonically increasing for BET data. Isotope ',count,' line ',i
                      stop
                   else
                      prev = betfiledata%iso_bet_data(count)%em_line(i)
                   end if


                end if
             end do
          end if
          
          
       end if
       
       
    end do
    betfiledata%numIsotopes = count
    write(*,*) 'Successfully processed BET file. Read isotopes: ',betfiledata%numIsotopes
   
    
  End Subroutine readBetFile
!----------------------------------------------------------------------
!----------------------------------------------------------------------


!=====================================================================================!
  subroutine  sortRadData ()
!=====================================================================================!
! Routine sorts  in ascending order.
!-------------------------------------------------------------------------------------!
!
! Revision record:
!--------------------------------------------------------------------------------------
! No.   Date          Programmer     Description of changes
!--------------------------------------------------------------------------------------
! 0     2013/11/11    R H Prinsloo    Original
!
    implicit none
! Dummy variables
! Local variables
    integer                     :: i,j,k, n_iso                   ! Array indices
!-------------------------------------------------------------------------------------!
    
    write(*,*) 'Sorting RAD data according to energy...'   
    n_iso = radfiledata%numIsotopes
    if (trim(adjustl(mode)) == 'MIX') n_iso = n_iso + 1

    do k=1, n_iso
       
       do i = 1, radfiledata%iso_rad_data(k)%num_entries
          do j = i+1, radfiledata%iso_rad_data(k)%num_entries
             if (radfiledata%iso_rad_data(k)%em_line(j) < radfiledata%iso_rad_data(k)%em_line(i)) then        ! If smaller array element is found
                call RSWAP (k,j,i)     ! Swap array elements
             end if
          end do
       end do
       

!       do i=1, radfiledata%iso_rad_data(k)%num_entries
!          write(unit_main_dump,*)   radfiledata%iso_rad_data(k)%em_line(i)
!       end do
    end do

    write(*,*) 'Successfully sorted RAD data according to energy'
     end subroutine SortRaddata

!=====================================================================================!

!=====================================================================================!
  subroutine RSWAP(k, &
       e1,          &
       e2)
!=====================================================================================!
! Routine swaps double precision real array elements.
!-------------------------------------------------------------------------------------!
!
! Revision record:
!--------------------------------------------------------------------------------------
! No.   Date          Programmer     Description of changes
!--------------------------------------------------------------------------------------
! 0     2013/11/11    R H Prinsloo    Original
!
    implicit none
! Dummy variables
    integer, intent(IN)    :: e1,e2,k  ! Array elements to be swapped
! Local variables
    real(kind = doubel)    :: d1,d2
    integer                :: d3
    character(len= 255)    :: d4
!-------------------------------------------------------------------------------------!
    
!    dummy = e1
    d1 = radfiledata%iso_rad_data(k)%em_line(e1)
    d2 = radfiledata%iso_rad_data(k)%em_yield(e1)
    d3 = radfiledata%iso_rad_data(k)%em_type(e1)
    d4 = radfiledata%iso_rad_data(k)%em_part(e1)

!    e1    = e2
    radfiledata%iso_rad_data(k)%em_line(e1) = radfiledata%iso_rad_data(k)%em_line(e2)
    radfiledata%iso_rad_data(k)%em_yield(e1) = radfiledata%iso_rad_data(k)%em_yield(e2) 
    radfiledata%iso_rad_data(k)%em_type(e1) = radfiledata%iso_rad_data(k)%em_type(e2)
    radfiledata%iso_rad_data(k)%em_part(e1) = radfiledata%iso_rad_data(k)%em_part(e2)
    
!    e2    = dummy
    radfiledata%iso_rad_data(k)%em_line(e2) = d1
    radfiledata%iso_rad_data(k)%em_yield(e2) = d2
    radfiledata%iso_rad_data(k)%em_type(e2) = d3
    radfiledata%iso_rad_data(k)%em_part(e2) = d4


  end subroutine RSWAP
!=====================================================================================!
!=====================================================================================!

!=====================================================================================!
  subroutine  writeMCNPSource (outputfile,iso)
!=====================================================================================!
! Routine dumps MCNP sourcedef card in output file
!-------------------------------------------------------------------------------------!
!
! Revision record:
!--------------------------------------------------------------------------------------
! No.   Date          Programmer     Description of changes
!--------------------------------------------------------------------------------------
! 0     2013/11/11    R H Prinsloo    Original
!
    implicit none
! Dummy variables
! Local variables
    integer                         :: i,j,k,append, counter   ! Array indices
    integer                         :: unit_output             ! Array indices
    character(len=*), INTENT(IN)    :: outputfile
    integer, INTENT(IN)             :: iso
    integer                         :: biso
    real(kind=doubel)               :: yield,tempyield
    character(len=255)              :: zaid,description
    character(len=255),dimension(9) :: legals
    logical                         :: found,ff,mix
!-------------------------------------------------------------------------------------!
    

    write(*,*) 'Dumping data to: ',trim(adjustl(outputfile))//'_'//trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name))
    unit_output = openFile(trim(adjustl(outputfile))//'_'//trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name))//'.txt','REPLACE','WRITE')


! writing generic header for file

    write(unit_output,'(A)') 'c'
    write(unit_output,'(A22,A10)') 'c == GENERAL SDEF FOR ',trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name))
    write(unit_output,'(A)') 'c'
    write(unit_output,'(A)') 'mode   N  P  E  A  #'
    write(unit_output,'(A)') 'c'
    write(unit_output,'(A)') 'c Cutoff Energies for Charged Particle Transport:'
    write(unit_output,'(A16,ES15.7,A)') 'c cut:E      J  ',e_cutoff,'  J  J     $ Cutoff energy for electron       transport '
    write(unit_output,'(A16,ES15.7,A)') 'c cut:A      J  ',e_cutoff,'  J  J     $ Cutoff energy for alpha-particle transport '
    write(unit_output,'(A16,ES15.7,A)') 'c cut:#      J  ',e_cutoff,'  J  J     $ Cutoff energy for heavy ion      transport '
    write(unit_output,'(A)') 'c'
    write(unit_output,'(A)') 'SDEF   POS = 0  0  0'
    write(unit_output,'(A)') '       PAR = D7'
    write(unit_output,'(A)') '       ERG = FPAR  = D10'
    write(unit_output,'(A)') 'c'
    write(unit_output,'(A)') '#       SI7      SP7               DS10'
    write(unit_output,'(A)') '         L        D                 S'
    
! find equivalent beta file iso index    
    biso = 0
    mix = .false.
    do i=1, betfiledata%numIsotopes
       if (trim(adjustl(betfiledata%iso_bet_data(i)%iso_name)) == trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name))) then
          biso = i
!          write(*,*) 'Found equivalent isotope: ',trim(adjustl(betfiledata%iso_bet_data(i)%iso_name)),biso
          exit
       end if
    end do
    if (iso == radfiledata%numIsotopes+1) then
       biso = betfiledata%numIsotopes + 1      ! mixed source is written to iso + 1 position in both files.
       mix = .true.
    end if

! writing particle distribution card
! Loop counter for particle distribution card, then energy distribution per particle

    do j=1, 2 

! neutron source
       
       append = 1
       ff = .false.
       yield = 0.0_doubel
       call totalYield(iso,'N',yield)
       if (yield > 0.0_doubel) ff = .true.
       
       description='Neutrons'
       yield = 0.0_doubel
       call totalYield(iso,'N',yield)
       if (yield > 0.0_doubel)  then
          if (j == 1) write(unit_output,'(A17,ES14.7,A5,I2,A8,A)') '         N       ',yield,'     ',20+append,'      $ ',trim(adjustl(description))
          if (j == 2) then
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A2,A10,A2,A20)') 'c ',trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name)),' ',trim(adjustl(description))
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A10,I2,A17,I2)') '#       SI',20+append,'               SP',20+append
             write(unit_output,*) '       L                  D'
!             do i = radfiledata%iso_rad_data(iso)%num_entries,1,-1
!                if ((trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(i))) == 'N').and.(radfiledata%iso_rad_data(iso)%em_line(i) > e_cutoff)) write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',radfiledata%iso_rad_data(iso)%em_line(i),'       ', radfiledata%iso_rad_data(iso)%em_yield(i)
!             end do
             i = radfiledata%iso_rad_data(iso)%num_entries
             do while (i > 0) 
                tempyield = 0.0_doubel
                counter = 1
                do k = i,1,-1
                   if ((radfiledata%iso_rad_data(iso)%em_line(k) > e_cutoff).and.(radfiledata%iso_rad_data(iso)%em_line(k)==radfiledata%iso_rad_data(iso)%em_line(i))) then
                      if (trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(i))) == 'N') tempyield = tempyield + radfiledata%iso_rad_data(iso)%em_yield(k)
                      counter = counter + 1
                   else
                      exit
                   end if
                end do
                if (tempyield > 0.0_doubel) write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',radfiledata%iso_rad_data(iso)%em_line(i),'       ', tempyield

                if (counter == 1) then
                   i=i - (counter)
                else
                   i=i - (counter - 1)
                end if


             end do

             if (ff) then
                  write(unit_output,'(A)') 'c'
                  write(unit_output,'(A)') 'Best available N-spectrum: Use Watt Fission neutron spectrum with parameter values a and b taken from MCNP manual'
                  write(unit_output,'(A)') 'c'
             end if
             
          end if
       end if
       
!  photon source
       append = 2
       yield = 0.0_doubel
       description='Photons'
       legals(1:4)=(/'X','G','DG','PG'/)
       do i=1,4
          call totalYield(iso,trim(adjustl(legals(i))),yield)
       end do
!   Ignore annihilation quantum photom
!   call totalYield(iso,'AQ',yield)     
       if (yield > 0.0_doubel)  then
          if (j == 1) write(unit_output,'(A17,ES14.7,A5,I2,A8,A)') '         P       ',yield,'     ',20+append,'      $ ',trim(adjustl(description))
          if (j == 2) then
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A2,A10,A2,A20)') 'c ',trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name)),' ',trim(adjustl(description))
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A10,I2,A17,I2)') '#       SI',20+append,'               SP',20+append
             write(unit_output,*) '       L                  D'
!             do i = radfiledata%iso_rad_data(iso)%num_entries,1,-1
!                if ((any(legals(1:4) == trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(i))))).and.(radfiledata%iso_rad_data(iso)%em_line(i) > e_cutoff)) write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',radfiledata%iso_rad_data(iso)%em_line(i),'       ', radfiledata%iso_rad_data(iso)%em_yield(i)
!             end do

             i = radfiledata%iso_rad_data(iso)%num_entries
             do while (i > 0) 
                tempyield = 0.0_doubel
                counter = 1
                do k = i,1,-1
                   if ((radfiledata%iso_rad_data(iso)%em_line(k) > e_cutoff).and.(radfiledata%iso_rad_data(iso)%em_line(k)==radfiledata%iso_rad_data(iso)%em_line(i))) then
                      if (any(legals(1:4) == trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(k))))) tempyield = tempyield + radfiledata%iso_rad_data(iso)%em_yield(k)
                      counter = counter + 1
                   else
                      exit
                   end if
                end do
                if (tempyield > 0.0_doubel) write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',radfiledata%iso_rad_data(iso)%em_line(i),'       ', tempyield
                if (counter == 1) then
                   i=i - (counter)
                else
                   i=i - (counter - 1)
                end if


             end do
             
          end if
       end if
       
! line electron source
       append = 3
       yield = 0.0_doubel
       description='Line-Electrons'
       legals(1:2)=(/'IE','AE'/)
       do i=1,2
          call totalYield(iso,trim(adjustl(legals(i))),yield)
       end do
       if (yield > 0.0_doubel)  then
          if (j == 1) write(unit_output,'(A17,ES14.7,A5,I2,A8,A)') '         E       ',yield,'     ',20+append,'      $ ',trim(adjustl(description))
          if (j == 2) then
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A2,A10,A2,A20)') 'c ',trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name)),' ',trim(adjustl(description))
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A10,I2,A17,I2)') '#       SI',20+append,'               SP',20+append
             write(unit_output,*) '       L                  D'
!             do i = radfiledata%iso_rad_data(iso)%num_entries,1,-1
!                if ((any(legals(1:2) == trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(i))))).and.(radfiledata%iso_rad_data(iso)%em_line(i) > e_cutoff)) write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',radfiledata%iso_rad_data(iso)%em_line(i),'       ', radfiledata%iso_rad_data(iso)%em_yield(i)
!             end do
             i = radfiledata%iso_rad_data(iso)%num_entries
             do while (i > 0) 
                tempyield = 0.0_doubel
                counter = 1
                do k = i,1,-1
                   if ((radfiledata%iso_rad_data(iso)%em_line(k) > e_cutoff).and.(radfiledata%iso_rad_data(iso)%em_line(k)==radfiledata%iso_rad_data(iso)%em_line(i))) then
                      if (any(legals(1:2) == trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(k))))) tempyield = tempyield + radfiledata%iso_rad_data(iso)%em_yield(k)
                      counter = counter + 1
                   else
                      exit
                   end if
                end do
                if (tempyield > 0.0_doubel) write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',radfiledata%iso_rad_data(iso)%em_line(i),'       ', tempyield
                if (counter == 1) then
                   i=i - (counter)
                else
                   i=i - (counter - 1)
                end if


             end do
             
          end if
       end if
       
! B- electron source
       append = 4
       yield = 0.0_doubel
       description='Beta-Electrons (-)'
       call totalYield(iso,'B-',yield)
       if (yield > 0.0_doubel)  then
          if (j == 1) write(unit_output,'(A17,ES14.7,A5,I2,A8,A)') '         E       ',yield,'     ',20+append,'      $ ',trim(adjustl(description))
          if (j == 2) then
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A2,A10,A2,A20)') 'c ',trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name)),' ',trim(adjustl(description))
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A10,I2,A17,I2)') '#       SI',20+append,'               SP',20+append
             write(unit_output,*) '       A                  D'
             if (biso /= 0) then
                do i = 1,betfiledata%iso_bet_data(biso)%num_entries
                   if (betfiledata%iso_bet_data(biso)%em_yield(i) > 0.0_doubel) &
                   write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',betfiledata%iso_bet_data(biso)%em_line(i),'       ', betfiledata%iso_bet_data(biso)%em_yield(i)
                end do
             else
                write(*,*) 'Beta data on isotope: ',trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name)),' not found.'
                stop
             end if
          end if
       end if
       
! B+ electron source
       append = 5
       yield = 0.0_doubel
       description='Beta-Positrons (+)'
       call totalYield(iso,'B+',yield)
       if (yield > 0.0_doubel)  then
          if (j == 1) write(unit_output,'(A17,ES14.7,A5,I2,A8,A)') '        -E       ',yield,'     ',20+append,'      $ ',trim(adjustl(description))
          if (j == 2) then
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A2,A10,A2,A20)') 'c ',trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name)),' ',trim(adjustl(description))
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A10,I2,A17,I2)') '#       SI',20+append,'               SP',20+append
             write(unit_output,*) '       A                  D'
             if (biso /= 0) then
                do i = 1,betfiledata%iso_bet_data(biso)%num_entries
                   if (mix) then
                   if (betfiledata%iso_bet_data(biso)%em_yield_bep(i) > 0.0_doubel) &
                      write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',betfiledata%iso_bet_data(biso)%em_line(i),'       ', betfiledata%iso_bet_data(biso)%em_yield_bep(i)
                   else
                      write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',betfiledata%iso_bet_data(biso)%em_line(i),'       ', betfiledata%iso_bet_data(biso)%em_yield(i)
                   end if
                end do
             else
                write(*,*) 'Beta data on isotope: ',trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name)),' not found.'
                stop
             end if
          end if
       end if
       
! Alpha source
       append = 6
       yield = 0.0_doubel
       description='Alpha-Particles'
       call totalYield(iso,'A',yield)
       if (yield > 0.0_doubel)  then
          if (j == 1) write(unit_output,'(A17,ES14.7,A5,I2,A8,A)') '         A       ',yield,'     ',20+append,'     $ ',trim(adjustl(description))
          if (j == 2) then
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A2,A10,A2,A20)') 'c ',trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name)),' ',trim(adjustl(description))
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A10,I2,A17,I2)') '#       SI',20+append,'               SP',20+append
             write(unit_output,*) '       L                  D'
!             do i = radfiledata%iso_rad_data(iso)%num_entries,1,-1
!                if (('A' == trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(i)))).and.(radfiledata%iso_rad_data(iso)%em_line(i) > e_cutoff)) write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',radfiledata%iso_rad_data(iso)%em_line(i),'       ', radfiledata%iso_rad_data(iso)%em_yield(i)
!             end do

             i = radfiledata%iso_rad_data(iso)%num_entries
             do while (i > 0) 
                tempyield = 0.0_doubel
                counter = 1
                do k = i,1,-1
                   if ((radfiledata%iso_rad_data(iso)%em_line(k) > e_cutoff).and.(radfiledata%iso_rad_data(iso)%em_line(k)==radfiledata%iso_rad_data(iso)%em_line(i))) then
                      if ('A' == trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(i)))) tempyield = tempyield + radfiledata%iso_rad_data(iso)%em_yield(k)
                      counter = counter + 1
                   else
                      exit
                   end if
                end do
                if (tempyield > 0.0_doubel) write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',radfiledata%iso_rad_data(iso)%em_line(i),'       ', tempyield
                if (counter == 1) then
                   i=i - (counter)
                else
                   i=i - (counter - 1)
                end if


             end do

             
          end if
       end if
       
!Recoil atom
       append = 7
       yield = 0.0_doubel
       description='Recoil Atom'
       call totalYield(iso,'AR',yield)
       call findZnumber(trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name)), 'AR',zaid)
       if (yield > 0.0_doubel) then
          if (j == 1) write(unit_output,'(A7,A5,A5,ES14.7,A5,I2,A8,A)') '       ',trim(adjustl(zaid)),'     ',yield,'     ',20+append,'     $ ',trim(adjustl(description))
          if (j == 2) then
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A2,A10,A2,A20)') 'c ',trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name)),' ',trim(adjustl(description))
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A10,I2,A17,I2)') '#       SI',20+append,'               SP',20+append
             write(unit_output,*) '       L                  D'
!             do i = radfiledata%iso_rad_data(iso)%num_entries,1,-1
!                if (('AR' == trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(i)))).and.(radfiledata%iso_rad_data(iso)%em_line(i) > e_cutoff)) write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',radfiledata%iso_rad_data(iso)%em_line(i),'       ', radfiledata%iso_rad_data(iso)%em_yield(i)
!             end do

             i = radfiledata%iso_rad_data(iso)%num_entries
             do while (i > 0) 
                tempyield = 0.0_doubel
                counter = 1
                do k = i,1,-1
                   if ((radfiledata%iso_rad_data(iso)%em_line(k) > e_cutoff).and.(radfiledata%iso_rad_data(iso)%em_line(k)==radfiledata%iso_rad_data(iso)%em_line(i))) then
                      if ('AR' == trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(i)))) tempyield = tempyield + radfiledata%iso_rad_data(iso)%em_yield(k)
                      counter = counter + 1
                   else
                      exit
                   end if
                end do
                if (tempyield > 0.0_doubel) write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',radfiledata%iso_rad_data(iso)%em_line(i),'       ', tempyield
                if (counter == 1) then
                   i=i - (counter)
                else
                   i=i - (counter - 1)
                end if


             end do
             
          end if
       end if
       
!HFF 
       append = 8
       yield = 0.0_doubel
       description='H Fission Frag'
       call totalYield(iso,'FF',yield)
       call findZnumber(trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name)), 'HFF',zaid)
       if (yield > 0.0_doubel)  then
          if (j == 1) write(unit_output,'(A7,A5,A5,ES14.7,A5,I2,A8,A)') '       ',trim(adjustl(zaid)),'     ',yield,'     ',20+append,'     $ ',trim(adjustl(description))
          if (j == 2) then
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A2,A10,A2,A20)') 'c ',trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name)),' ',trim(adjustl(description))
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A10,I2,A17,I2)') '#       SI',20+append,'               SP',20+append
             write(unit_output,*) '       L                  D'
!             do i = radfiledata%iso_rad_data(iso)%num_entries,1,-1
!                if (('FF' == trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(i)))).and.(radfiledata%iso_rad_data(iso)%em_line(i) > e_cutoff)) write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',radfiledata%iso_rad_data(iso)%em_line(i),'       ', radfiledata%iso_rad_data(iso)%em_yield(i)
!             end do
             i = radfiledata%iso_rad_data(iso)%num_entries
             do while (i > 0) 
                tempyield = 0.0_doubel
                counter = 1
                do k = i,1,-1
                   if ((radfiledata%iso_rad_data(iso)%em_line(k) > e_cutoff).and.(radfiledata%iso_rad_data(iso)%em_line(k)==radfiledata%iso_rad_data(iso)%em_line(i))) then
                      if ('FF' == trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(i)))) tempyield = tempyield + radfiledata%iso_rad_data(iso)%em_yield(k)
                      counter = counter + 1
                   else
                      exit
                   end if
                end do
                if (tempyield > 0.0_doubel) write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',radfiledata%iso_rad_data(iso)%em_line(i),'       ', tempyield
                if (counter == 1) then
                   i=i - (counter)
                else
                   i=i - (counter - 1)
                end if


             end do

             
          end if
       end if
       
!HFF 
       append = 9
       yield = 0.0_doubel
       description='L Fission Frag'
       call totalYield(iso,'FF',yield)
       call findZnumber(trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name)), 'LFF',zaid)
       if (yield > 0.0_doubel) then
          if (j == 1) write(unit_output,'(A7,A5,A5,ES14.7,A5,I2,A8,A)') '       ',trim(adjustl(zaid)),'     ',yield,'     ',20+append,'     $ ',trim(adjustl(description))
          if (j == 2) then
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A2,A10,A2,A20)') 'c ',trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name)),' ',trim(adjustl(description))
             write(unit_output,'(A)') 'c'
             write(unit_output,'(A10,I2,A17,I2)') '#       SI',20+append,'               SP',20+append
             write(unit_output,*) '       L                  D'
!             do i = radfiledata%iso_rad_data(iso)%num_entries,1,-1
!                if (('FF' == trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(i)))).and.(radfiledata%iso_rad_data(iso)%em_line(i) > e_cutoff)) write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',radfiledata%iso_rad_data(iso)%em_line(i),'       ', radfiledata%iso_rad_data(iso)%em_yield(i)
!             end do
             i = radfiledata%iso_rad_data(iso)%num_entries
             do while (i > 0) 
                tempyield = 0.0_doubel
                counter = 1
                do k = i,1,-1
                   if ((radfiledata%iso_rad_data(iso)%em_line(k) > e_cutoff).and.(radfiledata%iso_rad_data(iso)%em_line(k)==radfiledata%iso_rad_data(iso)%em_line(i))) then
                      if ('FF' == trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(i)))) tempyield = tempyield + radfiledata%iso_rad_data(iso)%em_yield(k)
                      counter = counter + 1
                   else
                      exit
                   end if
                end do
                if (tempyield > 0.0_doubel) write(unit_output,'(A8,ES13.7,A7,ES13.7)') '        ',radfiledata%iso_rad_data(iso)%em_line(i),'       ', tempyield
                if (counter == 1) then
                   i=i - (counter)
                else
                   i=i - (counter - 1)
                end if


             end do

             
          end if
       end if
       
    end do
    
    write(unit_output,'(A)') 'c'
    write(unit_output,'(A18,A10)') 'c END OF SDEF FOR ',trim(adjustl(radfiledata%iso_rad_data(iso)%iso_name))
    write(unit_output,'(A)') 'c'
    
    call closeFile(unit_output)
  end subroutine writeMCNPSource
!-------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------!


    
!=====================================================================================!
  subroutine  totalYield (iso,desig,yield)
!=====================================================================================!
! Routine dumps MCNP sourcedef card in output file
!-------------------------------------------------------------------------------------!
!
! Revision record:
!--------------------------------------------------------------------------------------
! No.   Date          Programmer     Description of changes
!--------------------------------------------------------------------------------------
! 0     2013/11/11    R H Prinsloo    Original
!
    implicit none
! Dummy variables
! Local variables
    integer                            :: i,j,k                ! Array indices
    integer, INTENT(IN)                :: iso                  ! isotope
    character(len=*), INTENT(IN)       :: desig                ! particle designator
    real(kind = doubel), INTENT(INOUT) :: yield                ! total yield
!-------------------------------------------------------------------------------------!

    do i=1,radfiledata%iso_rad_data(iso)%num_entries
!       write(*,*) 'Comparing ',trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(i))),trim(adjustl(desig))
       if (trim(adjustl(radfiledata%iso_rad_data(iso)%em_part(i))) == trim(adjustl(desig))) then
          if (radfiledata%iso_rad_data(iso)%em_line(i) > e_cutoff) yield = yield + radfiledata%iso_rad_data(iso)%em_yield(i)
!          write(901,*) 'Found: ',yield
       end if

    end do



end subroutine totalYield
!-------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------!




!=====================================================================================!
  subroutine  findZnumber (iso,part,zaid)
!=====================================================================================!
! Routine dumps MCNP sourcedef card in output file
!-------------------------------------------------------------------------------------!
!
! Revision record:
!--------------------------------------------------------------------------------------
! No.   Date          Programmer     Description of changes
!--------------------------------------------------------------------------------------
! 0     2013/11/11    R H Prinsloo    Original
!
    implicit none
! Dummy variables
! Local variables
    integer                            :: i,j,k,l                ! Array indices
    character(len=*), INTENT(IN)       :: iso                  ! isotope
    character(len=*), INTENT(IN)       :: part                 ! particle type
    character(len=*), INTENT(OUT)      :: zaid                 ! particle designator
    character(len=255)                 :: temp,temp2           ! temp string
    character(len=255)                 :: temp3,temp4           ! temp string
!-------------------------------------------------------------------------------------!
    
    temp2 = iso
    temp=temp2(1:index(trim(adjustl(temp2)),'-')-1)
    temp3=temp2(index(trim(adjustl(temp2)),'-')+1:len(trim(adjustl(temp2))))
    if ((index(trim(adjustl(temp3)),'m') > 0).or.((index(trim(adjustl(temp3)),'n') > 0))) temp3 = temp3(1:len(trim(adjustl(temp3)))-1)
!    write(*,*) 'Looking for: ',trim(adjustl(temp)),' ',trim(adjustl(temp3))
    
    zaid='XX'
    do i=1,num_elements
!       write(*,*) 'Comparing: ',trim(adjustl(temp)),' ',trim(adjustl(element_data(i)%element_name))
       if (trim(adjustl(element_data(i)%element_name)) == trim(adjustl(temp))) then
          if (trim(adjustl(part)) == 'AR') write(zaid,*) 1000 * (element_data(i)%z-2) + (CHINT(trim(adjustl(temp3)))-4)
          if (trim(adjustl(part)) == 'LFF') write(zaid,*) 1000 * (ceiling(0.4 * element_data(i)%z)) + ceiling(((CHINT(trim(adjustl(temp3)))-4) * 0.4_doubel))
          if (trim(adjustl(part)) == 'HFF') write(zaid,*) 1000 * (ceiling(0.6 * element_data(i)%z)) + ceiling(((CHINT(trim(adjustl(temp3)))-4) * 0.6_doubel))
          exit 
       end if
    end do



end subroutine findZNumber
!-------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------!

 !=====================================================================================!
  integer function CHINT (ch)
    !=====================================================================================!
    ! Function that transforms character variable values that are a subset of the string
    ! '012345679' into integer values.
    !-------------------------------------------------------------------------------------!
    !
    ! Revision record:
    !--------------------------------------------------------------------------------------
    ! No.   Date          Programmer     Description of changes
    !--------------------------------------------------------------------------------------
    ! 0     2008/04/11    G J Stander    Original
    ! 1     2008/08/22    G J Stander    Add negative integer processing
    !
    implicit none
    ! Dummy variables
    character(LEN=*), intent(IN) :: ch                  ! Character input
    ! Local variables
    character(LEN= &
         len(trim(adjustl(ch)))) :: lch                 ! Local ch
    integer                      :: i                   ! Array index
    integer                      :: n                   ! Power index
    character(LEN=10), parameter :: seq = '0123456789'  ! Sequence of subset of natural numbers
    logical                      :: negative            ! Flag to indicate negative number
    character(LEN=*), parameter  :: srtname = 'CHINT'   ! Routine name
    !-------------------------------------------------------------------------------------!

    CHINT = 0                                           ! Initialize CHINT
    negative = .false.                                  ! Initialize flag
    lch = trim(adjustl(ch))                             ! Assign lch

    ! If negative number
    if (index(lch,'-') == 1) then
       ! Multiply by 10
       lch      = trim(lch(2:))//'0'
       negative = .true.
    end if

    n = 0                                               ! Initialize n

    ! Loop is from last character of character string ch to first character
    do i = len(trim(lch)),1,-1

       ! If character from ch is not in seq then error message is generated
       if (index(seq,lch(i:i)) /= 0) then
          ! Index of character after other operations added to integer output
          CHINT = CHINT + (index(seq,lch(i:i))-1)*10**n
       else
          write(*,*) 'Error in integer conversion...',trim(adjustl(ch))   
       end if
       ! Raising of power n by 1
       n = n + 1
    end do

    ! If negative number
    if (negative) then
       CHINT = -int(real(CHINT,KIND=doubel) / 10.0_doubel)
    end if

  end function CHINT


End Program mcnpsource


