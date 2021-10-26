program ufsLandNoahMPReplaceRestart

  implicit none

  type noahmp_type
    double precision, allocatable :: swe                (:)
    double precision, allocatable :: snow_depth         (:)
    double precision, allocatable :: active_snow_layers (:)
    double precision, allocatable :: swe_previous       (:)
    double precision, allocatable :: snow_soil_interface(:,:)
    double precision, allocatable :: temperature_snow   (:,:)
    double precision, allocatable :: snow_ice_layer     (:,:)
    double precision, allocatable :: snow_liq_layer     (:,:)
    double precision, allocatable :: temperature_soil   (:)
  end type noahmp_type    

  type observation_type
    double precision, allocatable :: snow_depth (:)
  end type observation_type    

  type(noahmp_type)      :: noahmp
  type(observation_type) :: obs
  character*256          :: restart_filename
  character*256          :: obs_filename
  character*19           :: restart_date = ""
  character*128          :: restart_dir = ""
  integer                :: vector_length
  integer                :: yyyy,mm,dd,hh,nn,ss
  integer                :: option = -1
  character*32           :: arg = "blah"
  
  integer, parameter :: DO_NOTHING = 0
  integer, parameter :: TOP_LAYER = 1
  integer, parameter :: BOTTOM_LAYER = 2
  integer, parameter :: ALL_LAYERS = 3
  
  call get_command_argument(1, arg)
  if(arg == "") stop "add a valid option to the command line"
  read(arg,'(i1)') option
  
  call ReadNamelist(restart_dir, restart_date)

  read(restart_date( 1: 4),'(i4.4)') yyyy
  read(restart_date( 6: 7),'(i2.2)') mm
  read(restart_date( 9:10),'(i2.2)') dd
  read(restart_date(12:13),'(i2.2)') hh
  read(restart_date(15:16),'(i2.2)') nn
  read(restart_date(18:19),'(i2.2)') ss

  write(restart_filename,'(a17,i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a3)') &
    "ufs_land_restart.", yyyy, "-", mm, "-", dd, "_", hh, "-", nn, "-", ss, ".nc"

  restart_filename = trim(restart_dir)//trim(restart_filename)
  
  write(obs_filename,'(a22,i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') &
    "fake_snow/fake.sndpth.", yyyy, "-", mm, "-", dd, "_", hh, "-", nn, "-", ss
  
  call ReadVectorLength(restart_filename, vector_length)
  
  allocate(noahmp%swe                (vector_length))
  allocate(noahmp%snow_depth         (vector_length))
  allocate(noahmp%active_snow_layers (vector_length))
  allocate(noahmp%swe_previous       (vector_length))
  allocate(noahmp%snow_soil_interface(vector_length,7))
  allocate(noahmp%temperature_snow   (vector_length,3))
  allocate(noahmp%snow_ice_layer     (vector_length,3))
  allocate(noahmp%snow_liq_layer     (vector_length,3))
  allocate(noahmp%temperature_soil   (vector_length))

  allocate(obs%snow_depth(vector_length))

  call ReadObservation(obs_filename, vector_length, obs)
  
  call ReadRestartNoahMP(restart_filename, vector_length, noahmp)
  
  select case (option)
  
    case (DO_NOTHING)

      write(*,*) "Option 0: do nothing"
    
    case (TOP_LAYER)

      write(*,*) "Option 1: add/remove all snow from top"
      call UpdateTopLayer(vector_length, noahmp, obs)
    
    case (BOTTOM_LAYER)

      write(*,*) "Option 2: add/remove all snow from bottom"
      call UpdateBottomLayer(vector_length, noahmp, obs)
    
    case (ALL_LAYERS)

      write(*,*) "Option 3: distribute snow in all layers"
      call UpdateAllLayers(vector_length, noahmp, obs)
    
    case default
    
      write(*,*) "choose a valid partition option"
  
  end select 

  call WriteRestartNoahMP(restart_filename, vector_length, noahmp)
     
contains   

  subroutine UpdateTopLayer(vector_length, noahmp, obs)
  
  type(noahmp_type)      :: noahmp
  type(observation_type) :: obs
  integer                :: vector_length
  double precision       :: increment(vector_length), to_remove
  double precision       :: layer_density, swe_increment, liq_ratio, remove_ratio
  integer                :: iloc, ilayer, iinter, active_layers, vector_loc, pathway
  double precision       :: soil_interfaces(7) = (/0.0,0.0,0.0,0.1,0.4,1.0,2.0/)
  
  associate( &
      obs_snow_depth => obs%snow_depth            ,&
                 swe => noahmp%swe                ,&
          snow_depth => noahmp%snow_depth         ,&
  active_snow_layers => noahmp%active_snow_layers ,&
        swe_previous => noahmp%swe_previous       ,&
 snow_soil_interface => noahmp%snow_soil_interface,&
    temperature_snow => noahmp%temperature_snow   ,&
      snow_ice_layer => noahmp%snow_ice_layer     ,&
      snow_liq_layer => noahmp%snow_liq_layer     ,&
    temperature_soil => noahmp%temperature_soil )

  
  increment = obs_snow_depth - snow_depth  ! snow to add or remove [mm]
  
  do iloc = 1, vector_length
    
    pathway = 0
    
    if(obs_snow_depth(iloc) == 0.0) then

      swe                (iloc)   = 0.0
      snow_depth         (iloc)   = 0.0
      active_snow_layers (iloc)   = 0.0
      swe_previous       (iloc)   = 0.0
      snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
      temperature_snow   (iloc,:) = 0.0
      snow_ice_layer     (iloc,:) = 0.0
      snow_liq_layer     (iloc,:) = 0.0

    else
    
      active_layers = nint(active_snow_layers(iloc))  ! number of active layers (0,-1,-2,-3)

      if(active_layers < 0) then  ! in multi-layer mode

        if(increment(iloc) > 0.0) then  ! add snow in multi-layer mode

          pathway = 1

          vector_loc = 4 + active_layers

          layer_density = (snow_ice_layer(iloc,vector_loc)+snow_liq_layer(iloc,vector_loc)) / &
                            (-snow_soil_interface(iloc,vector_loc))
          swe_increment = increment(iloc) * layer_density / 1000.d0
          liq_ratio = snow_liq_layer(iloc,vector_loc) / &
                        ( snow_ice_layer(iloc,vector_loc) + snow_liq_layer(iloc,vector_loc) )
          snow_ice_layer(iloc,vector_loc) = snow_ice_layer(iloc,vector_loc) + &
                                              (1.0 - liq_ratio) * swe_increment
          snow_liq_layer(iloc,vector_loc) = snow_liq_layer(iloc,vector_loc) + &
                                              liq_ratio * swe_increment
          do ilayer = vector_loc, 3
            snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,ilayer) - increment(iloc)/1000.d0
          end do
          
        elseif(increment(iloc) < 0.0) then  ! remove snow in multi-layer mode

          pathway = 2
    
          to_remove = increment(iloc)  ! depth[mm] to be removed

          vector_loc = 4 + active_layers  ! location in vector of top layer
          
          layerloop: do ilayer = vector_loc, 3
          
            if(to_remove < 1000.d0*snow_soil_interface(iloc,ilayer)) then  ! this entire layer will be removed
              
              to_remove = to_remove - 1000.d0*snow_soil_interface(iloc,ilayer)
              snow_ice_layer(iloc,ilayer)      = 0.0
              snow_liq_layer(iloc,ilayer)      = 0.0
              temperature_snow(iloc,ilayer)    = 0.0
              do iinter = 3,ilayer, -1  ! remove snow from each snow layer
                snow_soil_interface(iloc,iinter) = snow_soil_interface(iloc,iinter) - snow_soil_interface(iloc,ilayer)
              end do
              
            else  ! this layer will be partially removed
              
              active_snow_layers(iloc) = ilayer - 4  ! new number of layers
              remove_ratio = 1.d0 - (to_remove/1000.d0/snow_soil_interface(iloc,ilayer)) ! fraction to remove from layer
              snow_ice_layer(iloc,ilayer) = remove_ratio * snow_ice_layer(iloc,ilayer)
              snow_liq_layer(iloc,ilayer) = remove_ratio * snow_liq_layer(iloc,ilayer)
              do iinter = ilayer, 3  ! remove snow from each snow layer
                snow_soil_interface(iloc,iinter) = snow_soil_interface(iloc,iinter) - to_remove/1000.d0
              end do
            
              exit layerloop

            end if

          end do layerloop
          
        end if  ! increment
        
        ! For multi-layer mode, recalculate interfaces and sum depth/swe

        do ilayer = 4, 7
          snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,3) - soil_interfaces(ilayer)
        end do

        snow_depth(iloc) = -snow_soil_interface(iloc,3) * 1000.d0

        swe(iloc) = 0.0

        do ilayer = 1, 3
          swe(iloc) = swe(iloc) + snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer)
        end do

        swe_previous(iloc) = swe(iloc)

        if(snow_depth(iloc) < 25.d0) then  ! go out of multi-layer mode
          active_snow_layers (iloc) = 0.d0
          snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
          temperature_snow   (iloc,:) = 0.0
          snow_ice_layer     (iloc,:) = 0.0
          snow_liq_layer     (iloc,:) = 0.0
        end if

      elseif(active_layers == 0) then  ! snow starts in zero-layer mode

        if(increment(iloc) > 0.0) then  ! add snow in zero-layer mode
    
          if(snow_depth(iloc) == 0) then   ! no snow present, so assume density based on soil temperature
            pathway = 3
            layer_density = max(80.0,min(120.,67.92+51.25*exp((temperature_soil(iloc)-273.15)/2.59)))
          else   ! use existing density
            pathway = 4
            layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
          end if
          snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
          swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
          swe_previous(iloc) = swe(iloc)

          active_snow_layers(iloc)      = 0.0
          snow_ice_layer(iloc,:)        = 0.0
          snow_liq_layer(iloc,:)        = 0.0
          temperature_snow(iloc,:)      = 0.0
          snow_soil_interface(iloc,1:3) = 0.0

          if(snow_depth(iloc) > 25.0) then  ! snow depth is > 25mm so put in a layer
            pathway = 5
            active_snow_layers(iloc) = -1.0
            snow_ice_layer(iloc,3)   = swe(iloc)
            temperature_snow(iloc,3) = temperature_soil(iloc)
            do ilayer = 3, 7
              snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,ilayer) - snow_depth(iloc)/1000.d0
            end do
          end if
          
        elseif(increment(iloc) < 0.0) then  ! remove snow in zero-layer mode

          pathway = 6
    
          if(snow_depth(iloc) <= 0.0) stop "inconsistency in snow_depth and increment"
          layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
          snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
          swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
          swe_previous(iloc) = swe(iloc)

          active_snow_layers(iloc)      = 0.0
          snow_ice_layer(iloc,:)        = 0.0
          snow_liq_layer(iloc,:)        = 0.0
          temperature_snow(iloc,:)      = 0.0
          snow_soil_interface(iloc,1:3) = 0.0

        end if  ! increment

      end if  ! active_layers
        
    end if  ! obs_snow_depth == 0
    
    ! do some gross checks

    if(abs(snow_soil_interface(iloc,7) - snow_soil_interface(iloc,3) + 2.d0) > 0.0000001) then
      print*, "Depth of soil not 2m"
      print*, pathway
      print*, snow_soil_interface(iloc,7), snow_soil_interface(iloc,3)
!      stop
    end if

    if(active_snow_layers(iloc) < 0.0 .and. abs(snow_depth(iloc) + 1000.d0*snow_soil_interface(iloc,3)) > 0.0000001) then
      print*, "snow_depth and snow_soil_interface inconsistent"
      print*, pathway
      print*, active_snow_layers(iloc), snow_depth(iloc), snow_soil_interface(iloc,3)
!      stop
    end if

    if(abs(obs_snow_depth(iloc) - snow_depth(iloc)) > 0.0000001) then
      print*, "observed snow and updated model snow inconsistent"
      print*, pathway
      print*, obs_snow_depth(iloc), snow_depth(iloc), snow_soil_interface(iloc,3)
!      stop
    end if

    if(snow_depth(iloc) < 0.0 .or. snow_soil_interface(iloc,3) > 0.0 ) then
      print*, "observed snow and updated model snow inconsistent"
      print*, pathway
      print*, snow_depth(iloc), snow_soil_interface(iloc,3)
!      stop
    end if

  end do
  
  end associate
   
  end subroutine UpdateTopLayer
  
  subroutine UpdateBottomLayer(vector_length, noahmp, obs)
  
  type(noahmp_type)      :: noahmp
  type(observation_type) :: obs
  integer                :: vector_length
  double precision       :: increment(vector_length), to_remove
  double precision       :: layer_density, swe_increment, liq_ratio, remove_ratio
  integer                :: iloc, ilayer, iinter, active_layers, vector_loc, pathway, removed
  double precision       :: soil_interfaces(7) = (/0.0,0.0,0.0,0.1,0.4,1.0,2.0/)
  double precision       :: temp_vector(3), layer_depths(3)
  
  associate( &
      obs_snow_depth => obs%snow_depth            ,&
                 swe => noahmp%swe                ,&
          snow_depth => noahmp%snow_depth         ,&
  active_snow_layers => noahmp%active_snow_layers ,&
        swe_previous => noahmp%swe_previous       ,&
 snow_soil_interface => noahmp%snow_soil_interface,&
    temperature_snow => noahmp%temperature_snow   ,&
      snow_ice_layer => noahmp%snow_ice_layer     ,&
      snow_liq_layer => noahmp%snow_liq_layer     ,&
    temperature_soil => noahmp%temperature_soil )

  
  increment = obs_snow_depth - snow_depth  ! snow to add or remove [mm]
  
  do iloc = 1, vector_length
    
    pathway = 0
    
    if(obs_snow_depth(iloc) == 0.0) then

      swe                (iloc)   = 0.0
      snow_depth         (iloc)   = 0.0
      active_snow_layers (iloc)   = 0.0
      swe_previous       (iloc)   = 0.0
      snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
      temperature_snow   (iloc,:) = 0.0
      snow_ice_layer     (iloc,:) = 0.0
      snow_liq_layer     (iloc,:) = 0.0

    else
    
      active_layers = nint(active_snow_layers(iloc))  ! number of active layers (0,-1,-2,-3)

      if(active_layers < 0) then  ! in multi-layer mode

        if(increment(iloc) > 0.0) then  ! add snow in multi-layer mode

          pathway = 1
    
          vector_loc = 3

          layer_density = (snow_ice_layer(iloc,vector_loc)+snow_liq_layer(iloc,vector_loc)) / &
                            (-snow_soil_interface(iloc,vector_loc))
          swe_increment = increment(iloc) * layer_density / 1000.d0
          liq_ratio = snow_liq_layer(iloc,vector_loc) / &
                        ( snow_ice_layer(iloc,vector_loc) + snow_liq_layer(iloc,vector_loc) )
          snow_ice_layer(iloc,vector_loc) = snow_ice_layer(iloc,vector_loc) + &
                                              (1.0 - liq_ratio) * swe_increment
          snow_liq_layer(iloc,vector_loc) = snow_liq_layer(iloc,vector_loc) + &
                                              liq_ratio * swe_increment
          snow_soil_interface(iloc,3) = snow_soil_interface(iloc,3) - increment(iloc)/1000.d0
          
        elseif(increment(iloc) < 0.0) then  ! remove snow in multi-layer mode

          pathway = 2
          
          layer_depths(1) = snow_soil_interface(iloc,1)  ! layer depth [m] (negative)
          layer_depths(2) = snow_soil_interface(iloc,2)-snow_soil_interface(iloc,1)
          layer_depths(3) = snow_soil_interface(iloc,3)-snow_soil_interface(iloc,2)

          removed = 0
    
          to_remove = increment(iloc)  ! depth[mm] to be removed

          vector_loc = 4 + active_layers  ! location in vector of top layer
          
          layerloop: do ilayer = 3, vector_loc, -1
          
            if(to_remove < 1000.d0*layer_depths(ilayer)) then  ! this entire layer will be removed
              
              to_remove = to_remove - 1000.d0*layer_depths(ilayer)
              removed = removed + 1

              snow_ice_layer(iloc,ilayer)      = 0.0
              snow_liq_layer(iloc,ilayer)      = 0.0
              temperature_snow(iloc,ilayer)    = 0.0
              snow_soil_interface(iloc,ilayer) = 0.0

            else  ! this layer will be partially removed
              
              remove_ratio = 1.d0 - (to_remove/1000.d0/snow_soil_interface(iloc,ilayer)) ! fraction to remove from layer
              snow_ice_layer(iloc,ilayer) = remove_ratio * snow_ice_layer(iloc,ilayer)
              snow_liq_layer(iloc,ilayer) = remove_ratio * snow_liq_layer(iloc,ilayer)
              snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,ilayer) - to_remove/1000.d0

              exit layerloop

            end if

          end do layerloop
          
          active_snow_layers(iloc) = active_snow_layers(iloc) + removed  ! new number of layers
          
          temp_vector = temperature_snow(iloc,:)
          temperature_snow(iloc,:) = 0.0
          temperature_snow(iloc,removed+1:3) = temp_vector(1:3-removed)

          temp_vector = snow_ice_layer(iloc,:)
          snow_ice_layer(iloc,:) = 0.0
          snow_ice_layer(iloc,removed+1:3) = temp_vector(1:3-removed)

          temp_vector = snow_liq_layer(iloc,:)
          snow_liq_layer(iloc,:) = 0.0
          snow_liq_layer(iloc,removed+1:3) = temp_vector(1:3-removed)

          temp_vector = snow_soil_interface(iloc,:)
          snow_soil_interface(iloc,:) = 0.0
          snow_soil_interface(iloc,removed+1:3) = temp_vector(1:3-removed)
          
        end if  ! increment
        
        ! For multi-layer mode, recalculate interfaces and sum depth/swe

        do ilayer = 4, 7
          snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,3) - soil_interfaces(ilayer)
        end do

        snow_depth(iloc) = -snow_soil_interface(iloc,3) * 1000.d0

        swe(iloc) = 0.0

        do ilayer = 1, 3
          swe(iloc) = swe(iloc) + snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer)
        end do

        swe_previous(iloc) = swe(iloc)
        
        if(snow_depth(iloc) < 25.d0) then  ! go out of multi-layer mode
          active_snow_layers (iloc) = 0.d0
          snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
          temperature_snow   (iloc,:) = 0.0
          snow_ice_layer     (iloc,:) = 0.0
          snow_liq_layer     (iloc,:) = 0.0
        end if

      elseif(active_layers == 0) then  ! snow starts in zero-layer mode

        if(increment(iloc) > 0.0) then  ! add snow in zero-layer mode
    
          if(snow_depth(iloc) == 0) then   ! no snow present, so assume density based on soil temperature
            pathway = 3
            layer_density = max(80.0,min(120.,67.92+51.25*exp((temperature_soil(iloc)-273.15)/2.59)))
          else   ! use existing density
            pathway = 4
            layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
          end if
          snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
          swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
          swe_previous(iloc) = swe(iloc)

          active_snow_layers(iloc)      = 0.0
          snow_ice_layer(iloc,:)        = 0.0
          snow_liq_layer(iloc,:)        = 0.0
          temperature_snow(iloc,:)      = 0.0
          snow_soil_interface(iloc,1:3) = 0.0

          if(snow_depth(iloc) > 25.0) then  ! snow depth is > 25mm so put in a layer
            pathway = 5
            active_snow_layers(iloc) = -1.0
            snow_ice_layer(iloc,3)   = swe(iloc)
            temperature_snow(iloc,3) = temperature_soil(iloc)
            do ilayer = 3, 7
              snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,ilayer) - snow_depth(iloc)/1000.d0
            end do
          end if
          
        elseif(increment(iloc) < 0.0) then  ! remove snow in zero-layer mode

          pathway = 6
    
          if(snow_depth(iloc) <= 0.0) stop "inconsistency in snow_depth and increment"
          layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
          snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
          swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
          swe_previous(iloc) = swe(iloc)

          active_snow_layers(iloc)      = 0.0
          snow_ice_layer(iloc,:)        = 0.0
          snow_liq_layer(iloc,:)        = 0.0
          temperature_snow(iloc,:)      = 0.0
          snow_soil_interface(iloc,1:3) = 0.0

        end if  ! increment

      end if  ! active_layers
        
    end if  ! obs_snow_depth == 0
    
    ! do some gross checks

    if(abs(snow_soil_interface(iloc,7) - snow_soil_interface(iloc,3) + 2.d0) > 0.0000001) then
      print*, "Depth of soil not 2m"
      print*, pathway
      print*, snow_soil_interface(iloc,7), snow_soil_interface(iloc,3)
!      stop
    end if

    if(active_snow_layers(iloc) < 0.0 .and. abs(snow_depth(iloc) + 1000.d0*snow_soil_interface(iloc,3)) > 0.0000001) then
      print*, "snow_depth and snow_soil_interface inconsistent"
      print*, pathway
      print*, active_snow_layers(iloc), snow_depth(iloc), snow_soil_interface(iloc,3)
!      stop
    end if

    if(abs(obs_snow_depth(iloc) - snow_depth(iloc)) > 0.0000001) then
      print*, "observed snow and updated model snow inconsistent"
      print*, pathway
      print*, obs_snow_depth(iloc), snow_depth(iloc)
!      stop
    end if

    if(snow_depth(iloc) < 0.0 .or. snow_soil_interface(iloc,3) > 0.0 ) then
      print*, "observed snow and updated model snow inconsistent"
      print*, pathway
      print*, snow_depth(iloc), snow_soil_interface(iloc,3)
!      stop
    end if

  end do
  
  end associate
   
  end subroutine UpdateBottomLayer
  
  subroutine UpdateAllLayers(vector_length, noahmp, obs)
  
  type(noahmp_type)      :: noahmp
  type(observation_type) :: obs
  integer                :: vector_length
  double precision       :: increment(vector_length), to_remove
  double precision       :: layer_density, swe_increment, liq_ratio, remove_ratio
  integer                :: iloc, ilayer, iinter, active_layers, vector_loc, pathway, removed
  double precision       :: soil_interfaces(7) = (/0.0,0.0,0.0,0.1,0.4,1.0,2.0/)
  double precision       :: partition_ratio, layer_depths(3)
  
  associate( &
      obs_snow_depth => obs%snow_depth            ,&
                 swe => noahmp%swe                ,&
          snow_depth => noahmp%snow_depth         ,&
  active_snow_layers => noahmp%active_snow_layers ,&
        swe_previous => noahmp%swe_previous       ,&
 snow_soil_interface => noahmp%snow_soil_interface,&
    temperature_snow => noahmp%temperature_snow   ,&
      snow_ice_layer => noahmp%snow_ice_layer     ,&
      snow_liq_layer => noahmp%snow_liq_layer     ,&
    temperature_soil => noahmp%temperature_soil )

  
  increment = obs_snow_depth - snow_depth  ! snow to add or remove [mm]
  
  do iloc = 1, vector_length
    
    pathway = 0
    
    if(obs_snow_depth(iloc) == 0.0) then

      swe                (iloc)   = 0.0
      snow_depth         (iloc)   = 0.0
      active_snow_layers (iloc)   = 0.0
      swe_previous       (iloc)   = 0.0
      snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
      temperature_snow   (iloc,:) = 0.0
      snow_ice_layer     (iloc,:) = 0.0
      snow_liq_layer     (iloc,:) = 0.0

    else
    
      active_layers = nint(active_snow_layers(iloc))  ! number of active layers (0,-1,-2,-3)

      if(active_layers < 0) then  ! in multi-layer mode
      
        layer_depths(1) = snow_soil_interface(iloc,1)
        layer_depths(2) = snow_soil_interface(iloc,2)-snow_soil_interface(iloc,1)
        layer_depths(3) = snow_soil_interface(iloc,3)-snow_soil_interface(iloc,2)

        if(increment(iloc) > 0.0) then  ! add snow in multi-layer mode

          pathway = 1
    
          vector_loc = 4 + active_layers  ! location in vector of top layer
          
          layerloop: do ilayer = vector_loc, 3
          
            partition_ratio = -layer_depths(ilayer)/snow_depth(iloc)*1000.d0
            layer_density = (snow_ice_layer(iloc,ilayer)+snow_liq_layer(iloc,ilayer)) / &
                              (-layer_depths(ilayer))
            swe_increment = partition_ratio * increment(iloc) * layer_density / 1000.d0
            liq_ratio = snow_liq_layer(iloc,ilayer) / &
                          ( snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer) )
            snow_ice_layer(iloc,ilayer) = snow_ice_layer(iloc,ilayer) + &
                                              (1.0 - liq_ratio) * swe_increment
            snow_liq_layer(iloc,ilayer) = snow_liq_layer(iloc,ilayer) + &
                                              liq_ratio * swe_increment
            do iinter = ilayer, 3  ! remove snow from each snow layer
              snow_soil_interface(iloc,iinter) = snow_soil_interface(iloc,iinter) - &
                                                   partition_ratio * increment(iloc)/1000.d0
            end do

          end do layerloop
            
        elseif(increment(iloc) < 0.0) then  ! remove snow in multi-layer mode

          pathway = 2
          
          vector_loc = 4 + active_layers  ! location in vector of top layer
          
          layerloop: do ilayer = vector_loc, 3
          
            partition_ratio = -layer_depths(ilayer)/snow_depth(iloc)*1000.d0
            layer_density = (snow_ice_layer(iloc,ilayer)+snow_liq_layer(iloc,ilayer)) / &
                              (-layer_depths(ilayer))
            swe_increment = partition_ratio * increment(iloc) * layer_density / 1000.d0
            liq_ratio = snow_liq_layer(iloc,ilayer) / &
                          ( snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer) )
            snow_ice_layer(iloc,ilayer) = snow_ice_layer(iloc,ilayer) + &
                                              (1.0 - liq_ratio) * swe_increment
            snow_liq_layer(iloc,ilayer) = snow_liq_layer(iloc,ilayer) + &
                                              liq_ratio * swe_increment
            do iinter = ilayer, 3  ! remove snow from each snow layer
              snow_soil_interface(iloc,iinter) = snow_soil_interface(iloc,iinter) - &
                                                   partition_ratio * increment(iloc)/1000.d0
            end do

          end do layerloop
          
        end if  ! increment
        
        ! For multi-layer mode, recalculate interfaces and sum depth/swe

        do ilayer = 4, 7
          snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,3) - soil_interfaces(ilayer)
        end do

        snow_depth(iloc) = -snow_soil_interface(iloc,3) * 1000.d0

        swe(iloc) = 0.0

        do ilayer = 1, 3
          swe(iloc) = swe(iloc) + snow_ice_layer(iloc,ilayer) + snow_liq_layer(iloc,ilayer)
        end do

        swe_previous(iloc) = swe(iloc)

        if(snow_depth(iloc) < 25.d0) then  ! go out of multi-layer mode
          active_snow_layers (iloc) = 0.d0
          snow_soil_interface(iloc,:) = (/0.0,0.0,0.0,-0.1,-0.4,-1.0,-2.0/)
          temperature_snow   (iloc,:) = 0.0
          snow_ice_layer     (iloc,:) = 0.0
          snow_liq_layer     (iloc,:) = 0.0
        end if

      elseif(active_layers == 0) then  ! snow starts in zero-layer mode

        if(increment(iloc) > 0.0) then  ! add snow in zero-layer mode
    
          if(snow_depth(iloc) == 0) then   ! no snow present, so assume density based on soil temperature
            pathway = 3
            layer_density = max(80.0,min(120.,67.92+51.25*exp((temperature_soil(iloc)-273.15)/2.59)))
          else   ! use existing density
            pathway = 4
            layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
          end if
          snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
          swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
          swe_previous(iloc) = swe(iloc)

          active_snow_layers(iloc)      = 0.0
          snow_ice_layer(iloc,:)        = 0.0
          snow_liq_layer(iloc,:)        = 0.0
          temperature_snow(iloc,:)      = 0.0
          snow_soil_interface(iloc,1:3) = 0.0

          if(snow_depth(iloc) > 25.0) then  ! snow depth is > 25mm so put in a layer
            pathway = 5
            active_snow_layers(iloc) = -1.0
            snow_ice_layer(iloc,3)   = swe(iloc)
            temperature_snow(iloc,3) = temperature_soil(iloc)
            do ilayer = 3, 7
              snow_soil_interface(iloc,ilayer) = snow_soil_interface(iloc,ilayer) - snow_depth(iloc)/1000.d0
            end do
          end if
          
        elseif(increment(iloc) < 0.0) then  ! remove snow in zero-layer mode

          pathway = 6
    
          if(snow_depth(iloc) <= 0.0) stop "inconsistency in snow_depth and increment"
          layer_density = swe(iloc) / snow_depth(iloc) * 1000.d0
          snow_depth(iloc) = snow_depth(iloc) + increment(iloc)
          swe(iloc) = swe(iloc) + increment(iloc) * layer_density / 1000.d0
          swe_previous(iloc) = swe(iloc)

          active_snow_layers(iloc)      = 0.0
          snow_ice_layer(iloc,:)        = 0.0
          snow_liq_layer(iloc,:)        = 0.0
          temperature_snow(iloc,:)      = 0.0
          snow_soil_interface(iloc,1:3) = 0.0

        end if  ! increment

      end if  ! active_layers
        
    end if  ! obs_snow_depth == 0
    
    ! do some gross checks

    if(abs(snow_soil_interface(iloc,7) - snow_soil_interface(iloc,3) + 2.d0) > 0.0000001) then
      print*, "Depth of soil not 2m"
      print*, pathway
      print*, snow_soil_interface(iloc,7), snow_soil_interface(iloc,3)
!      stop
    end if

    if(active_snow_layers(iloc) < 0.0 .and. abs(snow_depth(iloc) + 1000.d0*snow_soil_interface(iloc,3)) > 0.0000001) then
      print*, "snow_depth and snow_soil_interface inconsistent"
      print*, pathway
      print*, active_snow_layers(iloc), snow_depth(iloc), snow_soil_interface(iloc,3)
!      stop
    end if

    if(abs(obs_snow_depth(iloc) - snow_depth(iloc)) > 0.0000001) then
      print*, "observed snow and updated model snow inconsistent"
      print*, pathway
      print*, obs_snow_depth(iloc), snow_depth(iloc)
!      stop
    end if

    if(snow_depth(iloc) < 0.0 .or. snow_soil_interface(iloc,3) > 0.0 ) then
      print*, "observed snow and updated model snow inconsistent"
      print*, pathway
      print*, snow_depth(iloc), snow_soil_interface(iloc,3)
!      stop
    end if

  end do
  
  end associate
   
  end subroutine UpdateAllLayers
  
  subroutine WriteRestartNoahMP(filename, vector_length, noahmp)
  
  use netcdf

  type(noahmp_type) :: noahmp
  character*256     :: filename
  integer           :: vector_length
  integer           :: ncid, varid, status
  
  status = nf90_open(filename, NF90_WRITE, ncid)

! Start writing restart file
  
  status = nf90_inq_varid(ncid, "weasd", varid)
  status = nf90_put_var(ncid, varid , noahmp%swe   , &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "snwdph", varid)
  status = nf90_put_var(ncid, varid , noahmp%snow_depth  , &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "snowxy", varid)
  status = nf90_put_var(ncid, varid , noahmp%active_snow_layers  , &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "sneqvoxy", varid)
  status = nf90_put_var(ncid, varid , noahmp%swe_previous, &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "zsnsoxy", varid)
  status = nf90_put_var(ncid, varid , noahmp%snow_soil_interface , &
      start = (/1            , 1, 1/), &
      count = (/vector_length, 7, 1/))

  status = nf90_inq_varid(ncid, "tsnoxy", varid)
  status = nf90_put_var(ncid, varid , noahmp%temperature_snow  , &
      start = (/1            , 1, 1/), &
      count = (/vector_length, 3, 1/))

  status = nf90_inq_varid(ncid, "snicexy", varid)
  status = nf90_put_var(ncid, varid , noahmp%snow_ice_layer , &
      start = (/1            , 1, 1/), &
      count = (/vector_length, 3, 1/))

  status = nf90_inq_varid(ncid, "snliqxy", varid)
  status = nf90_put_var(ncid, varid , noahmp%snow_liq_layer , &
      start = (/1            , 1, 1/), &
      count = (/vector_length, 3, 1/))

  status = nf90_close(ncid)

  end subroutine WriteRestartNoahMP
  
  subroutine ReadRestartNoahMP(filename, vector_length, noahmp)
  
  use netcdf

  type(noahmp_type) :: noahmp
  character*256     :: filename
  integer           :: vector_length
  integer           :: ncid, dimid, varid, status
  
  status = nf90_open(filename, NF90_NOWRITE, ncid)

  status = nf90_inq_varid(ncid, "weasd", varid)
  status = nf90_get_var(ncid, varid , noahmp%swe   , &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "snwdph", varid)
  status = nf90_get_var(ncid, varid , noahmp%snow_depth  , &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "snowxy", varid)
  status = nf90_get_var(ncid, varid , noahmp%active_snow_layers  , &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "sneqvoxy", varid)
  status = nf90_get_var(ncid, varid , noahmp%swe_previous, &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_inq_varid(ncid, "tsnoxy", varid)
  status = nf90_get_var(ncid, varid , noahmp%temperature_snow  , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 3, 1/))

  status = nf90_inq_varid(ncid, "zsnsoxy", varid)
  status = nf90_get_var(ncid, varid , noahmp%snow_soil_interface , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 7, 1/))

  status = nf90_inq_varid(ncid, "snicexy", varid)
  status = nf90_get_var(ncid, varid , noahmp%snow_ice_layer , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 3, 1/))

  status = nf90_inq_varid(ncid, "snliqxy", varid)
  status = nf90_get_var(ncid, varid , noahmp%snow_liq_layer , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 3, 1/))
 
  status = nf90_inq_varid(ncid, "stc", varid)
  status = nf90_get_var(ncid, varid , noahmp%temperature_soil , &
      start = (/1            , 1, 1/)                , &
      count = (/vector_length, 1, 1/))

  status = nf90_close(ncid)

  end subroutine ReadRestartNoahMP

  subroutine ReadObservation(filename, vector_length, observation)
  
  use netcdf

  type(observation_type) :: observation
  character*256          :: filename
  integer                :: vector_length
  integer                :: ncid, dimid, varid, status
  
  status = nf90_open(filename, NF90_NOWRITE, ncid)

  status = nf90_inq_varid(ncid, "snwdph", varid)
  status = nf90_get_var(ncid, varid , observation%snow_depth  , &
      start = (/1,1/), count = (/vector_length, 1/))

  status = nf90_close(ncid)

  end subroutine ReadObservation

  subroutine ReadVectorLength(filename, vector_length)
  
  use netcdf

  character*256     :: filename
  integer           :: vector_length
  integer           :: ncid, dimid, varid, status
  
  status = nf90_open(filename, NF90_NOWRITE, ncid)

  status = nf90_inq_dimid(ncid, "location", dimid)
  status = nf90_inquire_dimension(ncid, dimid, len = vector_length)
  
  status = nf90_close(ncid)

  end subroutine ReadVectorLength

  subroutine ReadNamelist(restart_dir,restart_date)
  
    character*128  :: static_file = ""
    character*128  :: init_file = ""
    character*128  :: forcing_dir = ""
    character*128  :: output_dir = ""
    
    logical        :: separate_output = .false.
  
    integer        :: timestep_seconds = -999

    integer        :: restart_frequency_s = 0
    logical        :: restart_simulation = .false.
    character*19   :: restart_date
    character*128  :: restart_dir
  
    character*19   :: simulation_start = ""
    character*19   :: simulation_end = ""

    integer        :: run_days = -999
    integer        :: run_hours = -999
    integer        :: run_minutes = -999
    integer        :: run_seconds = -999
    integer        :: run_timesteps = -999
    
    integer        :: begloc = 1
    integer        :: endloc = 1
    
  
    namelist / run_setup  / static_file, init_file, forcing_dir, output_dir, timestep_seconds, &
                            simulation_start, simulation_end, run_days, run_hours, run_minutes, &
                            run_seconds, run_timesteps, separate_output, begloc, endloc, &
                            restart_dir, restart_frequency_s, restart_simulation, restart_date

    open(30, file="ufs-land.namelist", form="formatted")
     read(30, run_setup)
    close(30)

  end subroutine ReadNamelist

end program ufsLandNoahMPReplaceRestart
