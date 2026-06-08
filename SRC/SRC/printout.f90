! print out various parameters from the namelist, etc

	subroutine printout()

	use vars
        use tracers
	use params
        use buildings
	implicit none
        integer k
	
        select case (nrestart) 
	case(0) 
	  print*,'New Run.'
	case(1)
	  print*,'Continuation Run'
	case(2) 
          print*,'Branch run Restart Case:',trim(case_restart)
          print*,'Case-ID:',trim(caseid_restart)
        end select
	
	if(dt.gt.0) then
	   print*, 'Timestep (sec):',dt
	else
	   print*, 'Error: dt =',dt
	   call task_abort()
	endif   	

	print*,'day:',day

        print*,'doglobal = ',doglobal
        print*,'doglobalpresents = ',doglobalpresets
        print*,'donearglobal = ',donearglobal
        print*,'earth_factor = ',earth_factor
        print*,'gamma_RAVE = ',gamma_RAVE
        print*,'doflat = ',doflat
        print*,'dofcast =',dofcast
        print*
        print*,'doregion=',doregion
        print*,'dowallx=',dowallx
        print*,'dowally=',dowally
        print*,'bufferzonex=',bufferzonex
        print*,'bufferzoney=',bufferzoney
        print*,'doterrain=',doterrain
	print*
        print*,'docyclic=',docyclic
        if(docyclic) print*,'cycle_period = ',cycle_period 
        print*

        print*,'nadv_mom=',nadv_mom
	if(nadams.eq.2.or.nadams.eq.3) then
	   print*, 'Adams-Bashforth scheme order: nadams=',nadams
	else
	   print*, 'Error: nadams =',nadams
	   call task_abort()
	endif 
	  	
        print*,'dolatlon = ', dolatlon
        if(dolatlon) print*,'gmg_precision=',gmg_precision
        print*,'dofliplon = ',dofliplon
        if(dolatlon) then
           print*,' lon: min/max=',minval(lon_gl(1:nx_gl)),maxval(lon_gl(1:nx_gl))
           print*,' lat: min/max=',minval(lat_gl(1:ny_gl)),maxval(lat_gl(1:ny_gl))
           print*,' mu: cos(lat) min/max=',minval(mu_gl(1:ny_gl)),maxval(mu_gl(1:ny_gl))
           print*,' muv: cos(latv)min/max =',minval(muv_gl(1:ny_gl+YES3D)),maxval(muv_gl(1:ny_gl+YES3D))
        end if
	print*
	if(dx.gt.0.and.dy.gt.0.and.dz.gt.0. ) then
	   print*, 'Global Grid:',nx_gl,ny_gl,nz_gl
	   print*, 'Local Grid:',nx,ny,nzm
           print*,'doyvar=',doyvar
	   print*, 'Base grid spacing (km) dx (at equator for global grid):', 0.001*dx
	   print*, 'Base grid spacing (km) dy (at equator for global grid):', 0.001*dy*YES3D
	   print*, 'Grid spacing (km) in x at domain center:', 0.001*dx*mu_gl(ny_gl/(1+YES3D)) 
	   print*, 'Grid spacing (km) in y at domain center:', 0.001*(yv_gl(ny_gl/(1+YES3D)+YES3D)-yv_gl(ny_gl/(1+YES3D)))
	   print*, 'Min Grid spacing (km) in x (at pole for global grid):', 0.001*dx*minval(mu_gl(1:ny_gl)) 
	   print*, 'Max Grid spacing (km) in x (at equator for global grid):', 0.001*dx*maxval(mu_gl(1:ny_gl)) 
	   print*, 'Min Grid spacing (km) in y (at equator for global grid):', &
                       0.001*minval(yv_gl(1+YES3D:ny_gl+YES3D)-yv_gl(1:ny_gl)) 
	   print*, 'Max Grid spacing (km) in y (at pole for global grid):', &
                       0.001*maxval(yv_gl(1+YES3D:ny_gl+YES3D)-yv_gl(1:ny_gl)) 
           if(dolatlon) then
             print*,'grid spacing in lon (degrees):', (lonu_gl(nx_gl+1)-lonu_gl(1))/float(nx_gl)
             print*,'Min grid spacing in lat (degrees):',minval(latv_gl(1+YES3D:ny_gl+YES3D)-latv_gl(1:ny_gl))
             print*,'Max grid spacing in lat (degrees):',maxval(latv_gl(1+YES3D:ny_gl+YES3D)-latv_gl(1:ny_gl))
           end if
	   print*, 'Domain dimensions (at center, km):', 0.001*dx*nx_gl*mu_gl(ny_gl/(1+YES3D)), " x   ",&
                                                   0.001*(yv_gl(ny_gl+YES3D)-yv_gl(1))
           print*,'Domain top height (km):',0.001*z(nzm)
           print*,'Layer centers: (z, pres):'
           do k=1,nzm
            print*,k,z(k),pres(k)
           end do
           print*
           print*,'Layer interfaces: (zi, presi):'
           do k=1,nz
            print*,k,zi(k),presi(k)
           end do
	else
	   print*, 'Error: grid spacings dx, dy, dz:',dx, dy, dz
	   call task_abort()
	endif
	if(.not.(OCEAN.or.LAND.or.ISLAND)) then
	 print*, 'Neither OCEAN nor LAND nor ISLAND is set. Exitting...'
	 call task_abort()
	endif
	if(OCEAN.and.LAND) then
	 print*, 'Both OCEAN and LAND are set. Confused...'
	 call task_abort()
	endif
        if(OCEAN.and.ISLAND) then
         print*, 'Both OCEAN and LAND are set. Confused...'
         call task_abort()
        endif
        if(OCEAN.and.ISLAND) then
         print*, 'Both OCEAN and LAND are set. Confused...'
         call task_abort()
        endif
        print*
        print*,'npressure_iter:',npressure_iter
        print*
        print*,'dobuildings=',dobuildings
        print*,'doshadows=',doshadows
        print*,'doirgroundfromwalls=',doirgroundfromwalls
        print*,'doirwallsfromground=',doirwallsfromground
        print*
	print*, 'Finish at timestep: nstop=',nstop
	print*, 'Finish on day:',day+(nstop-nstep)*dt/3600./24.
	print*
	print*, 'Statistics file ouput frequency nstat=: ',nstat,' steps'
	print*, 'Statistics file sampling: every nstatfrq=',nstatfrq,' steps'
	print*, 'printouts frequency nprint=',nprint,' steps'
	print*, 'restart frequency nrestart_steps=',nrestart_steps,' steps'
	print*, 'nrestart_skip=',nrestart_skip,' steps'
        if(mod(nstop,nstat).ne.0.and.mod(nelapse,nstat).ne.0) then
          print*, 'Error: job will finish before statistics is done.'
          print*, 'Hint: nstop=',nstop, 'or/and nelapse=',nelapse, &
                  'should be divisible by nstat=',nstat
          call task_abort()
        endif
	print*,'cloud statistics type LES_S',LES_S
	print*,'cloud water/ice path threshold for cld. fraction>0: ',cwp_threshold
        print*,'domain translation velocities ug, vg:', ug, vg
        print*	
	print*, 'do column model: docolumn = ', docolumn
	print*, 'do convective parameterization: docup =', docup
	print*, 'frequency of calling (timesteps): n_cup', n_cup
        print*
	print*, 'dodamping = : ', dodamping
	print*, 'dodamping_poles = ', dodamping_poles
	print*, 'doupperbound = ',doupperbound
	print*, 'docloud = ',docloud
	print*, 'doprecip = ',doprecip
        if(docloud.and.dosmoke) then
          if(masterproc) print*,'docloud and dosmoke can not be true simultaneously'
	  call task_abort()
        end if
	 print*,'dosmoke = ',dosmoke
	print*, 'dosgs = ',dosgs
	print*, 'doimplicitdiff = ',doimplicitdiff
        print*
	print*, 'coriolis force is allowed: docoriolis = ',docoriolis	
	print*, 'vertical coriolis force is allowed: docoriolisz',docoriolisz
	if(docoriolis) then
            print*, 'do f-plane approximation:',dofplane
            if(dofplane) then
		print*, '   Coriolis parameter (1/s):',fcor
		print*, '   Vertical Coriolis parameter (1/s):',fcorz
            else
                print*, '   Coriolis parameter is the function of latitude'
            end if
	endif	
        if(doradforcing.and.(dolongwave.or.doshortwave)) then
          print*, 'prescribed rad. forcing and radiation '// &
          'calculations cannot be done at the same time.'
          call task_abort()
        endif
        if(dolongwave) then
            print*, 'dolongwave = ',dolongwave  
            print*,'nrad_ems = ',nrad_ems
        end if
        if(doshortwave) then
            print*, 'doshortwave = ',doshortwave
            print*, 'doseasons = ',doseasons
            print*, 'doperpetual = ',doperpetual
            print*, 'compute_reffc = :',compute_reffc
            print*, 'compute_reffi = :',compute_reffi
	endif
        if(dolongwave.or.doshortwave) print*,'nrad=',nrad
        print*,'doradavg = ',doradavg
        print*,'doradsimple = ',doradsimple
        print*,'doradhomo =',doradhomo
        print*,'doradhomozonal =',doradhomozonal
        print*,'dosolarconstant = ',dosolarconstant
        print*,'solar_constant = ',solar_constant
        print*,'zenith_angle = ',zenith_angle
        print*,'doradlon = ',doradlon
        print*,'doradlat = ',doradlat
        print*,'latitude0:',latitude0
        print*,'longitude0:',longitude0
        print*,'rundatadir = ',trim(rundatadir)
        print*,'doradforcing = ', doradforcing
        print*,'nxco2=',nxco2
	print*,'dosurface = ',dosurface
        print*,'dotc=',dotc
        print*,'lhf_fudge=',lhf_fudge
        print*,'shf_fudge=',shf_fudge
        print*,'tau_fudge=',tau_fudge
        if(dosurface) then
            if(LAND.or.ISLAND) then
               if(LAND)print*,'Surface type: LAND'
               if(ISLAND)then
                 print*,'Surface type: ISLAND'
                 if(.not.readlandmask) then
                     print*,'simple island:' 
                     if(doroundisland) then
                       print*,'doroundisland=',doroundisland
                       print*,'island_radius=',island_radius
                     else
                       print*,'island_x1=',island_x1
                       print*,'island_x2=',island_x2
                       print*,'island_y1=',island_y1
                       print*,'island_y2=',island_y2
                     end if
                 end if
               end if
               print*,'z0=',z0
            end if
            if((ISLAND).and.(SFC_FLX_FXD.or.SFC_TAU_FXD)) then
              print*, 'When ISLNAD is .true., surface fluxes cannot be prescribed. Quitting...'
              call task_abort()
            endif
            if(OCEAN) print*,'Surface type: OCEAN'
            print*, ' sensible heat flux prescribed SFC_FLX_FXD=',SFC_FLX_FXD
            if(SFC_FLX_FXD.and..not.dosfcforcing) print*, 'fluxt0 (W/m2)=',fluxt0*rhow(1)*cp
            print*, ' latent heat flux prescribed SFC_FLX_FXD=',SFC_FLX_FXD
            if(SFC_FLX_FXD.and..not.dosfcforcing) print*, 'fluxq0 (W/m2)=',fluxq0*rhow(1)*lcond
            print*, ' surface stress prescribed SFC_TAU_FXD=',SFC_TAU_FXD
            if(SFC_TAU_FXD.and..not.dosfcforcing) print*, 'tau0 (m2/s2)=',tau0
            print*,'dosfchomo = ',dosfchomo
        endif
        print*,'doisccp = ',doisccp
        print*,'domodis = ',domodis
        print*,'domisr = ',domisr
        print*,'dosimfilesout=',dosimfilesout
        print*	
        print*,'timeslargescale=',timelargescale
	print*, 'larger-scale subsidence is on:',dosubsidence
	print*, 'dolargescale = ',dolargescale
        if(dolargescale.or.dosubsidence) then
          if(    day.lt.dayls(1) &
            .or.day+(nstop-nstep)*dt/86400..gt.dayls(nlsf)) then
             print*,'Error: simulation time (from start to stop)'// &
              'can be beyond the l.s. forcing intervals'
             print*,'current day=',day
             print*,'stop day=',day+(nstop-nstep)*dt/86400.
             print*,'ls forcing: start =',dayls(1)
             print*,'ls forcing:   end =',dayls(nlsf)
             call task_abort()
          endif
        endif
        print*,'docap_snd_cu=',docap_snd_cu
        print*,'cap_snd_cu=',cap_snd_cu
        print*,'do_cap_wind=',do_cap_wind
        print*,'cap_wind=',cap_wind
        print*,'ocean_type =',ocean_type
        print*,'Initial SST tabs_s = ',tabs_s
        print*,'SST sin-amplitude: delta_sst = ',delta_sst
        print*,'SST shift in latitude: shift_sst = ',shift_sst
        print*,'mean ocean transport: Szero = ',Szero
        print*,'ocean transport linear max variation: deltaS = =',deltaS
        print*,'dodynamicocean =',dodynamicocean
        if(dodynamicocean) print*,'ncallocean =',ncallocean
        print*,'doslabocean =',doslabocean
        print*,'doequilocean =',doequilocean
        if(doslabocean.and.dodynamicocean) then
              print*, 'error: doslabocean and dodynamicocean cannot be both set true...'
              call task_abort()
        endif
        if(doslabocean.or.doequilocean) then
             print*,'depth_slab_ocean = ',depth_slab_ocean
             print*,'dossthomo =',dossthomo
             print*,'dossthomozonal =',dossthomozonal
             print*,'doseaice =',doseaice
             if(doseaice) print*,'seaicethickness =',seaicethickness
             print*,'doseaiceevol =',doseaiceevol
             print*,'dosstclimo =',dosstclimo
             print*,'sst_climo=',sst_climo
             print*,'tau_ocean (s)=',tau_ocean
        else
           if(dosstclimo.or.dossthomo) then
             print*,'dosstclimo or dossthome cannot be set to true '//&
                     'only when doslabcocean is true'
	     call task_abort()
           end if
        end if
        print*,'SLM=',SLM
        if(dosurface.and.dosfcforcing) then
          print*,'surface temperature prescribed: T'
          if(doslabocean) then
             print*,'ocean_type =',ocean_type
	     print*, 'doslabocean cannot be set to T'// &
                     'when dosfcforcing is also T'
	     call task_abort()
	  end if
          if(day.lt.daysfc(1).or.day+(nstop-nstep)*dt/86400..gt.daysfc(nsfc)) then
             print*,'Error: simulation time (from start to stop)'// &
              'can be beyond the sfc forcing intervals'
             print*,'current day=',day
             print*,'stop day=',day+(nstop-nstep)*dt/86400.
             print*,'sfc forcing:start =',daysfc(1)
             print*,'sfc forcing:  end =',daysfc(nsfc)
             call task_abort()
          end if
        endif
        if(doradforcing) then
          if ( day.lt.dayrfc(1) &
            .or.day+(nstop-nstep)/86400.*dt.gt.dayrfc(nrfc))then
             print*,'Error: simulation time (from start to stop)'// &
              'can be beyond the rad. forcing intervals'
             print*,'current day=',day
             print*,'stop day=',day+(nstop-nstep-1)*dt/86400.
             print*,'rad forcing:start =',dayrfc(1)
             print*,'rad forcing:  end =',dayrfc(nrfc)
             call task_abort()
          endif
        endif
        print*,'doseawater=',doseawater
        print*,'salt_factor =',salt_factor

        print*, 'nudge3Dstep_start:', nudge3Dstep_start
        print*, 'nudge3Dstep_end:', nudge3Dstep_end

        print*, 'donudging_uv:', donudging_uv
        print*, 'donudging_w:', donudging_w
        if(donudge3D.or.doregion) then
         print*,'nudge3D_tau = ',nudge3D_tau
        else
         print*,'tauls = ',tauls
        end if
        print*, 'donudging_tq:', donudging_tq
        print*, 'donudging_t:', donudging_t
        print*, 'donudging_q:', donudging_q
        print*,'nudging_uv_z1 = ', nudging_uv_z1
        print*,'nudging_uv_z2 = ', nudging_uv_z2
        print*,'nudging_t_z1 = ', nudging_t_z1
        print*,'nudging_t_z2 = ', nudging_t_z2
        print*,'nudging_q_z1 = ', nudging_q_z1
        print*,'nudging_q_z2 = ', nudging_q_z2

        print*,'dospectralnudging=',dospectralnudging
        if(dospectralnudging) then
         print*,'nx_spectral = ',nx_spectral
         print*,'ny_spectral = ',ny_spectral
         print*,'nstep_spectral = ',nstep_spectral
        end if

        print*,'dotracers =',dotracers
        print*,'ntracers =',ntracers
        if(.not.dotracers.and.ntracers.gt.0) then
         print*, 'ntracers should be set to 0 when dotracers=F'
         call task_abort()
        end if
        print*,'dotrsfcflux=',dotrsfcflux
        if(dotracers) then
           if(ntracers.eq.0) then
             print*,'dotracers is set to .true., yet ntracers = 0. Aborting ...'
             call task_abort()
           end if
           print*,ntracers, ' tracers are included'
           print*,'Tracer names:',tracername(1:ntracers)
        end if
        print*,'perturb_type = ',perturb_type
        print*,'doSAMconditionals =',doSAMconditionals
        print*,'dosatupdnconditionals =',dosatupdnconditionals
        print*,'doscamiopdata =',doscamiopdata
        print*,'iopfile:',trim(iopfile)
        print*,'dozero_out_day0 =',dozero_out_day0
        print*
        print*,'docolumn = ',docolumn
        print*,'doensemble = ',doensemble
        print*,'nensemble = ',nensemble
        print*
        print*,'save2Dnetcdf = ',save2Dnetcdf
        print*,'nsave2D = ',nsave2D
        print*,'nsave2Dstart = ',nsave2Dstart
        print*,'nsave2Dend = ',nsave2Dend
        print*,'save2Dbin = ',save2Dbin
        print*,'save2Dsep = ',save2Dsep
        print*,'save2Davg = ',save2Davg
        print*,'save2Drada = ',save2Drada
        print*,'save2Dradac = ',save2Dradac
        print*,'snow2Dout = ',snow2Dout
        print*,'dogzip2D = ',dogzip2D
        print*
        print*,'saveMnetcdf = ',saveMnetcdf
        print*,'nsaveM = ',nsaveM
        print*,'nsaveMstart = ',nsaveMstart
        print*,'nsaveMend = ',nsaveMend
        print*,'saveMbin = ',saveMbin
        print*,'saveMsep = ',saveMsep
        print*,'saveMavg = ',saveMavg
        print*,'dogzipM = ',dogzipM
        print*
        print*,'save2DLnetcdf = ',save2DLnetcdf
        print*,'nsave2DL = ',nsave2DL
        print*,'nsave2DLstart = ',nsave2DLstart
        print*,'nsave2DLend = ',nsave2DLend
        print*,'save2DLbin = ',save2DLbin
        print*,'save2DLsep = ',save2DLsep
        print*,'save2DLavg = ',save2DLavg
        print*,'dogzip2DL = ',dogzip2DL
        print*
        print*,'save2DZnetcdf = ',save2DZnetcdf
        print*,'nsave2DZ = ',nsave2DZ
        print*,'nsave2DZstart = ',nsave2DZstart
        print*,'nsave2DZend = ',nsave2DZend
        print*,'save2DZsep = ',save2DZsep
        print*,'dogzip2DZ = ',dogzip2DZ
        print*
        print*,'save3Dnetcdf = ',save3Dnetcdf
        print*,'nsave3D = ',nsave3D
        print*,'nsave3Dstart = ',nsave3Dstart
        print*,'nsave3Dend = ',nsave3Dend
        print*,'save3Dbin = ',save3Dbin
        print*,'save3Dsep = ',save3Dsep
        print*,'nfiles3D = ',nfiles3D
        print*,'dogzip3D = ',dogzip3D
        if(mod(nsubdomains,nfiles3D).ne.0) then
         print*,'total number of tasks should be divisible by nfiles3D. exiting...'
         call task_abort()
        end if
        print*,'qnsave3D = ',qnsave3D
        print*
        print*,'nstatmom = ',nstatmom
        print*,'nstatmomstart = ',nstatmomstart
        print*,'nstatmomend = ',nstatmomend
        print*,'savemombin = ',savemombin
        print*,'savemomsep = ',savemomsep
        print*
        print*,'dostatcoars = ',dostatcoars
        if(dostatcoars) then
          print*,'nstatcoars = ',nstatcoars
          print*,'nstatcoarsstart = ',nstatcoarsstart
          print*,'nstatcoarsend = ',nstatcoarsend
          print*,'savecoarsdbin = ',savecoarsbin
        end if
        print*
        print*,'nmovie = ',nmovie
        print*,'nmoviestart = ',nmoviestart
        print*,'nmovieend = ',nmovieend
        
	return
	end

