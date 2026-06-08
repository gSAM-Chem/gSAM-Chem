#include "fppmacros"

        subroutine write_all()

        use vars
        implicit none
        character *4 rankchar
        character *256 filename
        integer irank
        integer, external :: lenstr
        character(9):: restart0 = './RESTART'
        character(11) restart
        
        call t_startf ('restart_out')

        if(restart_number_w.eq.1) then
          restart=restart0//'/' 
        else
          restart=restart0//'1/' 
        end if


        if(masterproc) then
         print*,'Writing restart file ... '
         filename = trim(restart)//trim(case)//'_'//trim(caseid)//'_misc_restart.bin'
         open(66,file=trim(filename), status='unknown',form='unformatted', &
                 BUFFEREDYES ACTION='WRITE')
        end if


        if(restart_sep) then

          write(rankchar,'(i4)') rank

          filename = trim(restart)//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart.bin'


          open(65,file=trim(filename), status='unknown',form='unformatted', &
                  BUFFEREDYES ACTION='WRITE')
          write(65) nsubdomains, nsubdomains_x, nsubdomains_y

          call write_statement()


        else
          write(rankchar,'(i4)') nsubdomains
          filename = trim(restart)//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart.bin'

          do irank=0,nsubdomains-1

             call task_barrier()

             if(irank.eq.rank) then

               if(masterproc) then
      
                open(65,file=trim(filename), status='unknown',form='unformatted', &
                        BUFFEREDYES ACTION='WRITE')
                write(65) nsubdomains, nsubdomains_x, nsubdomains_y

               else

                open(65,file=trim(filename), status='unknown',form='unformatted',&
                   position='append', BUFFEREDYES ACTION='WRITE')

               end if

               call write_statement()

             end if
          end do

        end if ! restart_sep

        call task_barrier()

        call t_stopf ('restart_out')

        return
        end
 
 
 
 
     
        subroutine read_all()

        use vars
        implicit none
        character *4 rankchar
        character *256 filename
        integer irank, ii, nstep0, k
        integer, external :: lenstr
        character(9):: restart0 = './RESTART'
        character(11) restart



        if(nrestart.ne.2) then
          filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_restart_nstep.txt'
        else
          filename = './RESTART/'//trim(case_restart)//'_'//trim(caseid_restart)//'_restart_nstep.txt'
        end if
        if(masterproc) print*,'filename:',trim(filename)
        call task_barrier()
        open(67,file=trim(filename), status='unknown',form='formatted')
        read(67,*) nstep0
        read(67,*) restart_number_r
        close(67)
        if(masterproc) call systemf('cp '//trim(filename)//' '//trim(filename)//'_save')

        if(masterproc) print*,'Reading restart file ... restart_number_r=',restart_number_r
        if(dosaferestart) then
          if(restart_number_r.eq.1) then
            restart_number_w = 2
          else
            restart_number_w = 1
          end if
        end if
        if(restart_number_r.eq.1) then
          restart=restart0//'/'
        else
          restart=restart0//'1/'
        end if
      

        if(nrestart.ne.2) then
          filename = trim(restart)//trim(case)//'_'//trim(caseid)//'_misc_restart.bin'
        else
          filename = trim(restart)//trim(case_restart)//'_'//trim(caseid_restart)//'_misc_restart.bin'
        end if
        open(66,file=trim(filename), status='old',form='unformatted', BUFFEREDYES ACTION='READ')


        if(restart_sep) then

           write(rankchar,'(i4)') rank

           if(nrestart.ne.2) then
             filename = trim(restart)//trim(case)//'_'//trim(caseid)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart.bin'
           else
             filename = trim(restart)//trim(case_restart)//'_'//trim(caseid_restart)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart.bin'
           end if

           open(65,file=trim(filename), status='old',form='unformatted', BUFFEREDYES ACTION='READ')
           read(65)

           call read_statement()

        else

          write(rankchar,'(i4)') nsubdomains

          if(nrestart.ne.2) then
            filename=trim(restart)//trim(case)//'_'//trim(caseid)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart.bin'
          else
            filename=trim(restart)//trim(case_restart)//'_'//trim(caseid_restart)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart.bin'
          end if
          open(65,file=trim(filename), status='old',form='unformatted', BUFFEREDYES ACTION='READ')

          do irank=0,nsubdomains-1

             call task_barrier()

             if(irank.eq.rank) then

               read (65)
 
               do ii=0,irank-1 ! skip records
                 read(65)
                 read(65) k
                 if(k.eq.1) read(65)
                 read(65) k
                 if(k.eq.1) read(65)
               end do

               call read_statement()

             end if

          end do

        end if ! restart_sep

        call task_barrier()

        dtfactor = -1.

       ! update the boundaries 
       ! (just in case when some parameterization initializes and needs boundary points)

        call boundaries(1)
        call boundaries(4)

        return
        end
 

        subroutine write_statement()

        use vars
        use microphysics, only: micro_field, nmicro_fields
        use chemistry, only: gchem_field, ngchem_fields, aqchem_field, aqchem_gasprod_field, archem_field
        use chem_aqueous, only: naqchem_fields
        use chem_aerosol, only: narchem_fields

        use sgs, only: sgs_field, nsgs_fields, sgs_field_diag, nsgs_fields_diag, tk_factor, &
                       dochamber, t_top, t_bot, t_wall
        use tracers
        use params
        use movies, only: irecc
        use terrain
        use cup, only: dt_cu3D, dq_cu3D, du_cu3D, dv_cu3D, dcld_ls3D, dqn_ls3D, wd_cu, dwd_cu
        use buildings, only: doshadows, doirgroundfromwalls, &
                             doirwallsfromground, doirwallsfromwalls
        implicit none
        integer k

        write(65)  rank, &
           u, v, w, t, p, pp, qv, qcl, qci, qpl, qpi, cld, eis, dudt, dvdt, dwdt, u_w, u_e, v_s, v_n, & 
           tracer, micro_field, gchem_field, aqchem_field, aqchem_gasprod_field, archem_field, sgs_field, sgs_field_diag, sstxy, precinst, precac_xy, tr_ac_xy, &
           lwntac_xy, swntac_xy, lwnsac_xy, swnsac_xy, lwdsac_xy, swdsac_xy, &
           shf_ocn, lhf_ocn, taux_ocn, tauy_ocn, prec_ocn, lw_ocn, sw_ocn, u_ocn, v_ocn
        if(allocated(pfy)) then
           k=1
           write(65) k
           write(65) pfy
        else
           k=0
           write(65) k
        end if
        if(allocated(dt_cu3D)) then
           k=1
           write(65) k
           write(65) dt_cu3D, dq_cu3D, du_cu3D, dv_cu3D, dcld_ls3D, dqn_ls3D, wd_cu, dwd_cu
        else
           k=0
           write(65) k
        end if
        if(allocated(tend_in_u)) then
           k=1         
           write(65) k
           write(65) tend_in_u, tend_in_v
        else
           k=0
           write(65) k       
        end if
        close(65)
        if(masterproc) then
           write(66) dt,dx,dy
           write(66) version
           write(66) nx, ny, nz
           write(66) z, zi, adz, adzw, dz, pres, prespot, presi, prespoti, rho, rhow, bet 
           write(66) irecc, at, bt, ct, dtn, dt3, time, doconstdz, &
            day, day0, nstep, na, nb, nc, caseid, case, icycle, ncycle, &
            dodamping, doupperbound, docloud, doprecip, doradhomo, dosfchomo,&
            dolongwave, doshortwave, dosgs, dosubsidence, dotracers,  dosmoke, &
            docoriolis, dosurface, dolargescale,doradforcing, dossthomo, &
            dosfcforcing, doradsimple, donudging_uv, donudging_tq, donudging_w, &
            dowallx, dowally, doperpetual, doseasons, readlatlon, readlat, latlonfile, latlonfilebin, &
            docup, n_cup, docolumn, readterr, readlandmask, dodynamicocean, doslabocean, &
            ocean_type, delta_sst, depth_slab_ocean, Szero, deltaS, timesimpleocean, &
            pres0, ug, vg, fcor, fcorz, tabs_s, z0, n2ox, ch4x, cfc11x, cfc12x, &
            fluxt0, fluxq0, tau0, tauls, tautqls, timelargescale, epsv, ncallocean, &
            nudging_uv_z1, nudging_uv_z2, donudging_t, donudging_q, doisccp, domodis, domisr,  &
            les_s, dosimfilesout, dosolarconstant, solar_constant, zenith_angle, notracegases, &
            doSAMconditionals, dosatupdnconditionals, doirwallsfromwalls, &
            nudging_t_z1, nudging_t_z2, nudging_q_z1, nudging_q_z2, &
            ocean, land, island, sfc_flx_fxd, sfc_tau_fxd, readsst, dohs94, nub, tk_factor, &
            doirwallsfromground, doirgroundfromwalls, doshadows, dossthomozonal, &
            doseaice, seaicethickness, doradhomozonal, doseaiceevol, dotc, read_meters, &
            doimplicitdiff, dochamber, t_top, t_wall, t_bot, &
            nrad, nrad_ems, nxco2, latitude0, longitude0, dofplane, earth_factor, dodamping_w, &
            docoriolisz, doradlon, doradlat, doseawater, salt_factor, SLM, doequinox, &
            donudge3D, nudge3D_dir, nudge3D_file, nudge3D_tau, docap_snd_cu, cap_snd_cu, dobuildings, &
            doterrain, dobufferzonex, bufferzonex, dobufferzoney, bufferzoney, donobuoyancy, dometric,&
            doglobal, donearglobal, gmg_precision, dodamping_u, docyclic, cycle_period, &
            dodatefilename, nsaveM, nsaveMstart, nsaveMend, saveMbin, saveMsep, saveMavg, dogzipM, &
            sst_climo, tau_ocean, sst_mean, dosstclimo, dolatlon, reado3, doradavg, alpha_hybrid, &
            readlandmask, landmaskfile, readsst, sstfile, readterr, terrainfile, doequilocean, &
            dofixdynamics,fixdynamics_type, donodynamics, gamma_RAVE, nadv_mom, &
            ntracers, nmicro_fields, ngchem_fields, naqchem_fields, narchem_fields, nsgs_fields, nsgs_fields_diag
            close(66)
        end if
        if(rank.eq.nsubdomains-1) then
            print *,'Restart file was written at  nstep=',nstep
        endif
        if(dobuildings) call write_statement_buildings()
        return
        end


        subroutine read_statement()

        use vars
        use microphysics, only: micro_field, nmicro_fields
        use chemistry, only: gchem_field, ngchem_fields, aqchem_field, aqchem_gasprod_field, archem_field
        use chem_aqueous, only: naqchem_fields
        use chem_aerosol, only: narchem_fields

        use sgs, only: sgs_field, nsgs_fields, sgs_field_diag, nsgs_fields_diag, tk_factor, &
                       dochamber, t_top, t_bot, t_wall
        use tracers
        use params
        use cup, only: dt_cu3D, dq_cu3D, du_cu3D, dv_cu3D, dcld_ls3D, dqn_ls3D, wd_cu, dwd_cu
        use buildings, only: doshadows, doirgroundfromwalls, &
                             doirwallsfromground, doirwallsfromwalls
        use movies, only: irecc
        implicit none
        real tmp(nx,ny)

        integer  nx1, ny1, nz1, rank1, ntr, nmic, ngchem, naqchem, narchem, nsgs, nsgsd, k
        character(sizeof(case)) case1
        character(sizeof(caseid)) caseid1
        character(sizeof(version)) version1
        character(sizeof(landmaskfile)) landmaskfile1
        character(sizeof(sstfile)) sstfile1
        character(sizeof(terrainfile)) terrainfile1

        read(65)  rank1, &
           u, v, w, t, p, pp, qv, qcl, qci, qpl, qpi, cld, eis, dudt, dvdt, dwdt, u_w, u_e, v_s, v_n, &
           tracer, micro_field, gchem_field, aqchem_field, aqchem_gasprod_field, archem_field, sgs_field, sgs_field_diag, sstxy, precinst, precac_xy, tr_ac_xy, & 
           lwntac_xy, swntac_xy, lwnsac_xy, swnsac_xy, lwdsac_xy, swdsac_xy, &
           shf_ocn, lhf_ocn, taux_ocn, tauy_ocn, prec_ocn, lw_ocn, sw_ocn, u_ocn, v_ocn
        read(65) k
        if(k.eq.1) then
            allocate (pfy(nzm,ny_gl,nx_gl/nsubdomains))
            read(65) pfy
        end if
        read(65) k
        if(k.eq.1) then
            read(65) dt_cu3D, dq_cu3D, du_cu3D, dv_cu3D, dcld_ls3D, dqn_ls3D, wd_cu, dwd_cu
        end if
        read(65) k
        if(k.eq.1) then
            allocate (tend_in_u(nx,ny,nzm),tend_in_v(nx,ny,nzm))
            read(65) tend_in_u, tend_in_v
        end if
        close(65)
        read(66) dt,dx,dy
        read(66) version1
        if(version1.ne.version) then
          if(masterproc)print *,'Wrong restart file!'
          if(masterproc)print *,'Version of SAM that wrote the restart files:',version1
          if(masterproc)print *,'Current version of SAM',version
          call task_abort()
        end if
        read(66) nx1, ny1, nz1
        if(nz.eq.nz1) then
            read(66) z, zi, adz, adzw, dz, pres, prespot, presi, prespoti, rho, rhow, bet 
        else
            read(66)
        end if
        read(66) irecc, at, bt, ct, dtn, dt3, time, doconstdz, &
            day, day0, nstep, na, nb, nc, caseid1, case1, icycle, ncycle, &
            dodamping, doupperbound, docloud, doprecip, doradhomo, dosfchomo,&
            dolongwave, doshortwave, dosgs, dosubsidence, dotracers,  dosmoke, &
            docoriolis, dosurface, dolargescale,doradforcing, dossthomo, &
            dosfcforcing, doradsimple, donudging_uv, donudging_tq, donudging_w, &
            dowallx, dowally, doperpetual, doseasons, readlatlon, readlat, latlonfile, latlonfilebin, &
            docup, n_cup, docolumn, readterr, readlandmask, dodynamicocean, doslabocean, &
            ocean_type, delta_sst, depth_slab_ocean, Szero, deltaS, timesimpleocean, &
            pres0, ug, vg, fcor, fcorz, tabs_s, z0, n2ox, ch4x, cfc11x, cfc12x, &
            fluxt0, fluxq0, tau0, tauls, tautqls, timelargescale, epsv, ncallocean, &
            nudging_uv_z1, nudging_uv_z2, donudging_t, donudging_q, doisccp, domodis, domisr, &
            les_s, dosimfilesout, dosolarconstant, solar_constant, zenith_angle, notracegases, &
            doSAMconditionals, dosatupdnconditionals, doirwallsfromwalls, &
            nudging_t_z1, nudging_t_z2, nudging_q_z1, nudging_q_z2, &
            ocean, land, island, sfc_flx_fxd, sfc_tau_fxd, readsst, dohs94, nub, tk_factor, &
            doirwallsfromground, doirgroundfromwalls, doshadows, dossthomozonal,  &
            doseaice, seaicethickness, doradhomozonal, doseaiceevol, dotc, read_meters, &
            doimplicitdiff, dochamber, t_top, t_wall, t_bot, &
            nrad, nrad_ems, nxco2, latitude0, longitude0, dofplane, earth_factor, dodamping_w, &
            docoriolisz, doradlon, doradlat, doseawater, salt_factor, SLM, doequinox, &
            donudge3D, nudge3D_dir, nudge3D_file, nudge3D_tau, docap_snd_cu, cap_snd_cu, dobuildings, &
            doterrain, dobufferzonex, bufferzonex, dobufferzoney, bufferzoney, donobuoyancy, dometric, &
            doglobal, donearglobal, gmg_precision, dodamping_u, docyclic, cycle_period, &
            dodatefilename, nsaveM, nsaveMstart, nsaveMend, saveMbin, saveMsep, saveMavg, dogzipM, &
            sst_climo, tau_ocean, sst_mean, dosstclimo, dolatlon, reado3, doradavg, alpha_hybrid, &
            readlandmask, landmaskfile1, readsst, sstfile1, readterr, terrainfile1, doequilocean, &
            dofixdynamics,fixdynamics_type, donodynamics, gamma_RAVE, nadv_mom, &
            ntr, nmic, ngchem, naqchem, narchem, nsgs, nsgsd
        close(66)

        if(nstep.ge.nstop) then
           if(masterproc)print*,'Attempt to restart at  nstep greater or equal nstop!' 
           if(masterproc)print*,'nstep=',nstep,'nstop=',nstop,'Exit...'
           call task_abort()
        end if 
        if(nrestart.eq.1) then
         if(case.ne.case1) then
           if(masterproc)print*,'Case is wrong for the restart. Quitting...' 
           if(masterproc)print*,'Required case:',case 
           call task_abort()
         end if 
         if(caseid.ne.caseid1) then
           if(masterproc)print*,'Caseid is wrong for the restart. Quitting...' 
           if(masterproc)print*,'Required caseid:',caseid 
           call task_abort()
         end if 
         if(readlandmask.and.landmaskfile.ne.landmaskfile1) then
           if(masterproc)print*,'specified landmaskfile name is wrong for the restart. Quitting...' 
           if(masterproc)print*,'Required file:',sstfile1,' Quitting...' 
           call task_abort()
         end if 
         if(readterr.and.terrainfile.ne.terrainfile1) then
           if(masterproc)print*,'specified terrainfile name is wrong for the restart. ' 
           if(masterproc)print*,'Required file:',terrainfile1,' Quitting...' 
           call task_abort()
         end if 
!         if(readsst.and.sstfile.ne.sstfile1) then
!           if(masterproc)print*,'specified sstfile name is wrong for the restart. Quitting...' 
!           if(masterproc)print*,'Required file:',sstfile1,' Quitting...' 
!           call task_abort()
!         end if 
        end if
        if(nrestart.ne.2) then
          if(rank.ne.rank1) then
             if(masterproc)print *,'Error: rank of restart data is not the same as rank of the process'
             if(masterproc)print *,'rank1=',rank1,'   rank=',rank
          endif
          if(nx.ne.nx1.or.ny.ne.ny1.or.nz.ne.nz1) then
             if(masterproc)print *,'Error: domain dims (nx,ny,nz) set by grid.f'
             if(masterproc)print *,' not correspond to ones in the restart file.'
             if(masterproc)print *,'in executable:   nx, ny, nz:',nx,ny,nz
             if(masterproc)print *,'in restart file: nx, ny, nz:',nx1,ny1,nz1
             if(masterproc)print *,'Exiting...'
             call task_abort()
          endif
        end if
        if(nmic.ne.nmicro_fields) then
           if(masterproc)print*,'Error: number of micro_field in restart file is not the same as nmicro_fields'
           if(masterproc)print*,'nmicro_fields=',nmicro_fields,'   in reastart file:',nmic
           if(masterproc)print*,'Exiting...'
           call task_abort()
        end if
        if(ngchem.ne.ngchem_fields) then
           print*,'Error: number of gchem_field in restart file is not the same as ngchem_fields'
           print*,'ngchem_fields=',ngchem_fields,'   in file=',ngchem
           print*,'Exiting...'
        end if
        if(naqchem.ne.naqchem_fields) then
           print*,'Error: number of aqchem_field in restart file is not the same as nqachem_fields'
           print*,'naqchem_fields=',naqchem_fields,'   in file=',naqchem
           print*,'Exiting...'
        end if
        if(narchem.ne.narchem_fields) then
           print*,'Error: number of archem_field in restart file is not the same as narchem_fields'
           print*,'narchem_fields=',narchem_fields,'   in file=',narchem
           print*,'Exiting...'
        end if
        if(nsgs.ne.nsgs_fields.or.nsgsd.ne.nsgs_fields_diag) then
           if(masterproc)print*,'Error: number of sgs_field in restart file is not the same as nsgs_fields'
           if(masterproc)print*,'nsgs_fields=',nsgs_fields,'   in reastart file=',nsgs
           if(masterproc)print*,'nsgs_fields_diag=',nsgs_fields_diag,'   in file=',nsgsd
           if(masterproc)print*,'Exiting...'
             call task_abort()
        end if
        if(ntr.ne.ntracers) then
           if(masterproc)print*,'Error: number of tracers in restart file is not the same as ntracers.'
           if(masterproc)print*,'ntracers=',ntracers,'   ntracers(in file)=',ntr
           if(masterproc)print*,'Exiting...'
             call task_abort()
        end if
        close(65)
        if(rank.eq.nsubdomains-1) then
           print *,'Case:',caseid
           print *,'Restarting at step:',nstep
           print *,'Time(s):',nstep*dt
           print *,'day:',day
        endif
        if(masterproc) print*,'Done...'

        return
        end






