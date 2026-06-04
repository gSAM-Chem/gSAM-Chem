#include "./../fppmacros"
!==================================================================================
! Based on restart.f90 of SAM
! Contians write and read restart file subroutines for a certain SLM variables
!==================================================================================
subroutine write_statement_slm()
! Write slm variables in slm restart file
        use slm_vars
        use rad, only: coszrsxy, swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, lwdsxy, swdsxy
        use grid, only: restart_number_w

        implicit none
        character *4 rankchar
        character *256 filename
        integer irank
        integer, external :: lenstr
        character(9):: restart0 = './RESTART'
        character(11) restart

        if(restart_number_w.eq.1) then
          restart=restart0//'/'
        else
          restart=restart0//'1/'
        end if

        if(masterproc) print*,'Writing SLM restart file...'

        if(restart_sep) then

          write(rankchar,'(i4)') rank

          filename=trim(restart)//case(1:len_trim(case))//'_'// &
               caseid(1:len_trim(caseid))//'_'// &
               rankchar(5-lenstr(rankchar):4)//'_restart_slm.bin'
          open(176,file=trim(filename),status='unknown',form='unformatted',BUFFEREDYES ACTION='WRITE')
          write(176) nsubdomains
	  write(176) soilt, soilw, mw, mws, t_canop, t_cas, q_cas, soilt_obs, soilw_obs,&
		dosoilwnudging, dosoiltnudging,tausoil,landtype0,LAI0,clay0,sand0, &
                readland, landtypefile, LAIfile, readsoil, soilfile, readseaice, seaicefile, &
		ustar, tstar, snowt, snow_mass, nsoil, icemask, seaicemask, LAI, lhf_canop, shf_canop, &
                lhf_soil, shf_soil, evp_canop, evp_soil, coszrsxy, swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, &
                lwdsxy, swdsxy
          close(176)

        else  ! not restart_sep

          write(rankchar,'(i4)') nsubdomains

          do irank=0,nsubdomains-1

             call task_barrier()

             if(irank.eq.rank) then

               if(masterproc) then

                  open(176,file=trim(restart)//case(1:len_trim(case))//'_'// &
                      caseid(1:len_trim(caseid))//'_'// &
                      rankchar(5-lenstr(rankchar):4)//'_restart_slm.bin', &
                      status='unknown',form='unformatted',BUFFEREDYES ACTION='WRITE')
                  write(176) nsubdomains

               else

                  open(176,file=trim(restart)//case(1:len_trim(case))//'_'// & 
                      caseid(1:len_trim(caseid))//'_'// &
                      rankchar(5-lenstr(rankchar):4)//'_restart_slm.bin', &
                      status='unknown',form='unformatted', position='append',BUFFEREDYES ACTION='WRITE')

               end if
              write(176) soilt, soilw, mw, mws, t_canop, t_cas, q_cas, soilt_obs, soilw_obs,&
                dosoilwnudging, dosoiltnudging,tausoil,landtype0,LAI0,clay0,sand0, &
                readland, landtypefile, LAIfile, readsoil, soilfile, readseaice, seaicefile, &
                ustar, tstar, snowt, snow_mass, nsoil, icemask, seaicemask, LAI, lhf_canop, shf_canop, &
                lhf_soil, shf_soil, evp_canop, evp_soil, coszrsxy, swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, &
                lwdsxy, swdsxy
               close(176)
             end if
          end do

        end if ! restart_sep

        if(masterproc) then
           print *,'Saved SLM restart file. nstep=',nstep
        endif

        call task_barrier()

        return
        end
 
 
 
 
     
subroutine read_statement_slm()

use slm_vars
use rad, only: coszrsxy, swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, lwdsxy, swdsxy
use grid, only: restart_number_r

implicit none
        character *4 rankchar
        character *256 filename
        character(sizeof(landtypefile)) landtypefile1
        character(sizeof(seaicefile)) seaicefile1
        character(sizeof(LAIfile)) LAIfile1
        character(sizeof(soilfile)) soilfile1
        integer irank,ii,nsoil1
        integer it,jt,i,j,k
        integer, external :: lenstr
        character(9):: restart0 = './RESTART'
        character(11) restart

        if(restart_number_r.eq.1) then
          restart=restart0//'/'
        else
          restart=restart0//'1/'
        end if


        if(masterproc) print*,'Reading SLM restart file...'

        if(restart_sep) then

          write(rankchar,'(i4)') rank

          if(nrestart.ne.2) then
                filename =trim(restart)//case(1:len_trim(case))//'_'// &
                        caseid(1:len_trim(caseid))//'_'// &
                         rankchar(5-lenstr(rankchar):4)//'_restart_slm.bin'
          else
                filename =trim(restart)//case_restart(1:len_trim(case_restart))//'_'// &
                        caseid_restart(1:len_trim(caseid_restart))//'_'// &
                         rankchar(5-lenstr(rankchar):4)//'_restart_slm.bin'
           end if

          open(176,file=filename, status='unknown',form='unformatted',BUFFEREDYES ACTION='READ')
          read (176)
          read (176) soilt, soilw, mw, mws, t_canop, t_cas, q_cas, soilt_obs, soilw_obs,&
                dosoilwnudging, dosoiltnudging,tausoil,landtype0,LAI0,clay0,sand0, &
                readland, landtypefile1, LAIfile1, readsoil, soilfile1, readseaice, seaicefile1, &
                ustar, tstar, snowt, snow_mass, nsoil1, icemask, seaicemask, LAI, lhf_canop, shf_canop, &
                lhf_soil, shf_soil, evp_canop, evp_soil, coszrsxy, swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, &
                lwdsxy, swdsxy
          close(176)

        else

          write(rankchar,'(i4)') nsubdomains
          if(nrestart.ne.2) then
                filename =trim(restart)//case(1:len_trim(case))//'_'// &
                        caseid(1:len_trim(caseid))//'_'// &
                         rankchar(5-lenstr(rankchar):4)//'_restart_slm.bin'
          else
                filename =trim(restart)//case_restart(1:len_trim(case_restart))//'_'// &
                        caseid_restart(1:len_trim(caseid_restart))//'_'// &
                         rankchar(5-lenstr(rankchar):4)//'_restart_slm.bin'
           end if

          open(176,file=filename,status='unknown',form='unformatted',BUFFEREDYES ACTION='READ')

          do irank=0,nsubdomains-1

             call task_barrier()

             if(irank.eq.rank) then

               read (176)

               do ii=0,irank-1 ! skip records
                  read (176)
               end do
               read (176) soilt, soilw, mw, mws, t_canop, t_cas, q_cas, soilt_obs, soilw_obs,&
                dosoilwnudging, dosoiltnudging,tausoil,landtype0,LAI0,clay0,sand0, &
                readland, landtypefile1, LAIfile1, readsoil, soilfile1, readseaice, seaicefile1, &
                ustar, tstar, snowt, snow_mass, nsoil1, icemask, seaicemask, LAI, lhf_canop, shf_canop, &
                lhf_soil, shf_soil, evp_canop, evp_soil, coszrsxy, swdsvisxy, swdsnirxy, swdsvisdxy, swdsnirdxy, &
                lwdsxy, swdsxy
               close(176)
             end if

          end do

        end if ! restart_sep

        if(nsoil1.ne.nsoil) then
           if(masterproc)print*,'number of soil levels is wrong for the restart.'
           if(masterproc)print*,'Required nsoil:',nsoil,' Quitting...'
           call task_abort()
         end if
        if(nrestart.eq.1) then
         if(readland.and.landtypefile.ne.landtypefile1) then
           if(masterproc)print*,'specified landtypefile name is wrong for the restart. '
           if(masterproc)print*,'Required file:',landtypefile,' Quitting...'
           call task_abort()
         end if
         if(readland.and.LAIfile.ne.LAIfile1) then
           if(masterproc)print*,'specified LAIfile name is wrong for the restart. '
           if(masterproc)print*,'Required file:',LAIfile,' Quitting...'
           call task_abort()
         end if
         if(readsoil.and.soilfile.ne.soilfile1) then
           if(masterproc)print*,'specified soilfile name is wrong for the restart.'
           if(masterproc)print*,'Required file:',soilfile,' Quitting...'
           call task_abort()
         end if
        end if
        if(rank.eq.nsubdomains-1) then
             if(masterproc)print *,'Case:',caseid
             if(masterproc)print *,'Restart SLM at step:',nstep
             if(masterproc)print *,'Time:',nstep*dt
        endif


        call task_barrier()

        if(masterproc) print*,'done with reading SLM restart file...'

        end
