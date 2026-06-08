#include "./../fppmacros"

	subroutine write_rad()
	
	use rad
        use radae, only: abstot_3d, absnxt_3d, emstot_3d
	implicit none
	character *4 rankchar
	integer irank
        integer lenstr
        external lenstr
        character(9):: restart0 = './RESTART'
        character(11) restart

        if(masterproc) print*,'Writting radiation restart file...'

        if(restart_number_w.eq.1) then
          restart=restart0//'/'
        else
          restart=restart0//'1/'
        end if


        if(restart_sep) then

          write(rankchar,'(i4)') rank

          open(56,file=trim(restart)//case(1:len_trim(case))//'_'// &
               caseid(1:len_trim(caseid))//'_'// &
               rankchar(5-lenstr(rankchar):4)//'_restart_rad.bin', &
               status='unknown',form='unformatted',BUFFEREDYES ACTION='WRITE')
          write(56) nsubdomains
	  write(56) initrad,nradsteps,tabs_rad,qc_rad,qi_rad,qs_rad,qv_rad, &
                    cld_rad,rel_rad,rei_rad,res_rad,ozone,n2o,ch4,cfc11,cfc12, &
                    lwnsxy,swnsxy,lwntxy,swntxy,lwntmxy,swntmxy,lwnscxy,swnscxy,lwntcxy,swntcxy, &
                    swdsvisxy,swdsnirxy,swdsvisdxy,swdsnirdxy,swusvisdxy,swusnirdxy, &
                    coszrsxy,lwdsxy,swdsxy,solinxy, &
	            lwdscxy,swdscxy,qrad,absnxt_3d,abstot_3d,emstot_3d 
          close(56)

        else

          write(rankchar,'(i4)') nsubdomains

          do irank=0,nsubdomains-1

             call task_barrier()

             if(irank.eq.rank) then

               if(masterproc) then

                  open(56,file=trim(restart)//case(1:len_trim(case))//'_'// &
                      caseid(1:len_trim(caseid))//'_'// &
                      rankchar(5-lenstr(rankchar):4)//'_restart_rad.bin', &
                      status='unknown',form='unformatted',BUFFEREDYES ACTION='WRITE')
                  write(56) nsubdomains

               else

                  open(56,file=trim(restart)//case(1:len_trim(case))//'_'// & 
                      caseid(1:len_trim(caseid))//'_'// &
                      rankchar(5-lenstr(rankchar):4)//'_restart_rad.bin', &
                      status='unknown',form='unformatted', position='append',BUFFEREDYES ACTION='WRITE')

               end if

	       write(56) initrad,nradsteps,tabs_rad,qc_rad,qi_rad,qs_rad,qv_rad, &
                         cld_rad,rel_rad,rei_rad,res_rad,ozone,n2o,ch4,cfc11,cfc12, &
                         lwnsxy,swnsxy,lwntxy,swntxy,lwntmxy,swntmxy,lwnscxy,swnscxy,lwntcxy,swntcxy, &
                         swdsvisxy,swdsnirxy,swdsvisdxy,swdsnirdxy,swusvisdxy,swusnirdxy, &
                         coszrsxy,lwdsxy,swdsxy,solinxy, &
	      	         lwdscxy,swdscxy,qrad,absnxt_3d,abstot_3d,emstot_3d 
               close(56)
           end if
        end do

        end if ! restart_sep

	if(masterproc) then
           print *,'Saved radiation restart file. nstep=',nstep
	endif

        call task_barrier()

        return
        end
 
 
 
 
     
	subroutine read_rad()
	
	use rad
        use radae, only: abstot_3d, absnxt_3d, emstot_3d
	implicit none
	character *4 rankchar
	integer irank,ii
        integer lenstr
        external lenstr
        character(9):: restart0 = './RESTART'
        character(11) restart

        if(masterproc) print*,'Reading radiation restart file...'

        if(restart_number_r.eq.1) then
          restart=restart0//'/'
        else
          restart=restart0//'1/'
        end if

        if(restart_sep) then

          write(rankchar,'(i4)') rank

          if(nrestart.ne.2) then
            open(56,file=trim(restart)//trim(case)//'_'//trim(caseid)//'_'// &
              rankchar(5-lenstr(rankchar):4)//'_restart_rad.bin', &
              status='unknown',form='unformatted',BUFFEREDYES ACTION='READ')
          else
            open(56,file=trim(restart)//trim(case_restart)//'_'//trim(caseid_restart)//'_'// &
              rankchar(5-lenstr(rankchar):4)//'_restart_rad.bin', &
              status='unknown',form='unformatted',BUFFEREDYES ACTION='READ')
          end if
          read (56)
	  read(56) initrad,nradsteps,tabs_rad,qc_rad,qi_rad,qs_rad,qv_rad, &
                         cld_rad,rel_rad,rei_rad,res_rad,ozone,n2o,ch4,cfc11,cfc12, &
                         lwnsxy,swnsxy,lwntxy,swntxy,lwntmxy,swntmxy,lwnscxy,swnscxy,lwntcxy,swntcxy, &
                         swdsvisxy,swdsnirxy,swdsvisdxy,swdsnirdxy,swusvisdxy,swusnirdxy, &
                         coszrsxy,lwdsxy,swdsxy,solinxy, &
	                 lwdscxy,swdscxy,qrad,absnxt_3d,abstot_3d,emstot_3d 
          close(56)
        else

          write(rankchar,'(i4)') nsubdomains

          if(nrestart.ne.2) then
            open(56,file=trim(restart)//trim(case)//'_'//trim(caseid)//'_'// &
              rankchar(5-lenstr(rankchar):4)//'_restart_rad.bin', &
              status='unknown',form='unformatted',BUFFEREDYES ACTION='READ')
          else
            open(56,file=trim(restart)//trim(case_restart)//'_'//trim(caseid_restart)//'_'// &
              rankchar(5-lenstr(rankchar):4)//'_restart_rad.bin', &
              status='unknown',form='unformatted',BUFFEREDYES ACTION='READ')
          end if

          do irank=0,nsubdomains-1

             call task_barrier()

             if(irank.eq.rank) then

               read (56)

               do ii=0,irank-1 ! skip records
                 read(56)
               end do

	       read(56) initrad,nradsteps,tabs_rad,qc_rad,qi_rad,qs_rad,qv_rad, &
                         cld_rad,rel_rad,rei_rad,res_rad,ozone,n2o,ch4,cfc11,cfc12, &
                         lwnsxy,swnsxy,lwntxy,swntxy,lwntmxy,swntmxy,lwnscxy,swnscxy,lwntcxy,swntcxy, &
                         swdsvisxy,swdsnirxy,swdsvisdxy,swdsnirdxy,swusvisdxy,swusnirdxy, &
                         coszrsxy,lwdsxy,swdsxy,solinxy, &
	  	         lwdscxy,swdscxy,qrad,absnxt_3d,abstot_3d,emstot_3d 
               close(56)
             end if

          end do

        end if ! restart_sep

        call task_barrier()

        if(masterproc) print*,'Done with reading radiation restart file...'

        return
        end
