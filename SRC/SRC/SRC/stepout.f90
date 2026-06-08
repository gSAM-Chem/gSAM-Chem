subroutine stepout(nstatsteps)

use vars
use rad, only: qrad
use tracers
use microphysics, only: micro_print, micro_statistics, micro_field,index_water_vapor 
use chemistry, only: chem_print, chem_statistics
use sgs, only: tk, sgs_print, sgs_statistics, dodns
use stat_moments, only: compmoments
use movies, only: mvmovies
use params
use hbuffer
use instrument_diagnostics, only : isccp_write, modis_write, misr_write, zero_instr_diag
use terrain
use slm_vars, only: slm_stepout, slm_stat_2Dinit
use stat_coars, only: write_coars2D, write_coars3D
use buildings, only: buildings_stepout
use cup, only: precip_cu, precip_ls, wd_cu
implicit none	
	
integer i,j,k,ic,jc,nstatsteps
real qnmax(1), qnmax1(1)
real(8) buffer(6), buffer1(6)
real rbuf(1), rbuf1(1)
integer ibuf(1), ibuf1(1)
character(14) date
real, external :: qsatw,qsati
real(8) day1
real fld_mean, fldm(nzm)
real tmp(nx,ny,nzm)
character(256) filename
real(8), external :: utTimeSeconds
integer itime(6)

real(8) :: time_new_seconds, elapsed_time



if(dodatefilename) then
 day1 = day0 + (real(nstep,8)*dt+0.5_8)/86400._8
 call date_from_dayofyear(day1,year0,datechar)
 date_pr = datechar(1:4)//"-"//datechar(5:6)//"-"//datechar(7:8)//"-"//datechar(9:10)&
          //"-"//datechar(11:12)//"-"//datechar(13:14)
 read(datechar,'(i4,i2,i2,i2,i2,i2)') itime(:)
 if(itime(1).ge.1900) timeUTsec = &
                  utTimeSeconds(itime(1),itime(2),itime(3),itime(4),itime(5),itime(6))
 year = itime(1)
else
 date_pr = ''
end if

fld_mean = 0.

if(nstatsteps.eq.-1) goto 1000

! do k=1,nzm
!    tmp(1:nx,1:ny,k) = p(1:nx,1:ny,k,nb)*rho(k)
! end do
! call fminmax_print('p:',tmp(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)

if(mod(nstep,nstatis).eq.0) then
      call statistics()
      call micro_statistics()
      if(dochem) call chem_statistics()
      call sgs_statistics()
end if

if(mod(nstep,nmovie).eq.0.and.nstep.ge.nmoviestart &
                                   .and.nstep.le.nmovieend) then
      call mvmovies()
endif

if(mod(nstep,nstatmom).eq.0.and.nstep.ge.nstatmomstart &
                                   .and.nstep.le.nstatmomend) then
      call compmoments()
endif

if(dostatcoars.and.mod(nstep,nstatcoars).eq.0.and.nstep.ge.nstatcoarsstart &
                                   .and.nstep.le.nstatcoarsend) then
      call write_coars2D()
      call write_coars3D()
endif


if(mod(nstep,nstat).eq.0) then
  if(masterproc) print *,'Writting statistics:nstatsteps=',nstatsteps

  call t_startf ('stat_out')

  call hbuf_average(nstatsteps)
  call hbuf_write(nstatsteps)
  call hbuf_flush()  
  nstatsteps = 0
  call isccp_write()
  call modis_write()
  call misr_write()
  call zero_instr_diag()

  call t_stopf ('stat_out')

endif


if(.not.donorestart) then
 if(mod(nstep,nrestart_steps*(1+nrestart_skip)).eq.0.or.nstep.eq.nstop.or.nelapse.eq.0) then
  call write_all() ! save restart file
  if((LAND.or.ISLAND).and.SLM) call write_statement_slm()
  if(dobuildings) call write_statement_buildings()
 ! handle rafe restart by writing a file with the latest information
 ! where latest restart was weitten
  if(masterproc) then
      print*,'Done... restart_number_w =',restart_number_w
      filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_restart_nstep.txt'
      open(67,file=trim(filename), status='unknown',form='formatted')
      write(67,*) nstep
      write(67,*) restart_number_w
      close(67)
  end if

!  ! print out elapsed time since last restart file was written
  call get_system_time_real8(time_new_seconds)
  if(masterproc) write(*,999) time_new_seconds-time_old_seconds, real(nstep-oldstep)*dt/3600.
  999   format('CPU TIME = ',f12.4, ' OVER ', f8.2,' MODEL HOURS')

  if(nelapse_minutes.gt.0) then
     ! check to see if run should stop at the end of this time step
     !   since a clean set of restart files have been written.
     elapsed_time = time_new_seconds-time_init_seconds
     if ((elapsed_time + 2.*(time_new_seconds-time_old_seconds)) &
          .gt.60.*real(nelapse_minutes)) then
        nelapse=0 ! Job will stop when nelapse=0
        if(masterproc) then
           write(*,*) 'Job stopping at end of this step -- nelapse_minutes approaching...'
           call write_ready_to_restart_file('TRUE')
        end if
     end if

  end if

  time_old_seconds = time_new_seconds
  oldstep = nstep

 end if ! mod(nstep...)
end if ! .not.donorestart

call averageXY_SFC(prec_xy(:,:)/dtp*86400.,1,nx,1,ny,fld_mean)

if(mod(nstep,nsave2D).eq.0.and.nstep.ge.nsave2Dstart &
                                   .and.nstep.le.nsave2Dend) then
  call t_startf ('2D_out')
  call write_fields2D()
  call stat_2Dinit()
  call t_stopf ('2D_out')

endif

if(.not.save2Davg.or.nstep.eq.nsave2Dstart-nsave2D) call stat_2Dinit()

if(mod(nstep,nsaveM).eq.0.and.nstep.ge.nsaveMstart &
                                   .and.nstep.le.nsaveMend) then
  call t_startf ('MISC_out')
  call write_fields2DM()
  call stat_Minit()
  call t_stopf ('MISC_out')

endif

if(.not.saveMavg.or.nstep.eq.nsaveMstart-nsaveM) call stat_Minit()

if(SLM.and.mod(nstep,nsave2DL).eq.0.and.nstep.ge.nsave2DLstart &
                                   .and.nstep.le.nsave2DLend) then
  call t_startf ('2DL_out')
  call slm_write2D()
  call slm_stat_2Dinit()
  call t_stopf ('2DL_out')
endif

if(SLM.and..not.save2DLavg.or.nstep.eq.nsave2DLstart-nsave2DL) call slm_stat_2Dinit()

if(mod(nstep,nsave2DZ).eq.0.and.nstep.ge.nsave2DZstart.and.nstep.le.nsave2DZend) then
  call t_startf ('2DZ_out')
  call write_fields2DZ()
  call t_stopf ('2DZ_out')
endif


if(mod(nstep,nsave3D).eq.0.and.nstep.ge.nsave3Dstart.and.nstep.le.nsave3Dend) then
  ! determine if the maximum cloud water exceeds the threshold
  ! value to save 3D fields:
  call t_startf ('3D_out')
  qnmax(1)=0.
  do k=1,nzm
    do j=1,ny
      do i=1,nx
         qnmax(1) = max(qnmax(1),qcl(i,j,k))
         qnmax(1) = max(qnmax(1),qci(i,j,k))
      end do
    enddo
  enddo
  if(dompi) then
     call task_max_real(qnmax,qnmax1,1)
     qnmax(1) = qnmax1(1)
  end if
  if(qnmax(1).ge.qnsave3D) then 
    call write_fields3D()
  end if
  call t_stopf ('3D_out')
  call t_prf(rank)
endif

call t_startf ('stepout')

1000 continue

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Print stuff out:

if(NSTEP.ne.0.and.dompi) then
    rbuf(1) = error_max
    call task_barrier()
    call task_max_real(rbuf,rbuf1,1)
    error_max = rbuf1(1)
    ibuf(1) = niter_max
    call task_max_integer(ibuf,ibuf1,1)
    niter_max = ibuf1(1)
end if

if(nstatsteps.eq.-1) then
    if(masterproc) write(*,'(a,i7,a,i2,a,f6.3)') 'NSTEP = ',nstep,' NCYCLE= ',ncycle,' dtn=',dtn
else
    if(masterproc) write(*,'(a,i7,a,i2,a,f6.3,a,f4.2,a,f4.2,a,f4.2,a,f4.2,a,i2,a,g9.2)') &
      'NSTEP =',nstep,' NCYCLE=',ncycle,' dtn=',dtn, &
      ' CFL=',cfl_adv,' CFLH=',cfl_advh,' CFLZ=',cfl_advz,' CFLG=',cfl_sgs, ' niter=',niter_max,' err=',error_max
end if
ncycle = 1
!call fminmax_print('p:',p,0,nx,1-YES3D,ny,nzm)
if(nstatsteps.lt.0.or.mod(nstep,nprint).eq.0.and.icycle.eq.1) then

 if(nstatsteps.ne.-1.and.NSTEP.ne.0.and.dompi) then
   buffer(1) = total_water_before
   buffer(2) = total_water_after
   buffer(3) = total_water_evap
   buffer(4) = total_water_prec
   buffer(5) = total_water_ls
   buffer(6) = total_water_adv
   call task_sum_real8(buffer, buffer1,6)
   total_water_before = buffer1(1)
   total_water_after = buffer1(2)
   total_water_evap = buffer1(3)
   total_water_prec = buffer1(4)
   total_water_ls = buffer1(5)
   total_water_adv = buffer1(6)
 end if
!if(masterproc) then

!print*,'--->',tk(27,1,1)
!print*,'div->:'
!write(6,'(16f7.2)')((div(i,1,k)*1.e6,i=1,16),k=nzm,1,-1)
!print*,'tk->:'
!write(6,'(16f7.2)')((tk(i,1,k),i=1,16),k=nzm,1,-1)
!print*,'p->:'
!write(6,'(16f7.2)')((p(i,1,k),i=1,16),k=nzm,1,-1)
!print*,'u->:'
!write(6,'(16f7.2)')((u(i,1,k),i=1,16),k=nzm,1,-1)
!print*,'v->:'
!write(6,'(16f7.2)')((v(i,1,k),i=1,16),k=nzm,1,-1)
!print*,'w->:'
!write(6,'(16f7.2)')((w(i,1,k),i=1,16),k=nzm,1,-1)
!if(nstep.eq.3) stop
!print*,'qcl:'
!write(6,'(16f7.2)')((qcl(i,1,k)*1000.,i=16,31),k=nzm,1,-1)
!print*,'qpl->:'
!write(6,'(16f7.2)')((qpl(i,1,k)*1000.,i=16,31),k=nzm,1,-1)
!print*,'qcl:->'
!write(6,'(16f7.2)')((qci(i,1,k)*1000.,i=16,31),k=nzm,1,-1)
!print*,'qpl->:'
!write(6,'(16f7.2)')((qpi(i,1,k)*1000.,i=16,31),k=nzm,1,-1)
!print*,'qrad->:'
!write(6,'(16f7.2)')((qrad(i,1,k)*3600.,i=16,31),k=nzm,1,-1)
!print*,'qv->:'
!write(6,'(16f7.2)')((qv(i,1,k)*1000.,i=16,31),k=nzm,1,-1)
!print*,'t->:'
!write(6,'(16f7.2)')((t(i,1,k),i=1,16),k=nzm,1,-1)
!print*,'tabs->:'
!write(6,'(16f7.2)')((tabs(i,1,k),i=16,31),k=nzm,1,-1)

!end if

!--------------------------------------------------------
 if(masterproc) then
	
    print*,'time since start of run (secs) = ',time	
    print*,'time since start of run (days) = ',time/86400_8	
    print*,'DAY = ',day	
    if(dodatefilename) print*,trim(date_pr), '  timeUTsec=',timeUTsec
    write(6,*) 'NSTEP=',nstep
    write(6,*) 'surface pressure=',pres0

 endif
 call check_div(u(1:nx+1,1:ny,1:nzm),v(1:nx,1:ny+YES3D,1:nzm),w(1:nx,1:ny,1:nz))

 do k=1,nzm
    tmp(1:nx,1:ny,k) = p(1:nx,1:ny,k,nb)*rho(k)
 end do
 if(doterrain) then
  call fminmax_printm('u:',u(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,terrau(1:nx,1:ny,1:nzm).eq.1)
  call fminmax_printm('v:',v(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,terrav(1:nx,1:ny,1:nzm).eq.1)
  call fminmax_printm('w:',w(1:nx,1:ny,1:nz),1,nx,1,ny,nz,terraw(1:nx,1:ny,1:nz).eq.1)
  call fminmax_printm('p:',tmp(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,terra(1:nx,1:ny,1:nzm).eq.1)
  call fminmax_printm('pp:',pp(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,terra(1:nx,1:ny,1:nzm).eq.1)
  call fminmax_printm('t:',t(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,terra(1:nx,1:ny,1:nzm).eq.1)
  call fminmax_printm('tabs:',tabs(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,terra(1:nx,1:ny,1:nzm).eq.1)
  call fminmax_printm('qv:',qv(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,terra(1:nx,1:ny,1:nzm).eq.1)
  call fminmax_printm('u(mount):',u(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,terrau(1:nx,1:ny,1:nzm).eq.0)
  call fminmax_printm('v(mount):',v(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,terrav(1:nx,1:ny,1:nzm).eq.0)
  call fminmax_printm('w(mount):',w(1:nx,1:ny,1:nz),1,nx,1,ny,nz,terraw(1:nx,1:ny,1:nz).eq.0)
  call fminmax_printm('p(mount):',tmp(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,terra(1:nx,1:ny,1:nzm).eq.0)
  call fminmax_printm('t(mount):',t(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,terra(1:nx,1:ny,1:nzm).eq.0)
  call fminmax_printm('tabs(mount):',tabs(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,terra(1:nx,1:ny,1:nzm).eq.0)
  call fminmax_printm('qv(mount):',qv(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm,terra(1:nx,1:ny,1:nzm).eq.0)
 else
  call fminmax_print('u:',u(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('v:',v(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('w:',w(1:nx,1:ny,1:nz),1,nx,1,ny,nz)
  call fminmax_print('p:',tmp(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('pp:',pp(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('t:',t(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('tabs:',tabs(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('qv:',qv(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
 end if
! call averageXY_MPI(w(1:nx,1:ny,1:nz),1,nx,1,ny,nzm,fldm)
! if(masterproc) print*,'w mean:',fldm
 if(dosgs) call sgs_print()
 if(docloud) then
   call fminmax_print('qcl:',qcl,1,nx,1,ny,nzm)
   call fminmax_print('qci:',qci,1,nx,1,ny,nzm)
   call fminmax_print('cld:',cld,1,nx,1,ny,nzm)
   call micro_print()
 end if
 if(doprecip) then
   call fminmax_print('qpl:',qpl,1,nx,1,ny,nzm)
   call fminmax_print('qpi:',qpi,1,nx,1,ny,nzm)
 end if

 if(nstatsteps.eq.-1) goto 2000

 if(dolongwave.or.doshortwave) call fminmax_print('qrad(K/day):',qrad*86400.,1,nx,1,ny,nzm)
 if(dotracers) then
   do k=1,ntracers
    call fminmax_print(trim(tracername(k))//':',tracer(1:nx,1:ny,1:nzm,k),1,nx,1,ny,nzm)
   end do
 end if
 
 if(dochem) call chem_print()

 if(ISLAND) then
  call fminmax_print('shf_land:',shf_land,1,nx,1,ny,1)
  call fminmax_print('lhf_land:',lhf_land,1,nx,1,ny,1)
  call fminmax_print('shf_ocean:',shf_ocean,1,nx,1,ny,1)
  call fminmax_print('lhf_ocean:',lhf_ocean,1,nx,1,ny,1)
 else
  call fminmax_print('shf:',shf_all,1,nx,1,ny,1)
  call fminmax_print('lhf:',lhf_all,1,nx,1,ny,1)
  if(dodns) then
   call fminmax_print('shf_top:',shf_top,1,nx,1,ny,1)
   call fminmax_print('lhf_top:',lhf_top,1,nx,1,ny,1)
  end if
 end if
 call fminmax_print('taux:',fluxbu,1,nx,1,ny,1)
 call fminmax_print('tauy:',fluxbv,1,nx,1,ny,1)
 call fminmax_print('tau:',sqrt(fluxbu**2+fluxbv**2),1,nx,1,ny,1)
 call fminmax_print('skt:',real(sstxy(1:nx,1:ny)+t00),1,nx,1,ny,1)
 call fminmax_printm('sst:',real(sstxy(1:nx,1:ny)+t00),1,nx,1,ny,1,landmask(1:nx,1:ny).eq.0.)
 if(doseaiceevol) call fminmax_print('seaice_h:',seaice_h,1,nx,1,ny,1)
 call fminmax_print('precinst (mm/d):',precinst*86400.,1,nx,1,ny,1)
 if(masterproc) print*,'prec_xy (mm/d):',fld_mean  ! computed way above this line
 call averageXY_SFC(precinst(:,:)*86400.,1,nx,1,ny,fld_mean)
 if(masterproc) print*,'prec_mean (mm/d):',fld_mean
 call averageXY_SFC(evp_all(:,:)*86400.,1,nx,1,ny,fld_mean)
 if(masterproc) print*,'evap_mean (mm/d):',fld_mean
 call averageXY_SFC(shf_all(:,:),1,nx,1,ny,fld_mean)
 if(masterproc) print*,'shf_mean (W/m2):',fld_mean
 call averageXY_SFC(lhf_all(:,:),1,nx,1,ny,fld_mean)
 if(masterproc) print*,'lhf_mean (W/m2):',fld_mean
 if(docup) then
   call fminmax_print('precip_cu:',precip_cu,1,nx,1,ny,1)
   call averageXY_SFC(precip_cu(:,:),1,nx,1,ny,fld_mean)
   if(masterproc) print*,'prec_cu_mean (mm/d):',fld_mean
   call fminmax_print('precip_ls:',precip_ls,1,nx,1,ny,1)
   call averageXY_SFC(precip_ls(:,:),1,nx,1,ny,fld_mean)
   if(masterproc) print*,'prec_ls_mean (mm/d):',fld_mean
   call fminmax_print('wd_cu:',wd_cu,1,nx,1,ny,1)
   call averageXY_SFC(wd_cu(:,:),1,nx,1,ny,fld_mean)
   if(masterproc) print*,'wd_cu_mean (m/s):',fld_mean
 end if

 if(masterproc.and..not.dosmoke.and.nstatsteps.ne.-1) then
   
   print*,'total water budget:'
!   write(*,991) total_water_before !'before (mm):    ',total_water_before
!   write(*,992) total_water_after !'after (mm) :    ',total_water_after
   write(*,*) 'total_water (mm):  ',total_water_after
!   write(*,993) total_water_evap !'evap (mm/day):  ',total_water_evap
!   write(*,994) total_water_prec !'prec (mm/day):  ',total_water_prec
!   write(*,995) total_water_ls !'ls (mm/day):    ',total_water_ls
!   write(*,997) total_water_adv !'adv (mm/day):  ',total_water_adv
!   write(*,998) total_water_after-total_water_before !storage 
!   write(*,996) (total_water_after-(total_water_before+total_water_evap+ &
!                 total_water_ls-total_water_prec+total_water_adv))
!   991 format(' before (mm):       ',F16.11)
!   992 format(' after (mm):        ',F16.11)
!   993 format(' evaporation (mm):  ',F16.11)
!   994 format(' precipitation (mm):',F16.11)
!   995 format(' large-scale (mm):  ',F16.11)
!   996 format(' Imbalance (mm)     ',F16.11)
!   997 format(' advection (mm):    ',F16.11)
!   998 format(' storage (mm)       ',F16.11)
!   print*,' imbalance (rel error):', &
!     (total_water_after-(total_water_before+total_water_evap+total_water_ls &
!     -total_water_prec+total_water_adv))/ &
!     (total_water_after+1.e-15)
   print*,'evap (mm/day):',total_water_evap/dtn*86400.
   print*,'prec (mm/day):',total_water_prec/dtn*86400.
   print*,'ls (mm/day):',total_water_ls/dtn*86400.
   print*,'adv (mm/day):',total_water_adv/dtn*86400.
   print*,'storage (mm/day):',(total_water_after-total_water_before)/dtn*86400.
   print*,'imbalance (mm/day)', &
     (total_water_after-(total_water_before+total_water_evap+total_water_ls+total_water_adv &
             -total_water_prec))/dtn*86400.

 end if
 if(nstep.ne.0.and.SLM) call slm_stepout()
 if(dobuildings) call buildings_stepout()
 if(dodynamicocean) call ocn_stepout(1)

 if(docolumn) then
   print *,'  k      z      Tabs     t      qt      Qn      Qp      REL    u    v'
   do k = nzm,1,-1
    write(6,'(i4,1x,7f8.2,3f8.4,f8.2)')   k,z(k),tabs0(k), &
       t0(k),q0(k)*1.e3, qn0(k)*1.e3,qp0(k)*1.e3, &
       100.*qv0(k)/qsatw(tabs0(k),pres(k)),u0(k),v0(k)
   end do
   print *,'  k      z      Tabs     t      qt      Qn      Qp      REL    u    v'
 end if

end if ! (mod(nstep,nprint).eq.0)

call t_stopf ('stepout')

2000 continue

end
