program mmidentify

! Program to match triggered events in Mermaid SAC records with
! earthquakes in the catalogue obtained from IRIS

! Run this program after automaid from directory processed/
! Please adapt root_directory to where this can be found.
! Note that the root_directory name should end in '/'

! It reads all SAC files in directories such as 
! 452.020-P-06/20180628-19h20m09s/
! If it finds a match in the catalogue, will move the sac, html, png,
! and html files into into root_directory/identified/YEAR.DAY.HRMN
! for an earthquake on Julian DAY, YEAR at time HR:MN
! If there is no match, the files are moved into a directory such as
! root_directory/unidentified/452.020-P-06/20180628-19h20m09s/
! SAC files are named knetwkkstnm.sac etc (where knetwk and kstnm are 
! taken from the SAC file, e.g. MH06.sac)

! It needs:
! -- C-shell script rdeventcat
! -- C-shell script mksaclst
! -- program FetchEvent (from http://service.iris.edu/clients/)
! -- SAC library (normally in /usr/local/sac/lib/sacio.a) to be obtained
!   from IRIS: https://ds.iris.edu/ds/nodes/dmc/software/downloads/sac/

! compile (on my Mac): 
! g7 mmidentify ttak135 findfe march1days date2julian delaz sac bin

! compile (elsewhere):
! gfortran -c mod_ttak135.f90
! gfortran -c mod_findfe.f90
! gfortran -o ~/bin/mmidentify mmidentify.f90 march1days.f90 delaz.f90
!                date2julian.f90 mod_ttak135.o mod_findfe.o sacio.a 

use findfe

implicit none

! event catalogue
! largest number of catalogue entries expected is MC, MD for deep m>5,
! and MX for very large events (m>6.0)
integer, parameter  :: MC=100000, MD=1000, MX=2000
integer :: cyear(MC),cmonth,cday,cjday(MC),chour(MC),cmin(MC) 
integer :: dyear(MD),dmonth,dday,djday(MD),dhour(MD),dmin(MD) 
integer :: xyear(MX),xmonth(MX),xday(MX),xjday(MX),xhour(MX),xmin(MX) 
real*4 :: csec(MC),dsec(MD),xsec(MX)
character*32 :: cevnm(MC),devnm(MD),xevnm(MX)
real*4 :: cevlat(MC),cevlon(MC),cevdept(MC),cevmag(MC)
real*4 :: devlat(MD),devlon(MD),devdept(MD),devmag(MD)
real*4 :: xevlat(MX),xevlon(MX),xevdept(MX),xevmag(MX)

! SAC file parameters
! MS is largest number of processed SAC files expected
! MFL is largest number of different floats
integer, parameter :: NSAC=30000, MS=1000, MFL=100
real*4 :: B,E,O,A,t0,t1,t2,t3,t4,stla,stlo,stdp,gcarc,azis,azie
real*4 :: evdp,evla,evlo,evmag
real*4 :: tP,tp_P,tPP,tS,tPKPdf,tPKPab,tPKPbc
real*4 :: gpslat1(MFL),gpslat2(MFL),gpslon1(MFL),gpslon2(MFL)
real*4 :: gpst1(MFL),gpst2(MFL),gpst
real*4 :: d(NSAC),beg,del,dif1,dif2,x
integer :: fyear(MS),fmonth(MS),fday(MS),fhour(MS),fmin(MS),fsec(MS)
integer :: nsmp,nerr,ndttm1(6),ndttm2(6),nrfe
integer :: nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec
character*8 :: kstnm,knetwk,iztype,floatnm(MFL)
character*16 :: kevnm
character*32 :: kevnm2
character*120 :: sacfile(MS),shortsac(MS)
logical :: missed(MFL,MX)

! Other
character*50 :: floatdir(MFL)
character*180 :: fname, directory, root_directory, newname
character*280 :: syscommand
character*14 :: eventdirectory
character*8 :: today
logical :: hit,ruthere,back
integer :: year,month,day,mon(12)
integer :: i,ifloat,ios,j,k,khit,kmax,kmin,march1days,maxdays,mindays,  &
  n,ndays,ncat,ncatall,ncatc,ncatd,ncatx,mfloats,nfloats,nfiles
real :: tdays,t

data mon/31,28,31,30,31,30,31,31,30,31,30,31/

!----------------------------------------------------------------------
! Preliminaries (site-dependent, please edit)
root_directory = '/Users/guust/data/Mermaids/SPPIM/'

! No site-dependent edits should be needed beyond this line
!----------------------------------------------------------------------

back=.true.
missed=.true.
nfloats=0
gpst1=9999999.
gpst2=-9999999.

! open command file (copies files to identified or unidentified)
open(4,file='move2ident',action='write')
write(4,'(a)') '#!/bin/csh'
call date_and_time(today)
open(7,file='report_identify_'//today)
write(7,'(3a)') 'Year Day     Time      Lat      Lon Dep  Mag',  &
  '  Net  Sta       Lat      Lon  Delta       B       E  Directory', &
  '      Region'
! Get contents of current directory (processed/)
call system('mksaclst')
open(1,file='dumnr')
read(1,*) k
if(k.eq.0) then
  stop 'No SAC files available for identification'
else
  print *,k,' SAC files available for identification'
endif  
open(1,file='dumlstr',iostat=ios)       ! full file names
if(ios.ne.0) stop 'Error in result r of mksaclst'
k=0
do
  k=k+1
  if(k>MS) stop 'Nr of SAC files exceeds MS, increase MS'
  read(1,'(a)',iostat=ios) sacfile(k)
  if(is_iostat_end(ios)) exit
enddo
close(1)
nfiles=k-1

open(2,file='dumlsts',iostat=ios)       ! SAC file name w/o root
if(ios.ne.0) stop 'Error in result s of mksaclst'
do k=1,nfiles
  read(2,'(a)',iostat=ios) shortsac(k)
  if(is_iostat_end(ios)) stop 'Error reading dumlsts'
enddo
close(2)

open(3,file='dumlstt',iostat=ios)       ! times only
if(ios.ne.0) stop 'Error in result t of mksaclst'

! Find min,max date
k=0                       ! SAC file count
mindays=999999
maxdays=-999999
do k=1,nfiles
  read(3,iostat=ios,fmt='(i4,2i2,1x,3i2)') fyear(k),fmonth(k),fday(k),  &
    fhour(k),fmin(k),fsec(k)
  if(is_iostat_end(ios)) stop 'Error reading dumlsts'
  ndays=march1days(fyear(k),fmonth(k),fday(k))
  if(ndays.lt.mindays) then
    kmin=k
    mindays=ndays
  endif  
  if(ndays.gt.maxdays) then
    kmax=k
    maxdays=ndays
  endif  
enddo  
close(3)

year=fyear(kmin)
month=fmonth(kmin)
day=fday(kmin)
! to avoid misses near midnight, start catalogue search one day early
day=day-1
if(day<0) then
  month=month-1
  if(month<0) then
    year=year-1
  endif
  day=mon(month)
  if(mod(year,4).eq.0.and.month.eq.2) day=day+1         ! if leap year
endif  

! Use rdeventcat to search IRIS catalogue for events with magnitude > 2.5
! (this writes three catalogues: eventlist, eventlistdp and eventlist60)
write(syscommand,'(a,i4,a1,i2.2,a1,i2.2,1x,i4,a1,i2.2,a1,i2.2,a)')   &
  'rdeventcat ',year,'-',month,'-',day,fyear(kmax),'-',  &
  fmonth(kmax),'-',fday(kmax),' 2.5'
print *,'Calling ',trim(syscommand)
call system(syscommand)

! Read catalogue excerpts
open(3,file='eventlist')        ! list of all events
k=0
do
  k=k+1
  if(k>MC) stop 'Catalogue excerpt too large, increase MC'
  read(3,fmt='(i4,1x,i2,1x,i2,i3,1x,i2,1x,f6.3,2f11.4,2f7.1,2x,a)',  &
    iostat=ios)  cyear(k),cmonth,cday,chour(k),cmin(k),csec(k),  &
    cevlat(k),cevlon(k),cevdept(k),cevmag(k),cevnm(k)
  if(is_iostat_end(ios)) exit
  call date2julian(cyear(k),cmonth,cday,cjday(k))
enddo
ncatall=k-1
close(3)
open(3,file='eventlistdp')        ! list of strong deep events
k=0
do
  k=k+1
  if(k>MD) stop 'Catalogue excerpt too large, increase MD'
  read(3,fmt='(i4,1x,i2,1x,i2,i3,1x,i2,1x,f6.3,2f11.4,2f7.1,2x,a)',  &
    iostat=ios)  dyear(k),dmonth,dday,dhour(k),dmin(k),dsec(k),  &
    devlat(k),devlon(k),devdept(k),devmag(k),devnm(k)
  if(is_iostat_end(ios)) exit
  call date2julian(dyear(k),dmonth,dday,djday(k))
enddo
ncatd=k-1
close(3)
open(3,file='eventlist60')        ! list of events with M>6.0
k=0
do
  k=k+1
  if(k>MX) stop 'Catalogue excerpt too large, increase MX'
  read(3,fmt='(i4,1x,i2,1x,i2,i3,1x,i2,1x,f6.3,2f11.4,2f7.1,2x,a)',  &
    iostat=ios)  xyear(k),xmonth(k),xday(k),xhour(k),xmin(k),xsec(k),  &
    xevlat(k),xevlon(k),xevdept(k),xevmag(k),xevnm(k)
  if(is_iostat_end(ios)) exit
  call date2julian(xyear(k),xmonth(k),xday(k),xjday(k))
enddo
ncatx=k-1
close(3)

! Now loop through the SAC files
do n=1,nfiles
  ! open SAC file
  call rsac1(trim(sacfile(n)),d,nsmp,beg,del,NSAC,nerr)
  if(nerr.eq.0) then
    print *,'Processing ',trim(sacfile(n))
  else
    print *,'Error ',nerr,' reading ',trim(sacfile(n))
    stop
  endif  
  call getfhv('STLA',stla,nerr)
  call getfhv('STLO',stlo,nerr)
  call getfhv('STDP',stdp,nerr)
  call getfhv('B',B,nerr)
  ! get start time of this record:
  call getnhv('NZYEAR',nzyear,nerr)
  call getnhv('NZJDAY',nzjday,nerr)
  call getnhv('NZHOUR',nzhour,nerr)
  call getnhv('NZMIN',nzmin,nerr)
  call getnhv('NZSEC',nzsec,nerr)
  call getnhv('NZSEC',nzmsec,nerr)
  E=beg+(nsmp-1)*del
  call getkhv('KNETWK',knetwk,nerr)
  call getkhv('KSTNM',kstnm,nerr)
  kstnm(8:8)=' '                ! bug fix
  ifloat=0
  do j=1,nfloats
    if(floatnm(j).eq.kstnm) ifloat=j
  enddo
  if(ifloat.eq.0) then
    nfloats=nfloats+1
    ifloat=nfloats
    floatnm(ifloat)=kstnm
  endif  
  knetwk(8:8)=' '               ! removes non-ascii ending (Python bug?)
  kstnm(8:8)=' '
  ! store seismogram start time in ndttm1
  ndttm1(1)=nzyear
  ndttm1(2)=nzjday
  ndttm1(3)=nzhour
  ndttm1(4)=nzmin
  ndttm1(5)=nzsec
  ndttm1(6)=nzmsec
  ! store start- and end position of float [in this run of mmidentify]
  gpst=tdays(nzyear,nzjday,nzhour)      ! rough time in days since 2017
  if(gpst<gpst1(ifloat)) then
    gpst1(ifloat)=gpst
    gpslat1(ifloat)=stla
    gpslon1(ifloat)=stlo
  endif
  if(gpst>gpst2(ifloat)) then
    gpst2(ifloat)=gpst
    gpslat2(ifloat)=stla
    gpslon2(ifloat)=stlo
  endif

  ! Compare with catalogued events
  hit=.false.
  ncat=1                        ! catalogue no 1 (m>6.0)
  do k=1,ncatx                  ! try strongest events first
    ndttm2(1)=xyear(k)
    ndttm2(2)=xjday(k)
    ndttm2(3)=xhour(k)
    ndttm2(4)=xmin(k)
    ndttm2(5)=xsec(k)
    ndttm2(6)=1000*(xsec(k)-ndttm2(5))
    evla=xevlat(k)
    evlo=xevlon(k)
    evdp=xevdept(k)
    evmag=xevmag(k)
    kevnm2=xevnm(k)
    call getfe(kevnm2,kevnm,nrfe)
    ! find limits of SAC file in terms of seconds since origin time
    call ddttm(ndttm1,ndttm2,dif1)            ! start of file
    dif2=dif1+E-B                             ! end of file
    ! do rough check first
    if(dif2<0.) cycle           ! if quake starts after begin
    if(dif1>1250.) cycle        ! if PKP arrives before begin
    ! more precise check:
    call delaz(stla,stlo,evla,evlo,gcarc,azis,azie)
    call T_ak135(gcarc,evdp,tP,tp_P,tPP,tS,tPKPdf,tPKPab,tPKPbc)
    hit=(tP>dif1).and.(tP<dif2)
    hit=hit.or.(tp_P>dif1.and.tp_P<dif2)
    hit=hit.or.(tS>dif1.and.tS<dif2)
    hit=hit.or.(tPKPdf>dif1.and.tPKPdf<dif2)
    hit=hit.or.(tPKPab>dif1.and.tPKPab<dif2)
    hit=hit.or.(tPKPbc>dif1.and.tPKPbc<dif2)
    khit=k
    if(hit)exit
  enddo  
  if(.not.hit) then             ! try strong deep events next
    ncat=2                      ! deep event catalogue
    do k=1,ncatd                  ! try strongest events first
      ndttm2(1)=dyear(k)
      ndttm2(2)=djday(k)
      ndttm2(3)=dhour(k)
      ndttm2(4)=dmin(k)
      ndttm2(5)=dsec(k)
      ndttm2(6)=1000*(dsec(k)-ndttm2(5))
      evla=devlat(k)
      evlo=devlon(k)
      evdp=devdept(k)
      evmag=devmag(k)
      kevnm2=devnm(k)
      call getfe(kevnm2,kevnm,nrfe)
      ! find limits of SAC file in terms of seconds since origin time
      call ddttm(ndttm1,ndttm2,dif1)            ! start of file
      dif2=dif1+E-B                             ! end of file
      ! do rough check first
      if(dif2<0.) cycle         ! if quake starts after begin
      if(dif1>1250.) cycle      ! if PKP arrives before begin
      ! more precise check:
      call delaz(stla,stlo,evla,evlo,gcarc,azis,azie)
      call T_ak135(gcarc,evdp,tP,tp_P,tPP,tS,tPKPdf,tPKPab,tPKPbc)
      hit=(tP>dif1).and.(tP<dif2)
      hit=(tPP>dif1).and.(tPP<dif2)
      hit=hit.or.(tp_P>dif1.and.tp_P<dif2)
      hit=hit.or.(tS>dif1.and.tS<dif2)
      hit=hit.or.(tPKPdf>dif1.and.tPKPdf<dif2)
      hit=hit.or.(tPKPab>dif1.and.tPKPab<dif2)
      hit=hit.or.(tPKPbc>dif1.and.tPKPbc<dif2)
      khit=k
      if(hit)exit
    enddo  
  endif
  if(.not.hit) then             ! try smaller events next
    ncat=3                      ! small magnitude catalogue
    do k=1,ncatc                  ! try strongest events first
      ndttm2(1)=cyear(k)
      ndttm2(2)=cjday(k)
      ndttm2(3)=chour(k)
      ndttm2(4)=cmin(k)
      ndttm2(5)=csec(k)
      ndttm2(6)=1000*(csec(k)-ndttm2(5))
      evla=cevlat(k)
      evlo=cevlon(k)
      evdp=cevdept(k)
      evmag=cevmag(k)
      kevnm2=cevnm(k)
      call getfe(kevnm2,kevnm,nrfe)
      ! find limits of SAC file in terms of seconds since origin time
      call ddttm(ndttm1,ndttm2,dif1)            ! start of file
      dif2=dif1+E-B                             ! end of file
      ! do rough check first
      if(dif2<0.) cycle         ! if quake starts after begin
      if(dif1>1250.) cycle      ! if PKP arrives before begin
      ! more precise check:
      call delaz(stla,stlo,evla,evlo,gcarc,azis,azie)
      call T_ak135(gcarc,evdp,tP,tp_P,tPP,tS,tPKPdf,tPKPab,tPKPbc)
      hit=(tP>dif1).and.(tP<dif2)
      hit=hit.or.(tp_P>dif1.and.tp_P<dif2)
      hit=hit.or.(tPP>dif1.and.tPP<dif2)
      hit=hit.or.(tS>dif1.and.tS<dif2)
      hit=hit.or.(tPKPdf>dif1.and.tPKPdf<dif2)
      hit=hit.or.(tPKPab>dif1.and.tPKPab<dif2)
      hit=hit.or.(tPKPbc>dif1.and.tPKPbc<dif2)
      khit=k
      if(hit)exit
    enddo  
  endif

  if(hit) then
    B = dif1            ! reset zero time to origin time
    E = dif2
    iztype='IO'
    call setfhv('EVLA',evla,nerr)
    call setfhv('EVLO',evlo,nerr)
    call setfhv('EVDP',evdp,nerr)
    call setfhv('MAG',evmag,nerr)
    if(nrfe>0) call setkhv('KEVNM',kevnm,nerr)
    call setihv('IZTYPE','IO',nerr)
    if(nrfe>0) call setnhv('NEVID',nrfe,nerr)  ! IEVREG gives SAC error...
    call setfhv('O',0.0,nerr)
    ! store earthquake starttime as origin in (new) SAC file
    call setnhv('NZYEAR',ndttm2(1),nerr)
    call setnhv('NZJDAY',ndttm2(2),nerr)
    call setnhv('NZHOUR',ndttm2(3),nerr)
    call setnhv('NZMIN',ndttm2(4),nerr)
    call setnhv('NZSEC',ndttm2(5),nerr)
    call setnhv('NZMSEC',ndttm2(6),nerr)
    call setfhv('B',B,nerr)
    call setlhv('LCALDA',.true.,nerr)
    call setlhv('LEVEN',.true.,nerr)
    call setihv('IFTYPE','ITIME',nerr)
    call setihv('IDEP','iunkn',nerr)
    if(tP>0.) then
      call setfhv('A',tP,nerr)
      call setkhv('KA','P',nerr)
    endif
    if(tp_P>0.) then
      call setfhv('T0',tp_P,nerr)
      call setkhv('KT0','pP',nerr)
    endif
    if(tS>0.) then
      call setfhv('T1',tS,nerr)
      call setkhv('KT1','S',nerr)
    endif
    if(tPKPdf>0.) then
      call setfhv('T2',tPKPdf,nerr)
      call setkhv('KT2','PKPdf',nerr)
    endif
    if(tPKPab>0.) then
      call setfhv('T3',tPKPab,nerr)
      call setkhv('KT3','PKPab',nerr)
    endif
    if(tPKPbc>0.) then
      call setfhv('T4',tPKPbc,nerr)
      call setkhv('KT4','PKPbc',nerr)
    endif
    if(tPP>0.) then
      call setfhv('T5',tPP,nerr)
      call setkhv('KT5','PP',nerr)
    endif
    ! move file to identified directory
    write(eventdirectory,'(i4,a1,i3.3,a1,2i2.2,a1)') ndttm2(1),'.',  &
       ndttm2(2),'.',ndttm2(3),ndttm2(4),'/'
    write(7,'(2i4,i3,1h:,i2.2,1h:,i2.2,2f9.2,i4,f5.1,1x,2a5,2f9.2,f7.1,  &
      2f8.1,2x,a,1x,a)') (ndttm2(i),i=1,5),evla,evlo,int(evdp),evmag,  &
      knetwk(1:4),kstnm(1:4),stla,stlo,gcarc,dif1,dif2,  &
      trim(eventdirectory),trim(kevnm2)

    ! check if identified directory exists, if not create
    fname=trim(root_directory)//'identified/'
    inquire(file=fname,exist=ruthere)
    if(.not.ruthere) then
      call system('mkdir '//trim(fname))
      print *,'Created ',trim(fname)
    endif  

    ! check for eventdirectory and subdirectories
    fname=trim(fname)//trim(eventdirectory)
    inquire(file=fname,exist=ruthere)
    if(.not.ruthere) then
      call system('mkdir '//trim(fname))
      print *,'Created ',trim(fname)
    endif  
    inquire(file=trim(fname)//'sac/',exist=ruthere)
    if(.not.ruthere) then
      call system('mkdir '//trim(fname)//'sac/')
      print *,'Created ',trim(fname)//'sac/'
    endif  
    inquire(file=trim(fname)//'mseed/',exist=ruthere)
    if(.not.ruthere) then
      call system('mkdir '//trim(fname)//'mseed/')
      print *,'Created ',trim(fname)//'mseed/'
    endif  
    inquire(file=trim(fname)//'other/',exist=ruthere)
    if(.not.ruthere) then
      call system('mkdir '//trim(fname)//'other/')
      print *,'Created ',trim(fname)//'other/'
    endif  
    
    ! Move the sac file into the event directory
    newname=trim(fname)//'sac/'//trim(KNETWK)//trim(KSTNM)//'.sac'
    call wsac0(newname,x,d,nerr)
    missed(ifloat,khit)=.false.

    ! add 'quake' file to the sac directory (for clustertomo program)
    inquire(file=trim(fname)//'sac/quake',exist=ruthere)
    if(.not.ruthere) then
      open(3,file=trim(fname)//'sac/quake')
      write(3,*) ndttm2(1),ndttm2(2)
      write(3,*) (ndttm2(i),i=3,6)
      write(3,*) evla,evlo
      write(3,fmt='(f8.1,a,f6.1,2a)') evdp,'  magn:',evmag,' ',kevnm2
      close(3)
    endif  

    i=index(sacfile(n),'.sac')+1
    sacfile(n)(i:i+4)='mseed'
    newname=trim(fname)//'mseed/'//trim(KNETWK)//trim(KSTNM)//'.mseed'
    syscommand='mv '//trim(sacfile(n))//' '//trim(newname)
    write(4,'(a)') trim(syscommand)
    sacfile(n)(i:i+4)='*    '
    syscommand='mv '//trim(sacfile(n))//' '//trim(fname)//'other/'
    write(4,'(a)') trim(syscommand)
    i=index(sacfile(n),'/',back)
    syscommand='cp '//trim(sacfile(n)(1:i))//'*.env '//trim(fname)//'other/'
    write(4,'(a)') trim(syscommand)
    syscommand='cp '//trim(sacfile(n)(1:i))//'*.html '//trim(fname)//'other/'
    write(4,'(a)') trim(syscommand)
    syscommand='cp '//trim(sacfile(n)(1:i))//'*.LOG.h '//trim(fname)//'other/'
    write(4,'(a)') trim(syscommand)
  else
    ! move files to unidentified
    write(eventdirectory,'(i4,a1,i3.3,a1,i2.2,a2,a1)') NZYEAR,'.',  &
       NZJDAY,'.',NZHOUR,'XX','/'
    write(7,'(2i4,i3,1h:,i2.2,1h:,i2.2,1x,a16,11x,2a5,2f9.2,7x,  &
      2f8.1,2x,a)') (ndttm1(i),i=1,5),'    unidentified',  &
      knetwk(1:4),kstnm(1:4),stla,stlo,0.,(nsmp-1)*del,trim(eventdirectory)

    ! check if unidentified directory exists, if not create
    fname=trim(root_directory)//'unidentified/'
    inquire(file=fname,exist=ruthere)
    if(.not.ruthere) then
      call system('mkdir '//trim(fname))
      print *,'Created ',trim(fname)
    endif
    ! check for timed directory:
    fname=trim(fname)//trim(eventdirectory)
    inquire(file=fname,exist=ruthere)
    if(.not.ruthere) then
      call system('mkdir '//trim(fname))
      print *,'Created ',trim(fname)
    endif  
    ! move all files into unidentified/YEAR.JDAY.HRXX
    i=index(sacfile(n),'.sac')+1
    sacfile(n)(i:i+4)='*    '
    syscommand='mv '//trim(sacfile(n))//' '//trim(fname)
    write(4,'(a)') trim(syscommand)
    i=index(sacfile(n),'/',back)
    syscommand='cp '//trim(sacfile(n)(1:i))//'*.env '//trim(fname)
    write(4,'(a)') trim(syscommand)
    syscommand='cp '//trim(sacfile(n)(1:i))//'*.html '//trim(fname)
    write(4,'(a)') trim(syscommand)
    syscommand='cp '//trim(sacfile(n)(1:i))//'*.LOG.h '//trim(fname)
    write(4,'(a)') trim(syscommand)
  endif
enddo  

write(4,*)
write(4,'(a)') '# exit    # uncomment if not sure everything is OK'

write(4,*)
write(4,'(a)') '# Move *kml and *html files to ../processed.float/today'
open(3,file='dumflt',iostat=ios)
mfloats=0
do 
  mfloats=mfloats+1
  read(3,'(a)',iostat=ios) floatdir(mfloats)
  if(ios.ne.0) exit
enddo
mfloats=mfloats-1
if(mfloats.ne.nfloats) then
  print *,'WARNING! Nr of floats detected in SAC files is:',nfloats
  print *,'But nr of subdirectories for floats is:',mfloats,':'
  do i=1,mfloats
    write(6,'(i3,1x,a)') i,trim(floatdir(i))
  enddo
endif  
do i=1,mfloats
  inquire(file='../processed.'//trim(floatdir(i)),exist=ruthere)
  if(.not.ruthere) then
    write(4,'(a)') 'mkdir ../processed.'//trim(floatdir(i))
  endif  
  inquire(file='../processed.'//trim(floatdir(i))//trim(today),  &
    exist=ruthere)
  if(.not.ruthere) then
    write(4,'(a)') 'mkdir ../processed.'//trim(floatdir(i))//trim(today)
  endif
  write(4,'(2a)') 'mv '//trim(floatdir(i))//'*.kml ../processed.',  &
     trim(floatdir(i))//trim(today)
  write(4,'(2a)') 'mv '//trim(floatdir(i))//'*.html ../processed.',  &
     trim(floatdir(i))//trim(today)
enddo
write(4,*)
write(4,'(a)') '# Clean up directory processed/, leave only the report file'
write(4,'(a)') 'ls -d */ > dum'
write(4,'(a)') "sed -e 's/^/rm -fr /' dum > dum1"
write(4,'(a)') 'chmod a+x dum1'
write(4,'(a)') 'dum1'
write(4,'(a)') 'rm dum*'
write(4,'(a)') 'rm event*'
write(4,'(a)') 'rm sedcmds'
write(4,*)
write(4,'(a)') '# Clean up server directory by moving *LOG, *MEM and *vit'
write(4,'(a)') '# to a new directory server.'//trim(today)
write(4,'(a)') 'mkdir ../server.'//trim(today)
write(4,'(a)') 'mv ../server/*.LOG ../server.'//trim(today)
write(4,'(a)') 'mv ../server/*.MEM ../server.'//trim(today)
write(4,'(a)') 'mv ../server/*.vit ../server.'//trim(today)
write(4,'(a)') 'echo Cleanup complete, except move2ident and report files'
close(4)

! Check for possibly missed events. Note:
! It only checks for the current MEM files over their time span
! and it is possible events may already have been processed earlier.
! Linear interpolation for lat/lon between beginning and end of the
! time span is accurate enough to judge if a P or PKP wave could
! have been observed, but inaccurate for all other purposes,
! and the predicted onset times could easily be wrong by many
! seconds (depending on length of the time span).
! If, for some reason, you do the re-identification over all
! existing MEM files the time span is probably too long for this
! analysis to work.
write(7,'(/,a,i3)') 'Nr of floats analysed:',nfloats
write(7,'(/,a)') 'Possibly missed events (onset times very approximate):'
write(7,'(/,2a)') 'Sta  Date                 Time  Depth',  &
  '    Mag  Delta  onset  FE Region'
do k=1,ncatx
  j=0
  do i=1,nfloats
    if(missed(i,k)) then
      call getfe(xevnm(k),kevnm,nrfe)
      ! VERY rough interpolation of float position
      gpst=tdays(xyear(k),xjday(k),xhour(k))
      stla=gpslat1(i)+(gpst-gpst1(i))*(gpslat2(i)-gpslat1(i))/  &
         (gpst2(i)-gpst1(i))
      stlo=gpslon1(i)+(gpst-gpst1(i))*(gpslon2(i)-gpslon1(i))/  &
         (gpst2(i)-gpst1(i))
      call delaz(stla,stlo,xevlat(k),xevlon(k),gcarc,azis,azie)
      call T_ak135(gcarc,xevdept(k),tP,tp_P,tPP,tS,tPKPdf,tPKPab,tPKPbc)
      if(tP>0.) then
        t=tP
      else if(tPKPdf>0.) then
        t=tPKPdf
      else
        t=0.
      endif  
      write(7,'(a4,i5,1x,2i3,1h(,i3,1h),i3,1h:,i2.2,1h:,i2.2,  &
        3f7.1,i7,i4,1x,a)')   &
        floatnm(i)(1:4),xyear(k),xmonth(k),xday(k),  &
        xjday(k),xhour(k),xmin(k),nint(xsec(k)),xevdept(k),xevmag(k),  &
        gcarc,nint(t),nrfe,xevnm(k)
    endif
  enddo 
enddo  

call system('chmod a+x move2ident')
! call system('move2ident')             ! Uncomment for automatic cleanup
print *
print *,'You need to run move2ident to move the files'

end
subroutine T_ak135(dist,dep,tP,tp_P,tPP,tS,tPKPdf,tPKPab,tPKPbc)
! SAC header:               A  t0   t5  t1 t2     t3     t4

! gives AK135 travel time and slowness for P,PKP and S phases with 
! linear interpolation from the tables. Approximate but incredibly 
! robust, it interpolates the AK135 tables linearly. However,
! times can be wrong by several seconds at regional distances!
! If phase is absent, time is set to -12345.0 (SAC convention).

use ttak135
use findfe

implicit none
real*4, intent(in) :: dist, dep
real*4, intent(out) :: tP,tp_P,tPP,tS,tPKPdf, tPKPab,tPKPbc
real*4 :: t

tP=Ptime(dist,dep)
tPKPab=PKPabtime(dist,dep)
tPKPbc=PKPbctime(dist,dep)
tPKPdf=PKPdftime(dist,dep)
tPP=PPtime(dist,dep)
tS=Stime(dist,dep)
tp_P=-12345.
t=smallPtime(dist,dep)
if(tP>0.0 .and. t>0.0) tp_P=tP+t

return

end

subroutine ddttm(ndttm1,ndttm2,diff)

! difference in seconds betwee ndttm1 and ndttm2 (copied from SAC)

implicit none

integer, intent(in) :: ndttm1(6),ndttm2(6)
real*4, intent(out) :: diff
integer :: nday

if(4*(ndttm2(1)/4).eq.ndttm2(1))then
  nday=366
else
  nday=365
endif

diff=        0.001*float(ndttm1(6)-ndttm2(6))   &
    +            float(ndttm1(5)-ndttm2(5))   &
    +        60.000*float(ndttm1(4)-ndttm2(4))   &
    +      3600.000*float(ndttm1(3)-ndttm2(3))   &
    +     86400.000*float(ndttm1(2)-ndttm2(2))   &
    +nday*86400.000*float(ndttm1(1)-ndttm2(1))

return

end

real function tdays(nzyear,nzjday,nzhour)

! returns days (fractional, accurate to one hour) since Jan 1, 2017
! and can be used for GPS location interpolation

implicit none
integer, intent(in) :: nzyear,nzjday,nzhour

if(nzyear<2017) stop 'Error in function tdays'

! leap years every 4 years except century unless divisable by 400
tdays=(nzyear-2017)*365+nzjday+nzhour/24.+(nzyear-2017)/4-1-  &
      (nzyear-2001)/100+(nzyear-2001)/400

return
end


