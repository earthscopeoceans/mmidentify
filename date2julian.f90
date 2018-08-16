subroutine date2julian(year,month,day,jday)

! Converts year,month,day into Julian day jday

implicit none
character(len=12) :: arg

integer :: mday(12),kday(12),i,k,m,n,day,jday,month,year,leap
data mday/31,28,31,30,31,30,31,31,30,31,30,31/

leap=0
mday(2)=28
if(mod(year,4).eq.0.and.(mod(year,100).ne.0.or.mod(year,400).eq.0)) leap=1
if(leap.eq.1) mday(2)=29
kday(1)=0
do i=2,12
  kday(i)=kday(i-1)+mday(i-1)
enddo
jday=kday(month)+day

return
end

subroutine julian2date(jday,year,month,day)

! Converts Julian day jday and year into month,day

implicit none
character(len=12) :: arg

integer :: mday(12),kday(12),i,k,m,n,day,jday,month,year,leap
data mday/31,28,31,30,31,30,31,31,30,31,30,31/

leap=0
mday(2)=28
if(mod(year,4).eq.0.and.(mod(year,100).ne.0.or.mod(year,400).eq.0))  &
  mday(2)=29
m=1
day=jday
do while(m.lt.12.and.jday-mday(m)>0)
  day=day-mday(m)
  m=m+1
enddo
month=m

return
end
