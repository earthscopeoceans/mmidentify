integer function march1days(y,m,d)

! returns number of days since March 1, 1900.
! method: http://alcor.concordia.ca/~gpkatch/gdate-method.html

integer :: y,m,d,mm,yy,dsincemarch
mm=mod(m+9,12)                  ! set March=0, April=1,...,Feb=11
dsincemarch=(mm*306+5)/10       ! days since March 1
yy=y-1900-mm/10                 ! yy equals year-1900 unless Jan, Feb 
march1days=365*yy+yy/4-yy/100+dsincemarch+(d-1)
if( (y.eq.2000.and.mm.gt.2) .or. (y.gt.2000) ) march1days=march1days+1
return
end
