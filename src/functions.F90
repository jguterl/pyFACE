module functions
use header
implicit none

contains

real(kind=DP) function srcbin(j,k,l)
! Rate of induced-detrapping: detrappping rate from trap n+1 into 1 + n
! due to collisions of filled traps n with incoming ion
integer j, k, l, n
real(kind=DP):: s, rs, crsc
! why?
!      data crsc /3.d-20/
data crsc /0.d0/


srcbin=0.d0
do n=3,nspc,2
 if (l .eq. n) then
  s=crsc*inflx(1)
  rs=s*(1.d0-erf((x(j)-implantation_depth(n))/(sqrt2*implantation_width(n))))
 if ((k .eq. 1) .or. (k .eq. n-1)) then
  srcbin=+rs ! detrap from n -> create 1 free H
             ! detrap from n -> create n-1 empty trap
  return
  elseif (k .eq. n) then
  srcbin=-rs ! detrap from trap n -> filled trap n disapears
  return
  endif
 endif
enddo

return
end function srcbin
!
!






! pre-expoential factor of trapping process
! H + n->n+1
real(kind=16) function kbinar(k,l,m)
integer k, l, m, n
 kbinar=0.0
 do n=2,nspc-1,2
  if(m.eq.1) then
 if (l.eq.n) then
 if ((k.eq.1).or.(k.eq.n)) then
  kbinar=-nu(1)*lambda**3*cvlm
  return
  elseif (k.eq.n+1) then
  kbinar=+nu(1)*lambda**3*cvlm
  return
  endif
! this if for multiple H in traps
!       tmp(1  ,n+1,1)=+rk(n,2)
!       tmp(n  ,n+1,1)=+rk(n,2)
!       tmp(n+1,n+1,1)=-rk(n,2)
      endif
        endif
 enddo

return
end


!activation energy of trapping process
!    (1) + (n) -> (n+1)
real(kind=DP) function ebinar(k,l,m)
integer k, l, m, n

ebinar=0.0
 do n=2,nspc-1,2
 if (m.eq.1) then
 if (l.eq.n) then
 if ((k.eq.1).or.(k.eq.n).or.(k.eq.n+1)) then
  ebinar=etr(n)
  return
  endif
  endif
  endif
!        tmp(1  ,n+1,1)=edtr(n+1)
!        tmp(n  ,n+1,1)=edtr(n+1)
!        tmp(n+1,n+1,1)=edtr(n+1)
 enddo
return
end




!prex-exponential factor for thermal (activated) detrapping
!(n+1)-> (1) + (n)
real(kind=DP) function ktherm(k,l)
implicit none
integer k, l, n
 ktherm=0.d0
 do n=3,nspc,2
 if (l .eq. n) then
 if ((k .eq. 1) .or. (k .eq. n-1)) then
  ktherm=+nu(1)
  return
  elseif (k .eq. n) then
  ktherm=-nu(1)
  return
  endif
  endif
 enddo

return

end function ktherm


!> Activation energy for thermal (activated) detrapping
!! @param aggr information about the aggregates
!! @todo Handle special case
!! (n) -> (1) +(n-1)
real(kind=DP) function etherm(k,l)
implicit none

integer ::k,l !< index of species
integer::n
etherm=0
 do n=3,nspc,2
     if (l .eq. n) then
         if ((k .eq. 1) .or. (k .eq. n-1) .or. (k .eq. n)) then
  etherm=edtr(n)
  return
  endif
  endif
 enddo
return
end





       real(kind=DP) function integrale_src_T()
integrale_src_T=0d0
      end function integrale_src_T

      real function integrale_T()
integer ::j
real ::s

s=0d0
      do j=0,ngrd-1
          s=s+0.5d0*(temp(ndt,j)+temp(ndt,j+1))*dx(j)
      enddo
      integrale_T=s*rhocp
      end function integrale_T




 subroutine print_source(ifile)
integer j,k,ifile
write(ifile,*) "#source rate (function source(j,k))"
write(ifile,*)"xgrid ",(namespc(k)//" ",k=1,nspc)
do j=1,ngrd
write(ifile,*) x(j),(src_profile(j,k),k=1,nspc)
enddo
end subroutine print_source

end module functions
