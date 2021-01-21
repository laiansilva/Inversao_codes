Program inversao_GN
use aux
use inv
integer i, j
real x0, dx, dt
Allocatable x(:), p0(:), pv(:), pinv(:), d(:), f0(:), finv(:)
m=50;n=3
Allocate(x(m), p0(n), pv(n) , pinv(n), d(m), f0(m), finv(m))
open (10 , file ='inversao_GN.dat')
open (20 , file ='mod_GN.dat')
open (30 , file ='modpert_GN.dat')

Niter=15

!==============================================================
!               Parametros do modelo
!==============================================================
x0=-1200  ; dx=50 ; dt=1

do i=1, m
   x(i)= x0 + (i-1)*dx
end do

pv(1) = 600 ; pv(2) = 1500 ; pv(3) = 5
cte = acos(-1.0)/180
pv(3) = pv(3)*cte

call  forward (m, x, n, pv, d)! d = dados  sinteticos

!==============================================================


!================================= Modelo pertubado==============
p0(1) = pv(1)*0.4 ;  p0(2) = pv(2)*2 ;  p0(3) = pv(3)*3
!===============================================================

call inv_GN(m, x, d, n, p0, Niter, pinv, finv) ! inversão

call forward (m, x, n, p0, f0) ! calculo modelo pertubado

do i=1, m
  ! write (*, *) i, d(i), f0(i), finv(i)
   write (10, *) x(i), finv(i)
   write (20, *) x(i), d(i)
   write (30, *) x(i), f0(i)
end do 

print*, 'h_pert=', p0(1), 'vel_pert =', p0(2), 'alfa_pert=', p0(3)/cte 
print*,'parametros  reais' ,'  ', 'h=', pv(1), 'vel=' , pv(2), 'alfa=', pv(3)/cte
print*,'parametros  invertidos' ,'  ', 'h_inv=', pinv(1), 'vel_inv=' , pinv(2), 'alfa_inv=', pinv(3)/cte

deallocate (x , p0, pv, pinv, d, f0, finv)
stop
close (10)
close (20)
close (30)
end program



























