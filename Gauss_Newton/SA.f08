program  SMA
implicit none
real p, p0, pmin, pmax, x, theta, v
real T, deltT, E0, E, deltE, cte, d, din, pv, x0, prob,f
allocatable x(:), pv(:), p0(:), pmin(:), pmax(:), d(:), din(:),f(:), p(:)
integer i , j, m, n , ntemp , nmodel , k , xini , dx
!========================================
! Modelo Inicial- Parametros
!=========================================

        n=3 ; m=50
        theta = 5.0
        cte=acos(-1.0)/180
        x0=-1200
        dx=50
        allocate (pv(n), p0(n), pmin(n), pmax(n), d(m), din(m), x(m), f(n), p(n))
        pv(1)=600.; pv(2)=1500.; pv(3)=theta*cte

        open(unit=10, file="modelo.dat")
        open(unit=20, file="modelo_inicial.dat")
        open(unit=30, file="modelo_final.dat")
        open(unit=40, file="parametros.dat")

        do j=1, m
           x(j)=x0+(j-1)*dx
        end do
        
        call forward(m, x, n, pv,d)
        
        do i=1, m
           write(10,*) x(i), d(i)
        end do


!================= Modelo inicial =============================
p0(1) = pv(1)*0.5 ;  p0(2) = pv(2)*0.5 ;  p0(3) = pv(3)*0.5
!===============================================================
        
        call forward( m, x, n, p0, din)
        
        do i=1, m
            write(20,*) x(i), din(i)
        enddo
!===============================================================
!            Energia Inicial
!==============================================================
        
        call energia(m, x, n, p0, d, E0)
       
        pmin(1)=100   ; pmax(1) = 5000
        pmin(2)=p0(2) ; pmax(2) = 6000
        pmin(3)=p0(3) ; pmax(3) = 6*p0(3)
        
        T=1000
        
        do i=1, 250 
     
           T=0.95*T
           do j=1, 1000

              call get_newmodel( n, pmin, pmax, p)
              call energia(m, x, n, p, d, E)

              deltE= E-E0

              if (deltE.le. 0) then 
                 p0=p
                 E0=E
              else

              call random_number(v)
              prob=exp(-deltE/T)

              if (prob>v) then
                 p0=p
                 E0=E
              end if
              end if    
            end do

            write(*,*) i, "E0=", E0, "T=", T, "prob=", prob
            open(60, file='erro_sma.dat')
            write(60,*) i, E0
         
            if ((E0.lt.0.001)) then
               if(T.lt.0.01) then
                  write(*,*) " resfriou", T
                  exit
               end if
            end if  
        end do
        
        do i=1, n
           write(40,*) p0(i)
        end do
        
        write(*,*) " h=", p0(1), " velocidade=", p0(2),"theta=", p0(3)/cte
     
        call forward( m, x, n, p0, din)
        
        do i=1, m
           write(30, *) x(i), din(i)
        end do
        
        close(10)
        close(20)
        close(30)
        close(40)
        close(60)
        deallocate( pv, p0, pmin, pmax, d, din, x, f)
        end

!===============================================================
!       Subroutina forward
!===============================================================

        subroutine forward(m, x, n, p, f)
        dimension x(m), p(n), f(m)
        h=p(1); v=p(2); theta=p(3)
        a1=2.0*h*sin(theta)
        a2=(2.0*h*cos(theta))**2
        
        do i=1, m
           f(i)=sqrt((x(i)+a1)**2.0 + a2)/v
        end do
        return
        end
!===============================================================
!       Subroutina Calcula energia
!===============================================================

        subroutine energia (m, x, n, p, time, f)
        dimension x(m), p(n), time(m)
        allocatable tinv(:)
        allocate (tinv(m))
        call forward(m, x, n, p, tinv)
        f=dot_product(time-tinv, time-tinv)
        deallocate(tinv)
        return 
        end
!===============================================================
!       Subroutina Gera um Novo modelo
!===============================================================
        

        subroutine get_newmodel(n, pmin, pmax, p)
        dimension p(n), pmin(n), pmax(n)
        
        do j=1, n
           call random_number(v)
           p(j) = pmin(j)+v*(pmax(j)-pmin(j))
        end do 
        return
        end

















