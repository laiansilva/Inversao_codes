module inv
use aux
contains

!==========================================================================
!       Subrountina para inverter os dados
!==========================================================================      
        subroutine inv_GN (m, x, d, n, p0, Niter, pk, fk)
        dimension x(m), d(m), p0(n), pk(n), fk(m)
        Allocatable deltp(:), deltd(:), gk(:, :)
        Allocate (deltp(n), deltd(m), gk(m, n))
        
        call forward (m, x, n, p0, fk)
        deltd=d-fk    
        tol1=1.e-7
        tol2=150
        pk=p0
        ikey=1
        call dot (m, deltd, deltd, q0)
        write (*,*) 'q0', sqrt(q0)
        k=1
        open(60, file='erro_dado.dat')
        do while (k.le. Niter)
            !do while (ikey.eq.1)

                call sensimx(m, x, n, pk, fk, gk) 
                call cg_Xhy (m, n, gk, deltp, deltd) ! estou aqui 

                write(*,*) deltp
                pk = pk + deltp

                call forward(m, x, n, pk, fk)
                deltd = d-fk

               ! do i=1, m
                 ! write(*,*) i, deltd(i), fk(i)
               ! end do

                call dot(m, deltd, deltd, q1)
                  write(*, *)'iteracao=', k, 'erro=',sqrt(q1)
                  write(60, *) k, sqrt(q1)      
                  write(*, *) ' ============================'
                if (q1.lt.tol1) then
                    k=100
                else
                    if (sqrt(q1).gt.tol2) then
                       write(*, *) 'comecou a divergir'
                       stop             
                    else
                       q0=q1
                    end if
                end if
                k=k+1
            !end do
        end do
        close(60)
        deallocate (deltd, deltp, gk)
        return
        end

!
!==============================================================
!   Subroutina cg_Xhy
!==============================================================
!
        subroutine cg_Xhy(m, n, X, h, y)
        dimension X(m, n), y(m), h(n)
        allocatable A(:,:), b(:)
        allocate (A(n,n), b(n))

        call mult_xtx_xty(m, n, X, y, A, b)
        
        call cg_Ahb(n, A, h, b)

        deallocate(A, b)
        return
        end
        
!===============================================================
!       Subroutina para calculo do gradiente
!===============================================================  
        subroutine cg_Ahb(n, A, h, b)
        dimension A(n,n), h(n), b(n)
        allocatable g(:), p(:), v(:)
        allocate (g(n), p(n), v(n))

        h=0; g=-b ; v=g
        tol=1.e-10
        ikey=1 ; k=1

        do while (ikey.eq.1)
           call matmult (n, n, A, v, p)
           call dot (n, v, p, t)
           call dot (n, v, g, s)
           alfa=-s/t
           h=h+alfa*v
          ! do i=1, n
          !    write(*,*) i, h(i)
          ! end do
          ! write(*,*) "======================================="
           g=g+alfa*p
           call dot(n, p, g, pg)
           beta= -pg/t
           v= g+beta*v
           call dot(n, g, g, gg)
           if(gg.le.tol.or.k.eq.n) ikey=0
           ! write(*,*) k, gg
           k=k+1
        end do
        deallocate(g, v, p)
        return
        end
!=============================================================
!               Subroutina multiplicação de matriz
!=============================================================
        subroutine mult_xtx_xty(m, n, x, y, a, b)
        dimension x(m,n), y(m), a(n,n), b(n)
        do j=1, n
           do i=j, n
              call dot(m, x(:,j), x(:,i), a(i,j))
              a(j,i) = a(i,j)
            end do
        call dot (m, x(1,j), y, b(j))
        end do
        return
        end
end module inv
