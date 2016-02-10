c-----------------------------------------------------------------------
      subroutine cg(x,f,g,c,r,w,p,z,n,niter,flop_cg,n_shared,
     $     x_src,x_dest,y_src,y_dest,
     $     z_src,z_dest,
     $     npx,npy,npz,mx,my,mz)
      include 'SIZE'

c     Solve Ax=f where A is SPD and is invoked by ax()
c
c     Output:  x - vector of length n
c
c     Input:   f - vector of length n
c     Input:   g - geometric factors for SEM operator
c     Input:   c - inverse of the counting matrix
c
c     Work arrays:   r,w,p,z  - vectors of length n
c
c     User-provided ax(w,z,n) returns  w := Az,  
c
c     User-provided solveM(z,r,n) ) returns  z := M^-1 r,  
c

      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)
      parameter (lt=lx1*ly1*lz1*lelt)
      real ur(lt),us(lt),ut(lt),wk(lt)

      real x(n),f(n),r(n),w(n),p(n),z(n),g(1),c(n)
      integer count
      integer npx,npy, npz, mx, my,mz
      
      integer n_shared
      integer src,dest,ierr,status99
      character*1 ans
      integer x_src,x_dest,y_src,y_dest,z_src,z_dest

      
      pap = 0.0

c     set machine tolerances
      one = 1.
      eps = 1.e-20
      if (one+eps .eq. one) eps = 1.e-14
      if (one+eps .eq. one) eps = 1.e-7

      rtz1=1.0

      call rzero(x,n)
      call copy (r,f,n)
      call maskit (r,cmask,nx1,ny1,nz1) ! Zero out Dirichlet conditions

c     Testing cmask
      rnorm = sqrt(glsc3(r,c,r,n))
      iter = 0
      if (nid.eq.0)  write(6,6) iter,rnorm
      miter = niter
c     call tester(z,r,n)  
      call mytimer(3)
      do iter=1,miter
         call solveM(z,r,n)    ! preconditioner here

         rtz2=rtz1                                                       ! OPS
         rtz1=glsc3(r,c,z,n)   ! parallel weighted inner product r^T C z ! 3n

         beta = rtz1/rtz2
         if (iter.eq.1) beta=0.0
         call add2s1(p,z,beta,n)                                         ! 2n

         call ax(w,p,g,ur,us,ut,wk,n,n_shared,
     $        x_src,x_dest,y_src,y_dest,
     $        z_src,z_dest,
     $        npx,npy,npz,mx,my,mz)

         pap=glsc3(w,c,p,n)                                              ! 3n

         alpha=rtz1/pap
         alphm=-alpha
         call add2s2(x,p,alpha,n)                                        ! 2n
         call add2s2(r,w,alphm,n)                                        ! 2n

         rtr = glsc3(r,c,r,n)                                            ! 3n
         if (iter.eq.1) rlim2 = rtr*eps**2
         if (iter.eq.1) rtr0  = rtr
         rnorm = sqrt(rtr)
c        if (nid.eq.0.and.mod(iter,100).eq.0) 
c    $      write(6,6) iter,rnorm,alpha,beta,pap
    6    format('cg:',i4,1p4e12.4)
c        if (rtr.le.rlim2) goto 1001

      enddo
      call mytimer(2)
 1001 continue

      if (nid.eq.0) write(6,6) iter,rnorm,alpha,beta,pap

      flop_cg = flop_cg + iter*15.*n

      return
      end
c-----------------------------------------------------------------------        
      subroutine mytimer(flag)
      include 'mpif.h'
      include 'SIZE'
      include 'TOTAL'
      integer flag
      real*8 time0,time1
      real*8 tmp
      save time0,time1,tmp
      data time0/0.0/, time1/0.0/

      if(flag.eq.0) then
         tmp = mpi_wtime()
      endif
      if(flag.eq.1) then
         time1 = time1+ mpi_wtime() - tmp
      endif
      if(flag.eq.2) then
         if(nid.eq.0) then
            write(*,'(A,e12.4)') 'comm time ', time1
            write(777,'(A,1pe12.4)') 'comm time ',time1
         endif
      endif
      if(flag.eq.3) then
         time1 = 0.0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine solveM(z,r,n)
      include 'INPUT'
      real z(n),r(n)

      nn = n
      call h1mg_solve(z,r,nn)

      return
      end
c-----------------------------------------------------------------------
      subroutine ax(w,u,gxyz,ur,us,ut,wk,n,n_shared,
     $     x_src,x_dest,y_src,y_dest,
     $     z_src,z_dest,
     $     npx,npy,npz,mx,my,mz)
      include 'mpif.h'
      include 'SIZE'
      include 'TOTAL'
c      include 'omp_lib.h'
      integer n_shared
c     number of processes per ox, oy axes, number elements per ox,oy 
      integer npx,npy,npz,mx,my,mz
      real w(nx1*ny1*nz1,nelt),u(nx1*ny1*nz1,nelt)
      real gxyz(2*ldim,nx1*ny1*nz1,nelt)

      parameter (lt=lx1*ly1*lz1*lelt)
      real ur(lt),us(lt),ut(lt),wk(lt)
      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)
      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
      integer e

      integer src,dest,ierr1,status1,ierr2,status2
      integer status
      integer size
      real sum
      integer  e_x, e_y, e_z, e_right,e_left
      integer e_top, e_bottom,e_back,e_front
      integer i,j
      real*8 time1,time2
      
      integer x_src,x_dest,y_src,y_dest,z_src,z_dest


c     !$OMP PARALLEL DO PRIVATE(e,ur,us,ut,wk) SHARED(nelt,w,u,gxyz)  
c     !$OMP& SCHEDULE(STATIC)        
      do e=1,nelt                                ! ~
         call ax_e( w(1,e),u(1,e),gxyz(1,1,e)    ! w   = A  u
     $                             ,ur,us,ut,wk) !  L     L  L
      enddo                                      ! 
c     !$OMP END PARALLEL DO   

c      call rone(w,n)
c      time1 = mpi_wtime()
      call mytimer(0)
      call sync_xyz(w,n_shared,
     $     x_src,x_dest,y_src,y_dest,
     $     z_src,z_dest,
     $     npx,npy,npz,mx,my,mz)
      call mytimer(1)
c      time2 = mpi_wtime()-time1
c      if(nid.eq.0) then
c         write(*,'(A,e12.4)') 'time sync ', time2 
c      endif

c      call dssum(w)         ! Gather-scatter operation  ! w   = QQ  w
                                                           !            L
c      sum = glsum(w,n)
c      if(nid.eq.0) then
c      write(*,*) 'glsum = ', sum 
c      endif
c      call exit()
c      if(nid.eq.0) then 
c         do i=1,n
c            write(*,'(A,I3,A,e12.4)') 'w[ ',i,' ]= ', w(i,1)
c         enddo
c      endif

      call add2s2(w,u,.1,n)   !2n
      call maskit(w,cmask,nx1,ny1,nz1)  ! Zero out Dirichlet conditions

      nxyz=nx1*ny1*nz1
      flop_a = flop_a + (19*nxyz+12*nx1*nxyz)*nelt

      return
      end
c------------------------------------------------------------------------- 
      subroutine ax1(w,u,n)
      include 'SIZE'
      real w(n),u(n)
      real h2i
  
      h2i = (n+1)*(n+1)  
      do i = 2,n-1
         w(i)=h2i*(2*u(i)-u(i-1)-u(i+1))
      enddo
      w(1)  = h2i*(2*u(1)-u(2  ))
      w(n)  = h2i*(2*u(n)-u(n-1))

      return
      end
c-------------------------------------------------------------------------
      subroutine ax_e(w,u,g,ur,us,ut,wk) ! Local matrix-vector product
      include 'SIZE'
      include 'TOTAL'

      parameter (lxyz=lx1*ly1*lz1)
      real ur(lxyz),us(lxyz),ut(lxyz),wk(lxyz)
      real w(nx1*ny1*nz1),u(nx1*ny1*nz1),g(2*ldim,nx1*ny1*nz1)


      nxyz = nx1*ny1*nz1
      n    = nx1-1

      call local_grad3(ur,us,ut,u,n,dxm1,dxtm1)

      do i=1,nxyz
         wr = g(1,i)*ur(i) + g(2,i)*us(i) + g(3,i)*ut(i)
         ws = g(2,i)*ur(i) + g(4,i)*us(i) + g(5,i)*ut(i)
         wt = g(3,i)*ur(i) + g(5,i)*us(i) + g(6,i)*ut(i)
         ur(i) = wr
         us(i) = ws
         ut(i) = wt
      enddo

      call local_grad3_t(w,ur,us,ut,n,dxm1,dxtm1,wk)

      return
      end
c-------------------------------------------------------------------------
      subroutine local_grad3(ur,us,ut,u,n,D,Dt)
c     Output: ur,us,ut         Input:u,n,D,Dt
      real ur(0:n,0:n,0:n),us(0:n,0:n,0:n),ut(0:n,0:n,0:n)
      real u (0:n,0:n,0:n)
      real D (0:n,0:n),Dt(0:n,0:n)
      integer e

      m1 = n+1
      m2 = m1*m1

      call mxm(D ,m1,u,m1,ur,m2)
      do k=0,n
         call mxm(u(0,0,k),m1,Dt,m1,us(0,0,k),m1)
      enddo
      call mxm(u,m2,Dt,m1,ut,m1)

      return
      end
c-----------------------------------------------------------------------
      subroutine local_grad3_t(u,ur,us,ut,N,D,Dt,w)
c     Output: ur,us,ut         Input:u,N,D,Dt
      real u (0:N,0:N,0:N)
      real ur(0:N,0:N,0:N),us(0:N,0:N,0:N),ut(0:N,0:N,0:N)
      real D (0:N,0:N),Dt(0:N,0:N)
      real w (0:N,0:N,0:N)
      integer e

      m1 = N+1
      m2 = m1*m1
      m3 = m1*m1*m1

      call mxm(Dt,m1,ur,m1,u,m2)

      do k=0,N
         call mxm(us(0,0,k),m1,D ,m1,w(0,0,k),m1)
      enddo
      call add2(u,w,m3)

      call mxm(ut,m2,D ,m1,w,m1)
      call add2(u,w,m3)

      return
      end
c-----------------------------------------------------------------------
      subroutine maskit(w,pmask,nx,ny,nz)   ! Zero out Dirichlet conditions
      include 'SIZE'
      include 'PARALLEL'

      real pmask(-1:lx1*ly1*lz1*lelt)
      real w(1)
      integer e

      nxyz = nx*ny*nz
      nxy  = nx*ny
      if(pmask(-1).lt.0) then
        j=pmask(0)
        do i = 1,j
           k = pmask(i)
           w(k)=0.0
        enddo
      else
c         Zero out Dirichlet boundaries.
c
c                      +------+     ^ Y
c                     /   3  /|     |
c               4--> /      / |     |
c                   +------+ 2 +    +----> X
c                   |   5  |  /    /
c                   |      | /    /
c                   +------+     Z   
c

        nn = 0
        do e  = 1,nelt
          call get_face(w,nx,e)
          do i = 1,nxyz
             if(w(i).eq.0) then
               nn=nn+1
               pmask(nn)=i
             endif
          enddo
        enddo     
        pmask(-1) = -1.
        pmask(0) = nn
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine masko(w)   ! Old 'mask'
      include 'SIZE'
      real w(1)

      if (nid.eq.0) w(1) = 0.  ! suitable for solvability

      return
      end
c-----------------------------------------------------------------------
      subroutine masking(w,nx,e,x0,x1,y0,y1,z0,z1)
c     Zeros out boundary
      include 'SIZE'
      integer e,x0,x1,y0,y1,z0,z1
      real w(nx,nx,nx,nelt)
      
c       write(6,*) x0,x1,y0,y1,z0,z1
      do k=z0,z1
      do j=y0,y1
      do i=x0,x1
          w(i,j,k,e)=0.0
      enddo
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
c     conflict if only fortran
c      subroutine tester(z,r,n)
c     Used to test if solution to precond. is SPD
c      real r(n),z(n)
c
c      do j=1,n
c         call rzero(r,n)
c         r(j) = 1.0
c         call solveM(z,r,n)
c         do i=1,n
c            write(79,*) z(i)
c         enddo
c      enddo
c      call exitt0
c      return
c      end
c-----------------------------------------------------------------------
      subroutine get_face(w,nx,ie)
c     zero out all boundaries as Dirichlet
c     to change, change this routine to only zero out 
c     the nodes desired to be Dirichlet, and leave others alone.
      include 'SIZE'
      include 'PARALLEL'
      real w(1)
      integer nx,ie,nelx,nely,nelz
      integer x0,x1,y0,y1,z0,z1
      
      x0=1
      y0=1
      z0=1
      x1=nx
      y1=nx
      z1=nx
      
      nelxy=nelx*nely
      ngl = lglel(ie)        !global element number

      ir = 1+(ngl-1)/nelxy   !global z-count
      iq = mod1(ngl,nelxy)   !global y-count
      iq = 1+(iq-1)/nelx     
      ip = mod1(ngl,nelx)    !global x-count

c     write(6,1) ip,iq,ir,nelx,nely,nelz, nelt,' test it'
c  1  format(7i7,a8)

      if(mod(ip,nelx).eq.1.or.nelx.eq.1)   then  ! Face4
         x0=1
         x1=1
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(ip,nelx).eq.0)               then   ! Face2
         x0=nx
         x1=nx
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      x0=1
      x1=nx
      if(mod(iq,nely).eq.1.or.nely.eq.1) then    ! Face1
         y0=1
         y1=1
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(iq,nely).eq.0)              then    ! Face3
         y0=nx
         y1=nx
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      y0=1
      y1=nx
      if(mod(ir,nelz).eq.1.or.nelz.eq.1) then    ! Face5
         z0=1
         z1=1
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(ir,nelz).eq.0)              then    ! Face6
         z1=nx
         z0=nx
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      return
      end
c-----------------------------------------------------------------------
