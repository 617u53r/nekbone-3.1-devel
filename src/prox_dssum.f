c-----------------------------------------------------------------------
      subroutine dssum(f)
      include 'SIZE'
      include 'TOTAL'
      real f(1)

c     call nekgsync()
      call adelay
c      call gs_op(gsh,f,1,1,0)  ! Gather-scatter operation  ! w   = QQ  w

      return
      end
c-----------------------------------------------------------------------
      subroutine proxy_setupds(gs_handle,nx,glo_num,n_shared)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      integer gs_handle,dof
      integer*8 glo_num(lx1*ly1*lz1*lelt)
      integer n_shared
      
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      t0 = dnekclock()

      call set_vert_box(glo_num,nx,n_shared)  

c      call set_vert_box(glo_num,nx) ! Set global-to-local map
c      call outmat_glo_num(glo_num,nx)
c      call outmat_glo_num_general(glo_num, nx)

      ntot      = nx*nx*nx*nelt   ! assumes nx=ny=nz
c     disable c code
c      call gs_setup(gs_handle,glo_num,ntot,nekcomm,mp) ! Initialize gather-scatter
      dof = ntot *mp
      t1 = dnekclock() - t0
c     if (nid.eq.0) then
c        write(6,1) t1,gs_handle,nx,dof
c   1    format('   setupds time',1pe11.4,' seconds ',2i3,i12)
c     endif

      return
      end
c-----------------------------------------------------------------------
c      subroutine set_vert_box(glo_num,nx)
      subroutine set_vert_box(glo_num,nx,n_shared)

c     Set up global numbering for elements in a box

      include 'SIZE'
      include 'PARALLEL'

      integer*8 ii,kg,jg,ig ! The latter 3 for proper promotion

      integer e,ex,ey,ez,eg

      integer*8 glo_num(lx1*ly1*lz1*lelt)
      integer n_shared

      integer i_right,i_left,j_top,j_bottom,i_front,i_back
      i_right = 1
      i_left = 1
      j_top = 1 
      j_bottom = 1 
      i_front = 1 
      i_back = 1


      nn = nx-1  ! nn := polynomial order
      do e=1,nelt
        eg = lglel(e)                              
        call get_exyz(ex,ey,ez,eg,nelx,nely,nelz)  
        do k=0,nn
        do j=0,nn
        do i=0,nn
           kg = nn*(ez-1) + k                     
           jg = nn*(ey-1) + j                     
           ig = nn*(ex-1) + i
           ii = 1 + ig + jg*(nn*nelx+1) + kg*(nn*nelx+1)*(nn*nely+1) 
           ll = 1 + i + nx*j + nx*nx*k + nx*nx*nx*(e-1)
           glo_num(ll) = ii

           if ((i.eq.0).and.(e.eq.1)) then
              l_face_left(i_left) = ll
              i_left = i_left + 1
           endif
           if ((i.eq.nn).and.(e.eq.1)) then
              l_face_right(i_right)=ll
              i_right = i_right + 1
           endif

           if((j.eq.0).and.(e.eq.1)) then
              l_face_y_bottom(j_bottom) = ll
              if(nid.eq.0)then                                                         
              endif

              j_bottom = j_bottom + 1
           endif
           if((j.eq.nn).and.(e.eq.1))then
              l_face_y_top(j_top) = ll
              j_top = j_top + 1
           endif
           
           if((k.eq.0).and.(e.eq.1))then
              l_face_z_front(i_front) = ll
              if(nid.eq.0)then                                                                
              endif
              i_front = i_front + 1
           endif
           if((k.eq.nn).and.(e.eq.1))then
              l_face_z_back(i_back) = ll
              if(nid.eq.0) then
              endif
              i_back = i_back + 1
           endif
           
        enddo
        enddo
        enddo
      enddo
      n_shared = i_left - 1
c      write(*,'(A,I2)') 'setvertbox nshared = ',n_shared
      return
      end
c-----------------------------------------------------------------------        
      subroutine faces_ids(nx,mx,my,mz)
      include 'SIZE'
      include 'PARALLEL'

      integer nx
      integer i, j,k,m,q
      q = 1
      do m=1,mz
         do k = 1,my
            do i = 1,nx
               do j = 1,nx
                  face_left(q) = 1+(j-1)*nx+(i-1)*nx*nx+
     $                 (k-1)*nx*nx*nx*mx + (m-1)*nx*nx*nx*mx*my
                  face_right(q) = nx*nx*nx*(mx-1)+nx+
     $                 (j-1)*nx+(i-1)*nx*nx+
     $                 (k-1)*nx*nx*nx*mx + (m-1)*nx*nx*nx*mx*my
                  q = q+1
               enddo
            enddo
         enddo
      enddo

      q = 1
      do m=1,my
         do k = 1,mx
            do i = 1,nx
               do j = 1,nx
                  face_z_front(q)=j+(i-1)*nx+(k-1)*nx*nx*nx+
     $                 (m-1)*nx*nx*nx*mx
                  face_z_back(q) = nx*nx*nx*mx*my*(mz-1)+
     $                 nx*nx*(nx-1)+j+(i-1)*nx+(k-1)*nx*nx*nx+
     $                 (m-1)*nx*nx*nx*mx
                  q = q+1
               enddo
            enddo
         enddo
      enddo

      q = 1
      do m=1,mz
         do k = 1,mx
            do i = 1,nx
               do j = 1,nx
                  face_y_bottom(q) = j+(i-1)*nx*nx+(k-1)*nx*nx*nx+
     $                 (m-1)*nx*nx*nx*mx*my
                  face_y_top(q) = nx*nx*nx*mx*(my-1)+nx*(nx-1)+
     $                 j+(i-1)*nx*nx+(k-1)*nx*nx*nx+
     $                 (m-1)*nx*nx*nx*mx*my
                  q = q+1

               enddo
            enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_exyz(ex,ey,ez,eg,nelx,nely,nelz)
      integer ex,ey,ez,eg

      nelxy = nelx*nely
 
      ez = 1 +  (eg-1)/nelxy
      ey = mod1 (eg,nelxy)
      ey = 1 +  (ey-1)/nelx
      ex = mod1 (eg,nelx)
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_glo_num(glo_num,nx)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      integer*8 glo_num(lx1*ly1*lz1,lelt)

      integer e

      do e=1,nelt
         call outmat_e_i8(glo_num(1,e),e,nx)
      enddo
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_e_i8(gn,e,nx)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      integer*8 gn(lx1,ly1,lz1)

      integer e

      write(6,*)
      write(6,2) e
      write(6,*)

      do k0=3,1,-2

         k1=k0+1
            write(6,*) k0,k1
         do j=nx,1,-1
            write(6,1) ((gn(i,j,k),i=1,nx),k=k0,k1)
         enddo
         write(6,*)

      enddo
    1 format('gn:',4i7,3x,4i7)
    2 format('gn: element: ',i4)
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_glo_num_general(glo_num,nx)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      integer*8 glo_num(1)

      integer e

      io = 0
      do e=1,nelt
         write(6,*)
         write(6,2) e
         write(6,*)

         do k=nx,1,-1
              write(6,*) 'k = ',k
            do j=nx,1,-1
                write(6,1) (glo_num(igo),igo=1+io,io+nx)
                io = io+nx
            enddo
         enddo

      enddo
    1 format('gn:',6i10)
    2 format('gn: element: ',i4)
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_r(x,name5)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      character*5 name5

      real x(lx1*ly1*lz1,lelt)

      integer e

      do e=1,nelt
         call outmat_e_r(x(1,e),name5,e)
      enddo
 
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat_e_r(x,name5,e)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      character*5 name5

      real x(lx1,ly1,lz1)

      integer e

      write(6,*)
      write(6,2) e,name5
      write(6,*)

      do k0=3,1,-2

         k1=k0+1
         do j=ny1,1,-1
            write(6,1) ((x(i,j,k),i=1,4),k=k0,k1)
         enddo
         write(6,*)

      enddo
    1 format('mat: ',4f8.3,3x,4f8.3)
    2 format('mat: element: ',i4,2x,a5)
 
      return
      end
c-----------------------------------------------------------------------
