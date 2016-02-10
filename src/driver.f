c-----------------------------------------------------------------------
      program nekbone
      
      include 'SIZE'
      include 'TOTAL'
      include 'SEMHAT'
      include 'mpif.h'
c      include 'comm_mpi'
c      include 'omp_lib.h'
      
      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)
      parameter (lxyz = lx1*ly1*lz1)
      parameter (lt=lxyz*lelt)

      real x(lt),f(lt),r(lt),w(lt),p(lt),z(lt),c(lt)
      real g(6,lt)
      real mfloplist(1024), avmflop
      integer icount  

      logical ifbrick
      integer iel0,ielN,ielD   ! element range per proc.
      integer nx0,nxN,nxD      ! poly. order range
      integer npx,npy,npz      ! processor decomp
      integer mx ,my ,mz       ! element decomp

c      write to file                                                             
      character(len=1024) :: filename

c     topology                                                                  
      integer comm_top,ndims2,ierr,source,dest

      
      integer dim_size(3)
      logical periods(3)
      integer coords(3)
      logical reorder
      integer direction, displ
      integer dir_x, dir_y,dir_z

      integer x_src,x_dest,y_src,y_dest,z_src,z_dest

c     topology                                                                
c     local to global                                                          
      integer*8 glo_num(lx1*ly1*lz1*lelt)

      integer n_shared
      integer send_buf(2)
      integer recv_buf(2)
      integer status99
c     local to global                             


c     testing
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      integer,dimension(8) :: values
! using keyword arguments
      call date_and_time(date,time,zone,values)
      write(*,*) 'time', time
c     testing



c     set topology                                                                       
      direction =  1   ! shift along the 1st index (0 or 1)                              

c      dir_x = 1
c      dir_y = 0 
c      dir_z = 2

      dir_x = 2
      dir_y = 1
      dir_z = 0
c      dir_x = 0
c      dir_y = 1
c      dir_z = 2

      displ =  1                ! shift by  1                                            
      ndims = 3            ! 2D  grid                                                    
      periods(1) = .false.
      periods(2) = .false.
      periods(3) = .false.
      reorder = .false.
c     topology        



      call iniproc(mpi_comm_world)    ! has nekmpi common block
      call init_delay

      call read_param(ifbrick,iel0,ielN,ielD,nx0,nxN,nxD,
     $                               npx,npy,npz,mx,my,mz)

c!$OMP PARALLEL                                              
c      write(*,'(A,I)') 'numthreads= ', omp_get_num_threads()
c!$OMP END PARALLEL    

c     npy
c      dim_size(1) = npy           ! rows 
c     npx
c      dim_size(2) = npx
c     npz
c      dim_size(3) = npz

      dim_size(1) = npz 
      dim_size(2) = npy
      dim_size(3) = npx

c     topology                                                                           




      call MPI_Cart_create(mpi_comm_world, ndims, dim_size,
     &     periods, reorder, comm_top, ierr)
      call MPI_Cart_coords(comm_top, nid, ndims, coords, ierr)

      call MPI_Cart_shift(comm_top,dir_x,displ,x_src,x_dest,ierr)
      call MPI_Cart_shift(comm_top,dir_y,displ,y_src,y_dest,ierr)
      call MPI_Cart_shift(comm_top,dir_z,displ,z_src,z_dest,ierr)


c      write(filename,'(A3,I2,A4)') 'top',nid,'.txt'
c      write(*,'(A)') filename

c     GET PLATFORM CHARACTERISTICS
c     iverbose = 1
c     call platform_timer(iverbose)   ! iverbose=0 or 1
c      write(*,*) 'nid = ', nid
      icount = 0
      

c     SET UP and RUN NEKBONE
      do nx1=nx0,nxN,nxD
         call init_dim
         do nelt=iel0,ielN,ielD
           call init_mesh(ifbrick,cmask,npx,npy,npz,mx,my,mz)
c           write(*,'(A,I2,A,I2,A,I2)') 'npx',npx, 'npy',npy, 'npz',npz
c           write(*,'(A,I2,A,I2,A,I2)') 'mx',mx,'my',my,'mz',mz
           call proxy_setupds(gsh,nx1,glo_num,n_shared)
           call faces_ids(nx1,mx,my,mz)
c           call set_multiplicity (c)       ! Inverse of counting matrix

           call fset_multiplicity(c,n_shared,
     $          x_src,x_dest,y_src,y_dest,
     $          z_src,z_dest,
     $          npx,npy,npz,mx,my,mz)

           call proxy_setup(ah,bh,ch,dh,zh,wh,g) 
           call h1mg_setup

           niter = 100
           n     = nx1*ny1*nz1*nelt
c           call rone(f,n)
c           write(*,'(A,I2)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!n= ', n
c           call set_f(f,c,n)
c     fortran set_f
           
           call f_set_f(f, c, n,n_shared,
     $          x_src,x_dest,y_src,y_dest,
     $          z_src,z_dest,
     $          npx,npy,npz,mx,my,mz)

c           if(nid.eq.0) then
c              do i=1,n
c                 write(*,*) 'f ',i,' = ',f(i)
c              enddo
c           endif

c           if(nid.eq.0) then
c              do i=1,n
c                 write(*,*) 'f ',i,' ',f(i)
c              enddo
c           endif

           
           if(nid.eq.0) write(6,*)
c           write(*,*) 'testing'
           call cg(x,f,g,c,r,w,p,z,n,niter,flop_cg,n_shared,
     $          x_src,x_dest,y_src,y_dest,
     $          z_src,z_dest,
     $          npx,npy,npz,mx,my,mz)

           call nekgsync()

           call set_timer_flop_cnt(0)

           call cg(x,f,g,c,r,w,p,z,n,niter,flop_cg,n_shared,
     $          x_src,x_dest,y_src,y_dest,
     $          z_src,z_dest,
     $          npx,npy,npz,mx,my,mz)

           call set_timer_flop_cnt(1)
c           open(unit=nid,file=filename,status='replace')
c           do i=1,n
c              write(nid,'(A,I5,A,1pe12.4)') 'x[ ',i, ' ]= ',x(i)
c           enddo
c           close(unit=nid)

c           call gs_free(gsh)
           
           icount = icount + 1
           mfloplist(icount) = mflops*np
         enddo
      enddo

      avmflop = 0.0
      do i = 1,icount
         avmflop = avmflop+mfloplist(i)
      enddo
      if(icount.ne.0) then
         avmflop=avmflop/icount
      endif
      if(nid.eq.0) then
         write(6,1) avmflop
      endif
  1   format('Avg MFlops = ',1pe12.4)

c     TEST BANDWIDTH BISECTION CAPACITY
c     call xfer(np,cr_h)

      call exitt0

      end
c--------------------------------------------------------------
      subroutine set_f(f,c,n)
      real f(n),c(n)

      do i=1,n
         arg  = 1.e9*(i*i)
         arg  = 1.e9*cos(arg)
         f(i) = sin(arg)

      enddo

      call dssum(f)

      call col2 (f,c,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine f_set_f(f, c, n,n_shared,
     $     x_src,x_dest,y_src,y_dest,
     $     z_src,z_dest,
     $     npx,npy,npz,mx,my,mz)
      include 'mpif.h'
      include 'SIZE'
      include 'TOTAL'

      real f(n), c(n)
    
      integer n_shared
      
      integer x_src,x_dest,y_src,y_dest,z_src,z_dest
      integer npx,npy,npz,mx,my,mz

      do i=1,n
         arg  = 1.e9*(i*i)
         arg  = 1.e9*cos(arg)
         f(i) = sin(arg)
c         f(i) = 1.
      enddo
      call mytimer(0)
      call sync_xyz(f,n_shared,
     $     x_src,x_dest,y_src,y_dest,
     $     z_src,z_dest,
     $     npx,npy,npz,mx,my,mz)
      call mytimer(1)
      call col2 (f,c,n)

      return 
      end
c-----------------------------------------------------------------------

      subroutine init_dim

C     Transfer array dimensions to common

      include 'SIZE'
 
      ny1=nx1
      nz1=nx1
 
      ndim=ldim

      return
      end
c-----------------------------------------------------------------------
      subroutine init_mesh(ifbrick,cmask,npx,npy,npz,mx,my,mz)
      include 'SIZE'
      include 'TOTAL'
      real cmask(-1:nx1*ny1*nz1*nelt)
      logical ifbrick
      integer e,eg,offs,npx,npy,npz,mx,my,mz
 
c     Trigger reset of mask
      cmask(-1) = 1.0

      if(.not.ifbrick) then   ! A 1-D array of elements of length P*lelt
         nelx = nelt*np
         nely = 1
         nelz = 1

         npx=np
         npy=1
         npz=1

         mx=nelt
         my=1
         mz=1

         if(nid.eq.0) then 
           write(6,*)
           write(6,*) 'Processor Distribution:  npx,npy,npz=',npx,
     $                npy,npz
           write(6,*) 'Element Distribution: nelx,nely,nelz=',nelx,
     $                nely,nelz
           write(6,*) 'Local Element Distribution: mx,my,mz=',mx,
     $                my,mz
         endif
   
         do e=1,nelt
            eg = e + nid*nelt
            lglel(e) = eg
c            write(*,'(A,I,A,I)') 'nid= ', nid, ' eg = ', eg 
         enddo
      else              ! A 3-D block of elements 
         !xyz distribution of total proc if user-provided isn't valid
         if(npx*npy*npz.ne.np) then
            call cubic(npx,npy,npz,np) 
         endif

         !xyz distribution of total NELT if user-provided isn't valid
         if(mx*my*mz.ne.nelt)  then
            call cubic(mx,my,mz,nelt) 
         endif
      
         nelx = mx*npx
         nely = my*npy 
         nelz = mz*npz

         if(nid.eq.0) then 
           write(6,*)
           write(6,*) 'Processor Distribution:  npx,npy,npz=',npx,
     $                npy,npz
           write(6,*) 'Element Distribution: nelx,nely,nelz=',nelx,
     $                nely,nelz
           write(6,*) 'Local Element Distribution: mx,my,mz=',mx,
     $                my,mz
         endif

         e = 1
         offs = (mod(nid,npx)*mx) + npx*(my*mx)*(mod(nid/npx,npy)) 
     $       + (npx*npy)*(mx*my*mz)*(nid/(npx*npy))
         do k = 0,mz-1
         do j = 0,my-1
         do i = 0,mx-1
            eg = offs+i+(j*nelx)+(k*nelx*nely)+1
            lglel(e) = eg
c            write(*,'(A,I1,A,I3,A,I6)') 'id=',nid, 'e=', e, 'lglel=',eg
            e        = e+1
            
         enddo
         enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cubic(mx,my,mz,np)

      mx = np
      my = 1
      mz = 1
      ratio = np

      iroot3 = np**(1./3.) + 0.000001
      do i = iroot3,1,-1
        iz = i
        myx = np/iz
        nrem = np-myx*iz

        if (nrem.eq.0) then
          iroot2 = myx**(1./2.) + 0.000001
          do j=iroot2,1,-1
            iy = j
            ix = myx/iy
            nrem = myx-ix*iy
            if (nrem.eq.0) goto 20
          enddo
   20     continue

          if (ix.lt.iy) then
            it = ix
            ix = iy
            iy = it
          endif

          if (ix.lt.iz) then
            it = ix
            ix = iz
            iz = it
          endif

          if (iy.lt.iz) then
            it = iy
            iy = iz
            iz = it
          endif

          if ( REAL(ix)/iz.lt.ratio) then
            ratio = REAL(ix)/iz
            mx = ix
            my = iy
            mz = iz
          endif
        endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cubic0(mx,my,mz,np)

        rp = np**(1./3.)
        mz = rp*(1.01)

        do iz=mz,1,-1
           myx = np/iz
           nrem = np-myx*iz
           if (nrem.eq.0) goto 10
        enddo
   10   mz = iz
        rq = myx**(1./2.)
        my = rq*(1.01)
        do iy=my,1,-1
           mx = myx/iy
           nrem = myx-mx*iy
           if (nrem.eq.0) goto 20
        enddo
   20   my = iy

        mx = np/(mz*my)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_multiplicity (c)       ! Inverse of counting matrix
      include 'SIZE'
      include 'TOTAL'

      real c(1)

      n = nx1*ny1*nz1*nelt

      call rone(c,n)
      call adelay
c      call gs_op(gsh,c,1,1,0)  ! Gather-scatter operation  ! w   = QQ  w

      do i=1,n
         c(i) = 1./c(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine fset_multiplicity(c,n_shared,
     $     x_src,x_dest,y_src,y_dest,
     $     z_src,z_dest,
     $     npx,npy,npz,mx,my,mz)
      include 'mpif.h'
      include 'SIZE'
      include 'TOTAL'

      real c(1)

      integer n_shared

      integer x_src,x_dest,y_src,y_dest,z_src,z_dest
      integer npx,npy,npz,mx,my,mz


      n = nx1*ny1*nz1*nelt
      call rone(c,n)
      call adelay

      call mytimer(0)
      call sync_xyz(c,n_shared,
     $     x_src,x_dest,y_src,y_dest,
     $     z_src,z_dest,
     $     npx,npy,npz,mx,my,mz)
      call mytimer(1)

c      call gs_op(gsh,c,1,1,0)  ! Gather-scatter operation  ! w   = QQ  w                                   

      do i=1,n
         c(i) = 1./c(i)
      enddo



      return
      end
c----------------------------------------------------------------------- 
      subroutine set_timer_flop_cnt(iset)
      include 'SIZE'
      include 'TOTAL'

      real time0,time1
      save time0,time1

      if (iset.eq.0) then
         flop_a  = 0
         flop_cg = 0
         time0   = dnekclock()
      else
        time1   = dnekclock()-time0
        if (time1.gt.0) mflops = (flop_a+flop_cg)/(1.e6*time1)
c        write(6,2) mflops*np, mflops
c        write(*,*) 'nid ', nid, ' solve time ', time1
c        write(6,4) time1
        
        if (nid.eq.0) then
          write(6,*)
          write(6,1) nelt,np,nx1, nelt*np
          write(6,2) mflops*np, mflops
          write(6,3) flop_a,flop_cg
          write(6,4) time1
        endif
    1   format('nelt = ',i7, ', np = ', i9,', nx1 = ', i7,
     &         ', elements =', i10 )
    2   format('Tot MFlops = ', 1pe12.4, ', MFlops      = ', e12.4)
    3   format('Setup Flop = ', 1pe12.4, ', Solver Flop = ', e12.4)
    4   format('Solve Time = ', e12.4)
      endif

      return
      end
c-----------------------------------------------------------------------
c     conflict if only fortran
c      subroutine xfer(np,gsh)
c      include 'SIZE'
c      parameter(npts_max = lx1*ly1*lz1*lelt)

c      real buffer(2,npts_max)
c      integer ikey(npts_max)


c      nbuf = 800
c      npts = 1
c      do itest=1,200
c         npoints = npts*np

c         call load_points(buffer,nppp,npoints,npts,nbuf)
c         iend   = mod1(npoints,nbuf)
c         istart = 1
c         if(nid.ne.0)istart = iend+(nid-1)*nbuf+1
c         do i = 1,nppp
c            icount=istart+(i-1)
c            ikey(i)=mod(icount,np)
c         enddo

c         call nekgsync
c         time0 = dnekclock()
c         do loop=1,50
c            call crystal_tuple_transfer(gsh,nppp,npts_max,
c     $                ikey,1,ifake,0,buffer,2,1)
c         enddo
c         time1 = dnekclock()
c         etime = (time1-time0)/50
c
c         if (nid.eq.0) write(6,1) np,npts,npoints,etime
c   1     format(2i7,i10,1p1e12.4,' bandwidth' )
c         npts = 1.02*(npts+1)
c         if (npts.gt.npts_max) goto 100
c      enddo
c 100  continue
c
c      return
c      end
c-----------------------------------------------------------------------
      subroutine load_points(buffer,nppp,npoints,npts,nbuf)
      include 'SIZE'
      include 'PARALLEL'

      real buffer(2,nbuf)

      nppp=0
      if(nbuf.gt.npts) then
       npass = 1+npoints/nbuf

       do ipass = 1,npass
          if(nid.eq.ipass.and.ipass.ne.npass) then
            do i = 1,nbuf
             buffer(1,i)=i
             buffer(2,i)=nid
            enddo
            nppp=nbuf
          elseif (npass.eq.ipass.and.nid.eq.0) then
            mbuf=mod1(npoints,nbuf)
            do i=1,mbuf
               buffer(1,i)=i
               buffer(2,i)=nid
            enddo
            nppp=mbuf
          endif
       enddo
      else
       do i = 1,npts
          buffer(1,i)=i
          buffer(2,i)=nid
       enddo
       nppp=npts
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine read_param(ifbrick,iel0,ielN,ielD,nx0,nxN,nxD,
     $                                     npx,npy,npz,mx,my,mz)
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      logical ifbrick
      integer iel0,ielN,ielD,nx0,nxN,nxD,npx,npy,npz,mx,my,mz

      !open .rea
      ifbrick = .false.  
      ifmgrid = .false.  !initialize to false
      npx=0
      npy=0
      npz=0
      mx =0
      my =0
      mz =0

      if(nid.eq.0) then
         open(unit=9,file='data.rea',status='old') 
         read(9,*,err=100) ifbrick
         read(9,*,err=100) iel0,ielN,ielD
         read(9,*,err=100) nx0,nxN,nxD
         read(9,*,err=100,iostat=ii) npx,npy,npz !optional
         read(9,*,err=100,iostat=ii) mx,my,mz    !optional 
         close(9)
      endif
      call bcast(ifbrick,4)

      call bcast(iel0,4)! ELEMENT range
      call bcast(ielN,4)
      call bcast(ielD,4)

      call bcast(nx0,4) ! POLY Order range
      call bcast(nxN,4)
      call bcast(nxD,4)

      call bcast(npx,4) ! PROC decomp
      call bcast(npy,4)
      call bcast(npz,4)

      call bcast(mx,4)  ! NELT decomp
      call bcast(my,4)
      call bcast(mz,4)

#ifdef MGRID
      ifmgrid = .true.
#endif

      if(iel0.gt.ielN.or.nx0 .gt.nxN)    goto 200
      if(ielN.gt.lelt.or.nxN .gt.lx1)    goto 210
      if(ielD.gt.ielN.or.nxD.gt.nxN)     goto 220
      if(nx0.lt.4.and.ifmgrid)           goto 230

      if(nid.eq.0) write(6,*) "ifmgrid    :",ifmgrid
     $                   ,"    ifbrick    :",ifbrick

      return

  100 continue
      write(6,*) "ERROR READING    data.rea.....ABORT"
      write(6,*) "CHECK PARAMETERS data.rea.....ABORT"
      call exitt0

  200 continue
      write(6,*) "ERROR data.rea :: iel0 > ielN or nx0 > nxN :: ABORT"
      call exitt0

  210 continue
      write(6,*) "ERROR data.rea : ielN>lelt or nxN>lx1(SIZE) :: ABORT"
      call exitt0
  
  220 continue
      write(6,*) "WARNING data.rea : STRIDE   nxD>nxN or ielD>ielN !!!"
  
  230 continue
      write(6,*) "ERROR data.rea nx0 must be greater than or equal to 4"
     $          ," ...setting nx0=4"
      nx0=4
  
  
      return
      end
c-----------------------------------------------------------------------
