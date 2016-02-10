c-------------------------------------------------------------------------
      subroutine sync_xyz(w,n_shared,
     $     x_src,x_dest,y_src,y_dest,
     $     z_src,z_dest,
     $     npx,npy,npz,mx,my,mz)
      include 'mpif.h'
      include 'SIZE'
      include 'TOTAL'
      
      real w(nx1*ny1*nz1*nelt)
c     real w(nx1*ny1*nz1,nelt)
      integer n_shared
      integer x_src,x_dest,y_src,y_dest,z_src,z_dest
      integer npx,npy,npz,mx,my,mz

c      real*8 data_send_x_right(n_shared*my*mz)
c      real*8 data_send_x_left(n_shared*my*mz)
c      real*8 msg_x_from_left(n_shared*my*mz)
c      real*8 msg_x_from_right(n_shared*my*mz)

c      real*8 data_send_y_top(n_shared*mx*mz)
c      real*8 data_send_y_bottom(n_shared*mx*mz)
c      real*8 msg_y_from_top(n_shared*mx*mz)
c      real*8 msg_y_from_bottom(n_shared*mx*mz)

c      real*8 data_send_z_front(n_shared*mx*my)
c      real*8 data_send_z_back(n_shared*mx*my)
c      real*8 msg_z_from_front(n_shared*mx*my)
c      real*8 msg_z_from_back(n_shared*mx*my)
      
      
      call local_ox(w,mx,n_shared)
      call local_oy(w,mx,my,mz,n_shared)
      call local_oz(w,mx,my,mz,n_shared)
c      write(*,*) 'after local sync'
c      call packing_ox(w,mx,my,mz,npx,n_shared,
c     $     data_send_x_right,data_send_x_left)
c     packing_face(w,piter,pack_data,face)
      call packing_face(w,n_shared*my*mz,data_send_x_right,
     $     face_right)
      call packing_face(w,n_shared*my*mz,data_send_x_left,
     $     face_left)

c      call sending_ox(npx,my,mz,x_src,x_dest,n_shared,
c     $     data_send_x_right,data_send_x_left,
c     $     msg_x_from_left,msg_x_from_right)
c      write(*,*) 'after sending OX'

c     works but slow problems with flags
      call sending_ox_test1(npx,my,mz,x_src,x_dest,n_shared,
     $     data_send_x_right,data_send_x_left,
     $     msg_x_from_left,msg_x_from_right)

c      call sending_ox_test2(npx,my,mz,x_src,x_dest,n_shared,
c     $     data_send_x_right,data_send_x_left,
c     $     msg_x_from_left,msg_x_from_right)

c      call sum_ox(w,npx,mx,my,mz,n_shared,x_src,x_dest,
c     $     msg_x_from_left,msg_x_from_right)
c     sum_face(w,piter,msg_data,face,rank_rem)
      call sum_face(w,n_shared*my*mz,msg_x_from_left,
     $     face_left,x_src)
      call sum_face(w,n_shared*my*mz,msg_x_from_right,
     $     face_right,x_dest)

c      call packing_oy(w,mx,my,mz,npy,n_shared,
c     $     data_send_y_top,data_send_y_bottom)
c     face_y_top(lx1*lz1*lelt),face_y_bottom(lx1*lz1*lelt)
c     $     ,face_z_front(lx1*ly1*lelt),face_z_back(lx1*ly1*lelt)
      call packing_face(w,n_shared*mx*mz,data_send_y_top,
     $     face_y_top)
      call packing_face(w,n_shared*mx*mz,data_send_y_bottom,
     $     face_y_bottom)
      call sending_oy(npx,npy,mx,mz,y_src,y_dest,n_shared,
     $     data_send_y_top,data_send_y_bottom,
     $     msg_y_from_top,msg_y_from_bottom)
c      call sum_oy(w,npy,mx,my,mz,n_shared,y_src,y_dest,
c     $     msg_y_from_top,msg_y_from_bottom)
      call sum_face(w,n_shared*mx*mz,msg_y_from_top,
     $     face_y_top,y_dest)
      call sum_face(w,n_shared*mx*mz,msg_y_from_bottom,
     $     face_y_bottom,y_src)
      
c      call packing_oz(w,mx,my,mz,npz,n_shared,
c     $     data_send_z_front,data_send_z_back)
      call packing_face(w,n_shared*mx*my,data_send_z_front,
     $     face_z_front)
      call packing_face(w,n_shared*mx*my,data_send_z_back,
     $     face_z_back)
c      call exit(0)
      call sending_oz(npx,npy,npz,mx,my,z_src,z_dest,n_shared,
     $     data_send_z_front,data_send_z_back,
     $     msg_z_from_front,msg_z_from_back)
c      write(*,*) 'after sending OZ'

c      call check_arr(msg_z_from_back,n_shared*mx*my)
c      call check_arr(msg_z_from_front,n_shared*mx*my)
c      call exit(0)
c     write(*,*) 
c      call exit(0)
c      call sum_oz(w,npz,n_shared*mx*my,z_src,z_dest,
c     $     msg_z_from_front,msg_z_from_back)
      call sum_face(w,n_shared*mx*my,msg_z_from_front,
     $     face_z_front,z_src)
      call sum_face(w,n_shared*mx*my,msg_z_from_back,
     $     face_z_back,z_dest)
      
      return 
      end
c------------------------------------------------------------------------- 
      subroutine check_arr(array,iter)
      include 'SIZE'
      include 'TOTAL'

      real*8 array(1)
      integer iter
      
      integer i
      if(nid.eq.20) then
         do i=1,iter
            write(*,*) 'array i ', i, ' = ',array(i)
         enddo
      endif

      return
      end
c------------------------------------------------------------------------- 
      subroutine sending_ox_test1(npx,my,mz,x_src,x_dest,n_shared,
     $     d_send_x_right,d_send_x_left,
     $     m_x_from_left,m_x_from_right)
      include 'SIZE'
      include 'mpif.h'
      include 'TOTAL'
      
      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
      integer npx
      integer  my, mz
      integer x_src,x_dest
      integer n_shared
c      real*8 d_send_x_right(1)
c      real*8 d_send_x_left(1)
c      real*8 m_x_from_left(1)
c      real*8 m_x_from_right(1)

      real*8 d_send_x_right(n_shared*my*mz)
      real*8 d_send_x_left(n_shared*my*mz)
      real*8 m_x_from_left(n_shared*my*mz)
      real*8 m_x_from_right(n_shared*my*mz)
      
      if(mod(nid,2).eq.0) then
         if(x_dest.ge.0) then
            call mpi_send(d_send_x_right,8*n_shared*my*mz,
     $           mpi_byte,
     $           x_dest,1,nekcomm,ierr)

            call mpi_recv(m_x_from_right,8*n_shared*my*mz,
     $           mpi_byte,
     $           x_dest,0,nekcomm,status,ierr)
c            write(*,*) 'nid ',nid,' sent right to ', x_dest 
         endif
      endif
      
      if(mod(nid,2).eq.1) then
         if(x_src.ge.0) then
            call mpi_recv(m_x_from_left,8*n_shared*my*mz,
     $           mpi_byte,
     $           x_src,1,nekcomm,status,ierr)
c            write(*,*) 'nid ',nid,' sent left to ', x_src                                      
            call mpi_send(d_send_x_left,8*n_shared*my*mz,
     $           mpi_byte,
     $           x_src,0,nekcomm,ierr)
         endif
      endif

c      call mpi_barrier(nekcomm,ierr)

      if(mod(nid,2).eq.0) then   
         if(x_src.ge.0) then                                                               
            call mpi_send(d_send_x_left,8*n_shared*my*mz,                     
     $           mpi_byte,                                                      
     $           x_src,0,nekcomm,ierr)                                    
            call mpi_recv(m_x_from_left,8*n_shared*my*mz,                 
     $           mpi_byte, 
     $           x_src,1,nekcomm,status,ierr)              
         endif
      endif
      if(mod(nid,2).eq.1) then
         if(x_dest.ge.0) then
            call mpi_recv(m_x_from_right,8*n_shared*my*mz,
     $           mpi_byte,
     $           x_dest,0,nekcomm,status,ierr)
c            write(*,*) 'nid ',nid,' sent left to ', x_src        
            call mpi_send(d_send_x_right,8*n_shared*my*mz,
     $           mpi_byte,
     $           x_dest,1,nekcomm,ierr)
         endif
      endif
      return
      end
c-------------------------------------------------------------------------                            
      subroutine sending_ox_test2(npx,my,mz,x_src,x_dest,n_shared,
     $     d_send_x_right,d_send_x_left,
     $     m_x_from_left,m_x_from_right)
      include 'SIZE'
      include 'mpif.h'
      include 'TOTAL'

      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
      integer npx
      integer  my, mz
      integer x_src,x_dest
      integer n_shared
      real*8 d_send_x_right(n_shared,my,mz)
      real*8 d_send_x_left(n_shared,my,mz)
      real*8 m_x_from_left(n_shared,my,mz)
      real*8 m_x_from_right(n_shared,my,mz)

      if(mod(nid,2).eq.1) then
         if(x_dest.ge.0) then
            call mpi_send(d_send_x_right,8*n_shared*my*mz,
     $           mpi_byte,
     $           x_dest,1,nekcomm,ierr)
c     write(*,*) 'nid ',nid,' sent right to ', x_dest                                                 
            call mpi_recv(m_x_from_right,8*n_shared*my*mz,
     $           mpi_byte,
     $           x_dest,0,nekcomm,status,ierr)
         endif
      else if(x_src.ge.0) then
         call mpi_recv(m_x_from_left,8*n_shared*my*mz,
     $        mpi_byte,
     $        x_src,1,nekcomm,status,ierr)
c     write(*,*) 'nid ',nid,' sent left to ', x_src                                                   
         call mpi_send(d_send_x_left,8*n_shared*my*mz,
     $        mpi_byte,
     $        x_src,0,nekcomm,ierr)
      endif
      

      return
      end


c-------------------------------------------------------------------------   
      subroutine sending_ox(npx,my,mz,x_src,x_dest,n_shared,
     $     d_send_x_right,d_send_x_left,
     $     m_x_from_left,m_x_from_right)
c     send and recieve data along OX                                                  
      include 'SIZE'
      include 'mpif.h'
      include 'TOTAL'
      
      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
      integer npx
      integer  my, mz
      integer x_src,x_dest
      integer n_shared
      real*8 d_send_x_right(n_shared,my,mz)
      real*8 d_send_x_left(n_shared,my,mz)
      real*8 m_x_from_left(n_shared,my,mz)
      real*8 m_x_from_right(n_shared,my,mz)

      integer proc_step_x
      if(npx.gt.1)then
         do proc_step_x = 0,npx,2
            if(mod(nid,npx).eq.proc_step_x) then
               if(x_dest.ge.0) then
                  call mpi_send(d_send_x_right,8*n_shared*my*mz,
     $                 mpi_byte,
     $                 x_dest,1,nekcomm,ierr)
c                  write(*,*) 'nid ',nid,' sent right to ', x_dest
                  call mpi_recv(m_x_from_right,8*n_shared*my*mz,
     $                 mpi_byte,
     $                 x_dest,0,nekcomm,status,ierr)
               endif
            else if(mod(nid,npx).eq.(proc_step_x+1) ) then
               if(x_src.ge.0) then
                  call mpi_recv(m_x_from_left,8*n_shared*my*mz,
     $                 mpi_byte,
     $                 x_src,1,nekcomm,status,ierr)
c                  write(*,*) 'nid ',nid,' sent left to ', x_src
                  call mpi_send(d_send_x_left,8*n_shared*my*mz,
     $                 mpi_byte,
     $                 x_src,0,nekcomm,ierr)
               endif
            endif

            if(mod(nid,npx).eq.(proc_step_x+1)) then
               if(x_dest.ge.0) then
                  call mpi_send(d_send_x_right,8*n_shared*my*mz,
     $                 mpi_byte,
     $                 x_dest,1,nekcomm,ierr)
c     write(*,*) 'nid ',nid,' sent right to ', x_dest                                    
                  call mpi_recv(m_x_from_right,8*n_shared*my*mz,
     $                 mpi_byte,
     $                 x_dest,0,nekcomm,status,ierr)
               endif
            else if(mod(nid,npx).eq.(proc_step_x)) then
               if(x_src.ge.0) then
                  call mpi_recv(m_x_from_left,8*n_shared*my*mz,
     $                 mpi_byte,
     $                 x_src,1,nekcomm,status,ierr)
                  call mpi_send(d_send_x_left,8*n_shared*my*mz,
     $                 mpi_byte,
     $                 x_src,0,nekcomm,ierr)
c     write(*,*) 'nid ',nid,' sent left to ', x_src                                      
               endif
            endif
           
         enddo
      endif
c      if(npx.gt.1)then
c         do proc_step_x = 0,npx,2
c            if(mod(nid,npx).eq.(proc_step_x+1)) then
c               if(x_dest.ge.0) then
c                  call mpi_send(d_send_x_right,8*n_shared*my*mz,
c     $                 mpi_byte,
c     $                 x_dest,1,nekcomm,ierr)
c                  write(*,*) 'nid ',nid,' sent right to ', x_dest
c                  call mpi_recv(m_x_from_right,8*n_shared*my*mz,
c     $                 mpi_byte,
c     $                 x_dest,0,nekcomm,status,ierr)
c               endif
c            else if(mod(nid,npx).eq.(proc_step_x)) then
c               if(x_src.ge.0) then
c                  call mpi_recv(m_x_from_left,8*n_shared*my*mz,
c     $                 mpi_byte,
c     $                 x_src,1,nekcomm,status,ierr)
c                  call mpi_send(d_send_x_left,8*n_shared*my*mz,
c     $                 mpi_byte,
c     $                 x_src,0,nekcomm,ierr)
c                  write(*,*) 'nid ',nid,' sent left to ', x_src
c               endif
c            endif
c         enddo
c      endif
c      write(*,*) 'after sending ox nid ', nid
c     send and recieve data along OX          
      return
      end
c-------------------------------------------------------------------------  
      subroutine sending_oy(npx,npy,mx,mz,y_src,y_dest,n_shared,
     $     d_send_y_top,d_send_y_bottom,
     $     m_y_from_top,m_y_from_bottom)
c     send msg along OY                                                              
      include 'SIZE'
      include 'mpif.h'
      include 'TOTAL'

      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
      integer npx,npy
      integer  mx, mz
      integer y_src,y_dest
      integer n_shared
      
      real*8 d_send_y_top(1)
      real*8 d_send_y_bottom(1)
      real*8 m_y_from_top(1)
      real*8 m_y_from_bottom(1)

c      real*8 d_send_y_top(n_shared*mx*mz)
c      real*8 d_send_y_bottom(n_shared*mx*mz)
c      real*8 m_y_from_top(n_shared*mx*mz)
c      real*8 m_y_from_bottom(n_shared*mx*mz)

      if(npy.gt.1)then
         if(mod(nid,2*npx).lt.npx) then
            if(y_dest.ge.0) then
               call mpi_send(d_send_y_top,8*n_shared*mx*mz,mpi_byte,
     $              y_dest,1,nekcomm,ierr)
               call mpi_recv(m_y_from_top,8*n_shared*mx*mz,mpi_byte,
     $              y_dest,0,nekcomm,status,ierr)
            endif
         else
            if(y_src.ge.0) then
               call mpi_recv(m_y_from_bottom,8*n_shared*mx*mz,
     $              mpi_byte,
     $              y_src,1,nekcomm,status,ierr)
               call mpi_send(d_send_y_bottom,8*n_shared*mx*mz,
     $              mpi_byte,
     $              y_src,0,nekcomm,ierr)
            endif
         endif
      endif
c      call mpi_barrier(nekcomm,ierr)
      if(npy.gt.1)then
         if(mod(nid,2*npx).ge.npx) then
            if(y_dest.ge.0) then
               call mpi_send(d_send_y_top,8*n_shared*mx*mz,mpi_byte,
     $              y_dest,1,nekcomm,ierr)
               call mpi_recv(m_y_from_top,8*n_shared*mx*mz,mpi_byte,
     $              y_dest,0,nekcomm,status,ierr)
            endif
         else
            if(y_src.ge.0) then
               call mpi_recv(m_y_from_bottom,8*n_shared*mx*mz,
     $              mpi_byte,
     $              y_src,1,nekcomm,status,ierr)
               call mpi_send(d_send_y_bottom,8*n_shared*mx*mz,
     $              mpi_byte,
     $              y_src,0,nekcomm,ierr)
            endif
         endif
      endif
c     send msg along OY                      
c      write(*,*) 'after sending oy nid ', nid
      return 
      end
c-------------------------------------------------------------------------  
      subroutine sending_oz(npx,npy,npz,mx,my,z_src,z_dest,n_shared,
     $     d_send_z_front,d_send_z_back,
     $     m_z_from_front,m_z_from_back)
c     send msg along OZ                                             
      include 'SIZE'
      include 'mpif.h'
      include 'TOTAL'

      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
      integer npx,npy,npz
      integer  mx, my
      integer z_src,z_dest
      integer n_shared

      real*8 d_send_z_front(1)
      real*8 d_send_z_back(1)
      real*8 m_z_from_front(1)
      real*8 m_z_from_back(1)

c      real*8 d_send_z_front(n_shared*mx*my)
c      real*8 d_send_z_back(n_shared*mx*my)
c      real*8 m_z_from_front(n_shared*mx*my)
c      real*8 m_z_from_back(n_shared*mx*my)

      if(npz.gt.1) then
         if(mod(nid,2*npx*npy).lt.(npx*npy)) then
            if(z_dest.ge.0) then
               call mpi_send(d_send_z_back,8*n_shared*mx*my,mpi_byte,
     $              z_dest,1,nekcomm,ierr)
               call mpi_recv(m_z_from_back,8*n_shared*mx*my,mpi_byte,
     $              z_dest,0,nekcomm,status,ierr)
            endif
         else
            if(z_src.ge.0) then
               call mpi_recv(m_z_from_front,8*n_shared*mx*my,mpi_byte,
     $              z_src,1,nekcomm,status,ierr)
               call mpi_send(d_send_z_front,8*n_shared*mx*my,
     $              mpi_byte,
     $              z_src,0,nekcomm,ierr)
            endif
         endif
      endif
c      call mpi_barrier(nekcomm,ierr)
      if(npz.gt.1) then
         if(mod(nid,2*npx*npy).ge.(npx*npy)) then
            if(z_dest.ge.0) then
               call mpi_send(d_send_z_back,8*n_shared*mx*my,mpi_byte,
     $              z_dest,1,nekcomm,ierr)
               call mpi_recv(m_z_from_back,8*n_shared*mx*my,mpi_byte,
     $              z_dest,0,nekcomm,status,ierr)
            endif
         else
            if(z_src.ge.0) then
               call mpi_recv(m_z_from_front,8*n_shared*mx*my,mpi_byte,
     $              z_src,1,nekcomm,status,ierr)
               call mpi_send(d_send_z_front,8*n_shared*mx*my,
     $              mpi_byte,
     $              z_src,0,nekcomm,ierr)
            endif
         endif
      endif
c      write(*,*) 'after sending oz nid ', nid
c     send msg along OZ                                                 
      return
      end
c------------------------------------------------------------------------- 
