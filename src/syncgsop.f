c-------------------------------------------------------------------------
      subroutine sync_xyz(w,n_shared,
     $     x_src,x_dest,y_src,y_dest,
     $     z_src,z_dest,
     $     npx,npy,npz,mx,my,mz)
      use, intrinsic :: ISO_C_BINDING
      use gaspi
      include 'mpif.h'
      include 'SIZE'
      include 'TOTAL'
      
      real w(nx1*ny1*nz1*nelt)
c     real w(nx1*ny1*nz1,nelt)
      integer n_shared
      integer x_src,x_dest,y_src,y_dest,z_src,z_dest
      integer npx,npy,npz,mx,my,mz

      integer arr_size,cnst
      integer i,j,k

      integer(gaspi_queue_id_t):: queue_faces

      cnst = 16
      queue_faces = 0
      

      arr_size = 4*n_shared*(my*mz+mx*mz+mx*my)+12*cnst
      
      call local_ox(w,mx,n_shared)
      call local_oy(w,mx,my,mz,n_shared)
      call local_oz(w,mx,my,mz,n_shared)

c     sync along X dir
      call gpi_packing_face(w,arr_size,offset_right,n_shared*mz*my,
     $     face_right)
      call gpi_send_face(0,x_dest,8*off_from_left,
     $     8*my*mz*n_shared,1,1,queue_faces)      
c     gpi_send_face(offset_loc,rank_rem,offset_rem,
c     $     size,n_id,n_val,queue_arg)
      call gpi_wait(queue_faces)
      call gpi_packing_face(w,arr_size,offset_left,n_shared*mz*my,
     $     face_left)
      call gpi_send_face(8*offset_left,x_src,8*off_from_right,
     $     8*my*mz*n_shared,2,2,queue_faces)
      call gpi_wait(queue_faces)
      call gpi_waitsome_face(x_src,1)
      call gpi_sum_face(w,arr_size,off_from_left,x_src,
     $     my*mz*n_shared,face_left)
      call gpi_waitsome_face(x_dest,2)
      call gpi_sum_face(w,arr_size,off_from_right,x_dest,
     $     my*mz*n_shared,face_right)
c     sync along X dir


c     sync along Y dir
      call gpi_packing_face(w,arr_size,offset_top,n_shared*mx*mz,
     $     face_y_top)
      call gpi_send_face(8*offset_top,y_dest,8*off_from_bottom,
     $     8*mx*mz*n_shared,3,3,queue_faces)
      call gpi_wait(queue_faces)
      call gpi_packing_face(w,arr_size,offset_bottom,n_shared*mx*mz,
     $     face_y_bottom)
      call gpi_send_face(8*offset_bottom,y_src,8*off_from_top,
     $     8*mx*mz*n_shared,4,4,queue_faces)
      call gpi_wait(queue_faces)
      call gpi_waitsome_face(y_src,3)
      call gpi_sum_face(w,arr_size,off_from_bottom,y_src,
     $     mx*mz*n_shared,face_y_bottom)
      call gpi_waitsome_face(y_dest,4)
      call gpi_sum_face(w,arr_size,off_from_top,y_dest,
     $     mx*mz*n_shared,face_y_top)
c     sync along Y dir


c     sync along Z dir
      call gpi_packing_face(w,arr_size,offset_front,n_shared*mx*my,
     $     face_z_front)
      call gpi_send_face(8*offset_front,z_src,8*off_from_back,
     $     8*mx*my*n_shared,5,5,queue_faces)
      call gpi_wait(queue_faces)
      call gpi_packing_face(w,arr_size,offset_back,n_shared*mx*my,
     $     face_z_back)
      call gpi_send_face(8*offset_back,z_dest,8*off_from_front,
     $     8*mx*my*n_shared,6,6,queue_faces)
      call gpi_wait(queue_faces)
      call gpi_waitsome_face(z_dest,5)
      call gpi_sum_face(w,arr_size,off_from_back,z_dest,
     $     mx*my*n_shared,face_z_back)
      call gpi_waitsome_face(z_src,6)
      call gpi_sum_face(w,arr_size,off_from_front,z_src,
     $     mx*my*n_shared,face_z_front)
c     sync along Z dir
      
      call gpi_wait(queue_faces) 
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
