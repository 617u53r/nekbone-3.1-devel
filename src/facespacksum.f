      subroutine gpi_packing_face(w,arr_s,offset,pack_niter,face)
      use, intrinsic :: ISO_C_BINDING
      use gaspi
      include 'SIZE'
      
      real w(1)
      integer pack_niter,arr_s,offset
      integer*8 face(1)
      
      integer i
      real*8, pointer :: arr(:)
      integer(gaspi_return_t) :: ret
      integer(gaspi_size_t) :: seg_size,arr_size
      integer(gaspi_alloc_t) :: seg_alloc
      type(c_ptr) :: seg_ptr
      
      arr_size=arr_s
      
      ret = gaspi_segment_ptr(INT(0,1), seg_ptr)
      if(ret .ne. GASPI_SUCCESS) then
         write(*,*) "gaspi_segment_ptr failed"
         call exit(-1)
      endif
      call c_f_pointer(seg_ptr, arr, shape=[arr_size])
      
      do i=1,pack_niter
         arr(i+offset) = w(face(i))
      enddo

      return
      end
c-------------------------------------------------------------------------      
      subroutine gpi_sum_face(w,arr_s,offset,rank_remote,
     $     pack_niter,face)
      use, intrinsic :: ISO_C_BINDING
      use gaspi
      include 'SIZE'

      integer npdir
      integer*8 face(1)
      real w(nx1*ny1*nz1*nelt)
      integer offset
      integer pack_niter
      integer arr_s
      integer rank_remote

      real*8, pointer :: arr(:)
      integer(gaspi_return_t) :: ret
      integer(gaspi_size_t) :: seg_size,arr_size
      integer(gaspi_alloc_t) :: seg_alloc
      type(c_ptr) :: seg_ptr

      arr_size=arr_s

      ret = gaspi_segment_ptr(INT(0,1), seg_ptr)
      if(ret .ne. GASPI_SUCCESS) then
         write(*,*) "gaspi_segment_ptr failed"
         call exit(-1)
      endif
      call c_f_pointer(seg_ptr, arr, shape=[arr_size])

      if(rank_remote.ge.0) then
         do i=1,pack_niter
            w(face(i)) = w(face(i)) + arr(i+offset)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------
      subroutine gpi_send_face(offset_loc,rank_rem,offset_rem,
     $     size,n_id,n_val,queue_arg)
      use, intrinsic :: ISO_C_BINDING
      use gaspi
      include 'SIZE'
      include 'TOTAL'
      integer offset_loc,rank_rem,offset_rem,size
      integer     n_id,n_val,queue_arg

      integer(gaspi_return_t) :: ret
      integer(gaspi_size_t) :: msg_size
      integer(gaspi_offset_t) :: offset_local,offset_remote
      integer(gaspi_rank_t):: rank_remote
      integer(gaspi_notification_id_t):: notif_id
      integer(gaspi_notification_t)::notif_val
      integer(gaspi_queue_id_t):: queue

      offset_local = offset_loc
      offset_remote = offset_rem
      rank_remote = rank_rem
      msg_size = size
      notif_id = n_id
      notif_val = n_val
      queue=queue_arg

      if(rank_rem.ge.0) then
         ret = gaspi_write_notify(INT(0,1),offset_local,rank_remote,
     $        INT(0,1),offset_remote,msg_size,notif_id,
     $        notif_val,queue, 10000)
         if(ret.eq.GASPI_TIMEOUT) then
            write(*,*) 'gaspi_write_notify timeout'
            call exit(-1)
         endif
         if(ret .ne. GASPI_SUCCESS) then
            write(*,*) "gaspi_barrier failed"
            call exit(-1)
         endif
      endif

      return
      end
c-------------------------------------------------------------------------      
      subroutine gpi_wait(queue_arg)
      use, intrinsic :: ISO_C_BINDING
      use gaspi

      include 'SIZE'
      include 'TOTAL'
      integer queue_arg
      integer(gaspi_queue_id_t):: queue
      integer(gaspi_return_t) :: ret
      queue = queue_arg

      ret = gaspi_wait(queue,10000)
      if(ret.eq.GASPI_TIMEOUT) then
         write(*,*) 'gaspi_wait timeout'
         call exit(-1)
      endif

      return
      end
c-------------------------------------------------------------------------      
      subroutine gpi_waitsome_face(rank_rem,notification_begin)
      use, intrinsic :: ISO_C_BINDING
      use gaspi

      include 'SIZE'
      include 'TOTAL'
      integer np_dir,notification_begin,rank_rem
      integer(gaspi_return_t) :: ret
      integer(gaspi_notification_id_t):: notif_begin,
     $     first_id
      integer(gaspi_number_t)::notif_num
      integer(gaspi_notification_t)::recv_notif_val
      notif_num = 1
      recv_notif_val = 0

      if(rank_rem.ge.0) then
         notif_begin = notification_begin
         ret = gaspi_notify_waitsome(INT(0,1),notif_begin,notif_num,
     $        first_id,10000)
         if(ret.eq.GASPI_TIMEOUT) then
            write(*,*)'notif_begin ',notif_begin,
     $           ' gaspi_notify_waitsome timeout'
            call exit(-1)
         endif
         ret = gaspi_notify_reset(INT(0,1), first_id,recv_notif_val)
         if(recv_notif_val.ne.0) then
c     write(*,*) 'recv data in gpi_waitsome_face'                     
         endif
      endif
      return
      end
c-------------------------------------------------------------------------      
      subroutine packing_face(w,piter,pack_data,face)
      include 'SIZE'

      real w(1)
c      real w(nx1*ny1*nz1*nelt)                                                                      
      integer piter
      real*8 pack_data(1)
      integer*8 face(1)

      integer i
      do i=1,piter
         pack_data(i) = w(face(i))
c         if(nid.eq.0) then                                                                          
c            write(*,*) ' pack data ', pack_data(i),' face(i) ',face(i)                              
c         endif                                                                                      
      enddo


      return
      end
c-------------------------------------------------------------------------                           
      subroutine sum_face(w,piter,msg_data,face,rank_rem)
      include 'SIZE'

      real w(1)
c      real w(nx1*ny1*nz1*nelt)                                                                      
      integer piter
      real*8 msg_data(1)
      integer*8 face(1)
      integer rank_rem

      integer i

      if(rank_rem.ge.0) then
         do i=1,piter
            w(face(i)) = w(face(i)) + msg_data(i)
         enddo
      endif


      return
      end
c-------------------------------------------------------------------------                           
      subroutine local_ox(w,mx,n_shared)

      include 'SIZE'

      integer mx, n_shared
      real w(nx1*ny1*nz1,nelt)

      integer e,i

      if(mx.gt.1) then
         do e=1,nelt-1
            if(mod(e,mx).ne.0)then
               do i=1,n_shared
                  w(l_face_right(i),e) = w(l_face_right(i),e)+
     $                 w(l_face_left(i),e+1)
                  w(l_face_left(i),e+1) = w(l_face_right(i),e)
               enddo
            endif
         enddo
      endif

      return
      end
c-------------------------------------------------------------------------                           
      subroutine local_oy(w,mx,my,mz,n_shared)
      include 'SIZE'

      integer mx, my, mz
      integer n_shared
      real w(nx1*ny1*nz1,nelt)

      integer e_x,e_y,e_z
      integer e, e_top

      if(my.gt.1) then
         do e_z =1,mz
            do e_x = 1, mx
               do e_y = 1,my-1
                  e = e_x + mx*(e_y - 1) + mx*my*(e_z-1)
                  e_top = e + mx
                  do i=1,n_shared
                     w(l_face_y_top(i),e) = w(l_face_y_top(i),e)+
     $                    w(l_face_y_bottom(i),e_top)
                     w(l_face_y_bottom(i),e_top) = w(l_face_y_top(i),e)
                  enddo
               enddo
            enddo
         enddo
      endif

      return
      end
c-------------------------------------------------------------------------                           
      subroutine local_oz(w,mx,my,mz,n_shared)
      include 'SIZE'

      integer mx, my, mz
      integer n_shared
      real w(nx1*ny1*nz1,nelt)

      integer e_x,e_y,e_z
      integer e_back,e_front


      if(mz.gt.1)then
         do e_z =1,mz-1
            do e_y = 1,my
               do e_x = 1,mx
                  e_back = e_x + mx*(e_y-1)+mx*my*(e_z-1)
                  e_front = e_back+mx*my
                  do i=1,n_shared
                     w(l_face_z_back(i),e_back)=
     $                    w(l_face_z_back(i),e_back)+
     $                    w(l_face_z_front(i),e_front)
                     w(l_face_z_front(i),e_front) =
     $                    w(l_face_z_back(i),e_back)
                  enddo
               enddo
            enddo
         enddo
      endif

      return
      end
c-------------------------------------------------------------------------                           
