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
