!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                        ²ú Éú(0,1) ¾ùÔÈ·Ö²¼Ëæ»úÊý×Ó³ÌÐò
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      function ran1(idum)
      integer idum,ia,im,iq,ir,ntab,ndiv
       real ran1,am,eps,rnmx
      parameter(ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,&
                ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)
      integer j,k,iv(ntab),iy
      save iv,iy
      data iv/ntab*0/,iy/0/
      if (idum.le.0.or.iy.eq.0) then
          idum=max(-idum,1)
          do j=ntab+8,1,-1
             k=idum/iq
             idum=ia*(idum-k*iq)-ir*k
             if (idum.lt.0) idum=idum+im
             if (j.le.ntab) iv(j)=idum
          end do
          iy=iv(1)
      end if
      k=idum/iq
        idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=min(am*iy,rnmx)
      end
!     -----------------------------------------------------------------------