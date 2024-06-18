c
c...................................................................
c
      subroutine distaz(ep_la,ep_lo,s_la,s_lo,dist,azi,baz,id_coord)
c
c     subroutine distaz(epla,eplo,sla,slo,del,dist,azi,baz)
c     a subroutine to calculate great circle distances and
c     azimuth between two points on the earth'h surface. the
c     earth is assummed to be an ellipsoid of revolution.
c     this routine is from a program by lynn peseckis, 1979.
c
c     input parameters :
c
c     epla, eplo ........latitude and longitude of first point
c                        on earth's surface. north latitude
c                        and east longitude is positive, south
c                        latitude and west longitude is negative.
c     sla, slo ..........latitude and longitude of second point
c                        on earth's surface.
c
c     returned parameters :
c
c     del ...............distance in degrees between two points.
c     dist ..............distance in kilometers between two points.
c     az ................azimuth from first point to second.
c     baz ...............back azimuth from second point to first.
c
c                                  Compiled by Yuehua Zeng at USC
c
      parameter (pi=3.141592654,pi5=0.5*pi,dr=180./pi)
      parameter (re=6378.155, geocen = 0.99327733)
      real*8 c,bp,ap,gamma,cbp,sbp,cap,sap,abtem,altem,baztem
      if(id_coord .eq. 2) then
        dx=s_la-ep_la
        dy=s_lo-ep_lo
        dist=sqrt(dx**2+dy**2)
        azi=acos(abs(dx)/dist)*dr
        if(dx.lt.0.0) azi=180.0-azi
        if(dy.lt. 0.0) azi=360.-azi
        baz=azi+180.0
        return
      endif

      write(*,*) ep_la,ep_lo,s_la,s_lo        

      epla= ep_la/dr
      eplo= ep_lo/dr
      sla= s_la/dr
      slo= s_lo/dr
c
c  begin calculation of distance
      dlon=abs(eplo-slo)
      if(abs(epla)-pi5 .eq. 0.0)then
        bp=pi5-epla
        ap=pi5-(atan(geocen*tan(sla)))
        if(epla) 170,170,180
      else
        bp=pi5-(atan(geocen*tan(epla)))
      endif
      if(abs(sla)-pi5 .eq. 0.0) then
        ap=pi5-sla
        if(sla) 180,170,170
      else
        ap=pi5-(atan(geocen*tan(sla)))
      endif
      if(dlon-0.00001.gt.0.0) goto 200
      if(epla-sla)170,170,180
170   azi=0.0
      baz=pi
      goto 190
180   azi=pi
      baz=0.0
190   del=abs(ap-bp)
      goto 210
200   gamma=slo-eplo
      cbp=cos(bp)
      sbp=sin(bp)
      cap=cos(ap)
      sap=sin(ap)
      abtem=cbp*cap+sbp*sap*cos(gamma)
      c=atan(sqrt((1.0-abtem)*(1.+abtem))/abtem)
      if(abtem.lt.0.0)c=pi+c
      del=c
      altem=(cap-abtem*cbp)/(sin(c)*sbp)
      azi=atan(sqrt((1.0-altem)*(1.+altem))/altem)
      if(altem.lt.0.0)azi=pi+azi
      baztem=(cbp-abtem*cap)/(sin(c)*sap)
      baz=atan(sqrt((1.0-baztem)*(1.+baztem))/baztem)
      if(baztem.lt.0.0)baz=pi+baz
      if(sin(gamma).lt.0.0)then
        azi=2.0*pi-azi
      else
        baz=2.0*pi-baz
      endif
210   dist = del*re
      azi=azi*dr
      baz=baz*dr
c
      return
      end

