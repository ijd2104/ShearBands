c     subroutine getsbfl, evaluate flow law function g {{{1    
      subroutine getsbfl(d,g,strs_eff,Tem,eqps,dgds,dgdT,dgdp,matlaw
     &                    ,fluc,pmult,isw)
c-----------------------------------------------------------------------------c
c                                                                             c
c     Inputs:                                                                 c
c         - d        : material parameters                                    c
c         - strs_eff : effective stress                                       c
c         - Tem      : temperature at the current gp                          c
c         - eqps     : Equivalent Plastic Strain at the current gp            c
c         - matlaw   : material law (LT, ML, LD, etc...)                      c
c         - fluc     : material properties fluctuation (for imperfections)    c
c         - isw      : element task switch                                    c
c                                                                             c
c     Outputs:                                                                c
c         - g        : plastic flow at the current gp                         c
c         - dgds     : derivative of the flow law wrt the stress              c
c         - dgdT     : derivative of the flow law wrt the temperature         c
c         - dgdp     : derivative of the flow law wrt the eqps                c
c         - pmult    : plastic multiplier at the current gp                   c
c                                                                             c
c-----------------------------------------------------------------------------c

      !Evaluate flowlaw and if isw = 3, derivatives of flowlaw 
      implicit none

      real*8  d(*),strs_eff,Tem,eqps,dgds,dgdT,dgdp,g,matlaw,Ts,fluc
      real*8  sigy0,epsy0,pmult,shr,dnm,ths,rprm,bigf,gh,smlf,rpp
      real*8  strs_effd,Temd,kab,ddnmdT,ddnmdp,g1
      real*8  dfdx,dqdp,drdt,dxdp,dxds,dxdt,q,r,x,tf,tfab,tfml
      integer isw
      logical naninfck
      real*8  gr0,c,B,N,T0,Tm,m

      !get unscaled stress and temperature for flow law
      strs_effd = strs_eff
      Temd      = Tem
      
          gr0   = d(8)
          c     = d(9)
          sigy0 = fluc*d(10) !=A
          B     = fluc*d(11)
          N     = d(12)
          T0    = d(13)
          Tm    = d(14)
          m     = d(15)

          Ts    = (Temd-T0)/(Tm-T0)
          ths   = 1.d0 - Ts**m
          shr   = sigy0+B*eqps**N
          dnm   = ths*shr

          g = gr0*exp((strs_effd/dnm-1.d0)/c)

          if (g.lt.2.2204e-16) then
             g     = 0.d0
             pmult = 0.d0
             dgds  = 0.d0
             dgdT  = 0.d0
             dgdp  = 0.d0
          else
              if (strs_effd.lt. 2.2204e-16)then
                  pmult = 0.d0
              else
                  pmult = 1.5*g/strs_effd
              endif
             if (isw.eq.3) then
                dgds =g/(dnm*c)

                dgdT =-(g*d(15)*strs_effd*Ts**(d(15)-1.d0))/((d(13)-
     &               d(14))*d(9)*shr*(Ts**d(15)-1.d0)**2.d0)
                
                if (eqps .eq. 0.d0) then
                   dgdp = 0.d0
                else
                   dgdp=(g*d(11)*d(12)*strs_effd*eqps**(d(12)-1.d0))/
     &                  (d(9)*((Ts**d(15))-1.d0)*
     &                  (d(11)*eqps**d(12)+sigy0)**2)

                endif
             endif
          endif
     
      end
