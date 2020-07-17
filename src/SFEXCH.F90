!-----------------------------------------------------------------------
! Surface exchange coefficients
!-----------------------------------------------------------------------
subroutine SFEXCH(fsnow,gs1,KH,KHa,KHg,KHv,KWg,KWv,Uze)

#include "OPTS.h"

use CONSTANTS, only: &
  grav,              &! Acceleration due to gravity (m/s^2)
  vkman               ! Von Karman constant

use DRIVING, only: &
  Ta,                &! Air temperature (K)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Ua,                &! Wind speed (m/s)
  zT,                &! Temperature and humidity measurement height (m)
  zU                  ! Wind measurement height (m)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only : &
  bstb,              &! Atmospheric stability parameter
  cden,              &! Dense canopy turbulent transfer coefficient
  cveg,              &! Vegetation turbulent transfer coefficient
  gsnf,              &! Snow-free vegetation moisture conductance (m/s)
  rchd,              &! Ratio of displacement height to canopy height
  rchz,              &! Ratio of roughness length to canopy height
  z0sn,              &! Snow roughness length (m)
  wcan,              &! Canopy wind decay coefficient
  zsub,              &! Sub-canopy reference height (m)
  zgf,               &! z0g canopy dependence factors (-)
  zgr,               &! z0g canopy dependence range (-)
  khcf                ! diffusivity correction factor 

use PARAMMAPS, only: &
  hcan,              &! Canopy height (m)
  fveg,              &! Canopy cover fraction
  VAI,               &! Vegetation area index
  fsky,              &! Sky view fraction
  trcn,              &! Canopy transmissivity
  z0sf,              &! Snow-free surface roughness length (m)
  fves                ! stand-scale veg fraction

use STATE_VARIABLES, only : &
  Qcan,              &! Canopy air space humidity
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sveg,              &! Snow mass on vegetation (kg/m^2)
  Tcan,              &! Canopy air space temperature (K)
  Tsrf,              &! Surface skin temperature (K)
  Tveg                ! Vegetation temperature (K)

implicit none

real, intent(in) :: &
  fsnow(Nx,Ny),      &! Ground snowcover fraction
  gs1(Nx,Ny)          ! Surface moisture conductance (m/s)

real, intent(out) :: &
  KH(Nx,Ny),         &! Eddy diffusivity for heat to the atmosphere (m/s)
  KHa(Nx,Ny),        &! Eddy diffusivity from the canopy air space (m/s)
  KHg(Nx,Ny),        &! Eddy diffusivity for heat from the ground (m/s)
  KHv(Nx,Ny),        &! Eddy diffusivity for heat from vegetation (m/s)
  KWg(Nx,Ny),        &! Eddy diffusivity for water from the ground (m/s)
  KWv(Nx,Ny),        &! Eddy diffusivity for water from vegetation (m/s)
  Uze(Nx,Ny)          ! Local wind speed at reference height z1 (m/s)

integer :: &
  i,j                 ! Point counters

real :: &
  CD,                &! Drag coefficient
  dh,                &! Displacement height (m)
  fh,                &! Stability factor
  Qs,                &! Saturation humidity
  RiB,               &! Bulk Richardson number
  Ric,               &! Sub-canopy Richardson number
  Tint,              &! Interpolated canopy - ground temperature (K)
  ustar,             &! Friction velocity (m/s)
  zT1,               &! Temperature measurement height with offset (m)
  zU1,               &! Wind measurement height with offset (m)
  z0,                &! Roughness length for momentum (m)
  z0g,               &! Ground surface roughness length (m)
  z0h,               &! Roughness length for heat (m)
  z0v,               &! Vegetation roughness length (m)
  Uso,               &! Wind speed at sub-canopy reference height in open (m/s)
  Usf,               &! Wind speed at sub-canopy reference height in dense forest (m/s)
  Usub,              &! Wind speed at sub-canopy reference height (m/s)
  Uh,                &! Wind velocity at canopy top (m/s)
  KHh,               &! Eddy diffusivity at canopy top (m/s)
  rad,               &! Aerodynamic resistance between canopy air space and atmosphere in dense canopy, at reference height (s/m)
  rao,               &! Theoretical aerodynamic resistance at reference height for an open site (s/m)  
  rgd,               &! Aerodynamic resistance between canopy air space and ground in dense canopy, at reference height (s/m)
  rgo,               &! Aerodynamic resistance between canopy air space and ground in open, at reference height (s/m)
  Uc                  ! Wind speed in canopy layer (at height of turbulent flux from veg to cas) (m/s)

do j = 1, Ny
do i = 1, Nx
  
! ground roughness length
  z0g = (z0sn**(fsnow(i,j)**1)) * (z0sf(i,j)**(1 - (fsnow(i,j)**1))) 
                            
#if ZOFFST == 0
! Heights specified above ground
  zU1 = zU
  zT1 = zT
#endif
#if ZOFFST == 1
! Heights specified above canopy top
  zU1 = zU + hcan(i,j)
  zT1 = zT + hcan(i,j)
#endif

#if EXCHNG != 2
! Roughness lengths and friction velocity
  z0v = rchz*hcan(i,j)
  z0  = (z0v**fveg(i,j)) * (z0g**(1 - fveg(i,j)))
  z0h = 0.1*z0
  dh = fveg(i,j)*rchd*hcan(i,j)
  CD = (vkman / log((zU1 - dh)/z0))**2
  ustar = sqrt(CD)*Ua(i,j)
#endif

#if EXCHNG == 2
  ! Open site calculations
  z0h = 0.1 * z0g 
  Uso = Ua(i,j)*log(zsub/z0g)/log(zU/z0g)
  ustar = vkman*Ua(i,j)/log(zU/z0g) ! open site friction velocity
  rgo = log(zT/z0h)/(vkman*ustar) 

  ! Forest site calculations
  if (fveg(i,j) > 0) then
    z0g = (zgf + zgr*fveg(i,j)) * z0g    
    z0h = 0.1 * z0g
    dh = rchd * hcan(i,j) 
    z0v = rchz * hcan(i,j) 
    ustar = vkman*Ua(i,j)/log((zU1 - dh)/z0v)
    Uh = (ustar/vkman)*log((hcan(i,j) - dh)/z0v)
    KHh = vkman*ustar*(hcan(i,j) - dh)
    Usf = exp(wcan*(zsub/hcan(i,j) - 1))*Uh
  end if
#endif 

#if EXCHNG == 0
! No stability adjustment
  fh = 1
  Ric = 0
#endif
#if EXCHNG == 1
! Stability adjustment (Louis et al. 1982, quoted by Beljaars 1992)
  Tint = fveg(i,j)*Tveg(i,j) + (1 - fveg(i,j))*Tsrf(i,j)
  RiB = grav*(Ta(i,j) - Tint)*(zU1 - dh)**2 / ((zT1 - dh)*Ta(i,j)*Ua(i,j)**2)
  if (RiB > 0) then 
    fh = 1/(1 + 3*bstb*RiB*sqrt(1 + bstb*RiB))
  else
    fh = 1 - 3*bstb*RiB / (1 + 3*bstb**2*CD*sqrt(-RiB*zU1/z0))
  end if
  Ric = grav*(Tcan(i,j) - Tsrf(i,j))*hcan(i,j) / (Tcan(i,j)*ustar**2)
  Ric = max(min(Ric,10.),0.)
#endif

! Eddy diffusivities
  if (fveg(i,j) == 0) then
    KH(i,j) = fh*vkman*ustar / log(zT1/z0h)
    call QSAT(Ps(i,j),Tsrf(i,j),Qs)
    if (Sice(1,i,j) > 0 .or. Qa(i,j) > Qs) then
      KWg(i,j) = KH(i,j)
    else
      KWg(i,j) = gs1(i,j)*KH(i,j) / (gs1(i,j) + KH(i,j))
    end if
    KHv(i,j) = 0
    KWv(i,j) = 0
    KHa(i,j) = 0
    KHg(i,j) = 0
  else
#if EXCHNG != 2
    KHa(i,j) = fh*vkman*ustar / log((zT1 - dh)/z0)
    KHg(i,j) = vkman*ustar*((1 - fveg(i,j))*fh/log(z0/z0h) + fveg(i,j)*cden/(1 + 0.5*Ric))
    KHv(i,j) = sqrt(ustar)*VAI(i,j)/cveg
#endif
#if EXCHNG == 2
    rad = (log((zT1 - dh)/(hcan(i,j)- dh))/(vkman*ustar) +  & 
       hcan(i,j)*(exp(wcan*(1 -(z0v + dh)/hcan(i,j))) - 1)/(wcan*KHh))/khcf
    KHa(i,j) = sqrt(fves(i,j))/ rad 
    Usub = sqrt(fves(i,j))*Usf + (1 - sqrt(fves(i,j)))*Uso
    rgd = 1 / (vkman**2 * Usub) * log(zsub/z0h) * log(zsub/z0g) 
    KHg(i,j)  = 1/rgd
    Uc = exp(wcan*((z0v + dh)/hcan(i,j) - 1))*Uh
    KHv(i,j) = VAI(i,j)*sqrt(Uc)/cveg 
#endif
    call QSAT(Ps(i,j),Tsrf(i,j),Qs)
    if (Qcan(i,j) > Qs) then
      KWg(i,j) = KHg(i,j)
    else
      KWg(i,j) = gs1(i,j)*KHg(i,j) / (gs1(i,j) + KHg(i,j))
    end if
    call QSAT(Ps(i,j),Tveg(i,j),Qs)
    if (Sveg(i,j) > 0 .or. Qcan(i,j) > Qs) then
      KWv(i,j) = KHv(i,j)
    else
      KWv(i,j) = gsnf*KHv(i,j) / (gsnf + KHv(i,j))
    end if
    KH(i,j) = 0 
  end if 
  Uze(i,j) = Usub

#if CANMOD == 0
! Combined resistances for 0-layer canopy model
    KH(i,j) = KHg(i,j)*(KHa(i,j) + KHv(i,j)) / (KHa(i,j) + KHg(i,j) + KHv(i,j))
    KWg(i,j) = KWg(i,j)*(KHa(i,j) + KWv(i,j)) / (KHa(i,j) + KWg(i,j) + KWv(i,j))
#endif

end do
end do

end subroutine SFEXCH
