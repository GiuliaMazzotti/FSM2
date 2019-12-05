!-----------------------------------------------------------------------
! Cumulate fluxes
!-----------------------------------------------------------------------
subroutine CUMULATE(alb,G,Gsoil,H,Hsrf,LE,LEsrf,Melt,Rnet,Roff,Rsrf, &
LWsci,LWsrf,LWveg,SWsci,SWsrf,SWveg,intcpt,unload,Sbsrf,Sbveg, &
KH,KHa,KHg,KHv,KWg,KWv)

use CONSTANTS, only: &
  Lf                  ! Latent heat of fusion (J/kg)

use DIAGNOSTICS, only: &
  Nave,              &! Number of timesteps in average outputs
  diags,             &! Cumulated diagnostics
  SWin,              &! Cumulated incoming solar radiation (J/m^2)
  SWout               ! Cumulated reflected solar radiation (J/m^2)

use DRIVING, only: &
  dt,                &! Timestep (s)
  SW                  ! Incoming shortwave radiation (W/m^2)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

implicit none

real, intent(in) :: &
  alb(Nx,Ny),        &! Albedo
  G(Nx,Ny),          &! Heat flux into surface (W/m^2)
  Gsoil(Nx,Ny),      &! Heat flux into soil (W/m^2)
  H(Nx,Ny),          &! Sensible heat flux to the atmosphere (W/m^2)
  Hsrf(Nx,Ny),       &! Sensible heat flux from the surface (W/m^2)
  LE(Nx,Ny),         &! Latent heat flux to the atmosphere (W/m^2)
  LEsrf(Nx,Ny),      &! Latent heat flux from the surface (W/m^2)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  Rnet(Nx,Ny),       &! Net radiation (W/m^2)
  Roff(Nx,Ny),       &! Runoff from snow (kg/m^2)
  Rsrf(Nx,Ny),       &! Net radiation absorbed by the surface (W/m^2)
  LWsci(Nx,Ny),      &! Subcanopy incoming LWR (W/m^2)
  LWsrf(Nx,Ny),      &! Net LW radiation absorbed by the surface (W/m^2)
  LWveg(Nx,Ny),      &! Net LW radiation absorbed by vegetation (W/m^2)
  SWsci(Nx,Ny),      &! Subcanopy incoming SWR (W/m^2)
  SWsrf(Nx,Ny),      &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny),      &! Net SW radiation absorbed by vegetation (W/m^2)
  intcpt(Nx,Ny),     &! Canopy interception (kg/m^2)
  unload(Nx,Ny),     &! Snow mass unloaded from canopy (kg/m^2) 
  Sbsrf(Nx,Ny),      &! Sublimation from the surface (kg/m^2)
  Sbveg(Nx,Ny),      &! Sublimation from the vegetation (kg/m^2)
  KH(Nx,Ny),         &! Eddy diffusivity for heat to the atmosphere (m/s)
  KHa(Nx,Ny),        &! Eddy diffusivity from the canopy air space (m/s)
  KHg(Nx,Ny),        &! Eddy diffusivity for heat from the ground (m/s)
  KHv(Nx,Ny),        &! Eddy diffusivity for heat from vegetation (m/s)
  KWg(Nx,Ny),        &! Eddy diffusivity for water from the ground (m/s)
  KWv(Nx,Ny)          ! Eddy diffusivity for water from vegetation (m/s)


SWin(:,:) = SWin (:,:)+ SW(:,:)*dt
SWout(:,:) = SWout(:,:) + alb(:,:)*SW(:,:)*dt
diags(:,:,1) = diags(:,:,1) + G(:,:)
diags(:,:,2) = diags(:,:,2) + Gsoil(:,:)
diags(:,:,3) = diags(:,:,3) + H(:,:)
diags(:,:,4) = diags(:,:,4) + Hsrf(:,:)
diags(:,:,5) = diags(:,:,5) + LE(:,:)
diags(:,:,6) = diags(:,:,6) + LEsrf(:,:)
diags(:,:,7) = diags(:,:,7) + Lf*Melt(:,:)
diags(:,:,8) = diags(:,:,8) + Rnet(:,:)
diags(:,:,9) = diags(:,:,9) + Roff(:,:) * Nave
diags(:,:,10) = diags(:,:,10) + Rsrf(:,:)
diags(:,:,11) = diags(:,:,11) + LWsci(:,:)
diags(:,:,12) = diags(:,:,12) + LWsrf(:,:)
diags(:,:,13) = diags(:,:,13) + LWveg(:,:)
diags(:,:,14) = diags(:,:,14) + SWsci(:,:)
diags(:,:,15) = diags(:,:,15) + SWsrf(:,:)
diags(:,:,16) = diags(:,:,16) + SWveg(:,:)
diags(:,:,17) = diags(:,:,17) + intcpt(:,:) * Nave
diags(:,:,18) = diags(:,:,18) + unload(:,:) * Nave
diags(:,:,19) = diags(:,:,19) + Sbsrf(:,:)
diags(:,:,20) = diags(:,:,20) + Sbveg(:,:)
diags(:,:,21) = diags(:,:,21) + KH(:,:)
diags(:,:,22) = diags(:,:,22) + KHa(:,:)
diags(:,:,23) = diags(:,:,23) + KHg(:,:)
diags(:,:,24) = diags(:,:,24) + KHv(:,:)
diags(:,:,25) = diags(:,:,25) + KWg(:,:)
diags(:,:,26) = diags(:,:,26) + KWv(:,:)

! snow albedo: calculate offline from SWsci / SWsrf 
! Hveg, LEveg: calculate as H-Hsrf, LE-LEsrf 

end subroutine CUMULATE
