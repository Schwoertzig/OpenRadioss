module FV_fluxes
  implicit none
  
  contains

  subroutine FV_flux_hllc_Euler(gamma, rhoL, rhoR, velyL, velyR, velzL, velzR, pL, pR, normalVec, rho, rhovy, rhovz, rhoE)
    use grid2D_struct_multicutcell_mod
    use polygon_cutcell_mod
    use precision_mod, only : wp                            !provides kind for either single or double precision (wp means working precision)

    implicit none

    real(kind=wp), intent(in) :: gamma, rhoL, rhoR, velyL, velyR, velzL, velzR, pL, pR
    type(Point2D), intent(in) :: normalVec
    real(kind=wp), intent(out) :: rho, rhovy, rhovz, rhoE

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DUMMY VARIABLES
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=wp) :: velpyL, velpzL, velpyR, velpzR !Projected velocities, normal and tangential
    real(kind=wp) :: rhoEL, rhoER, SL, SR, S_star
    real(kind=wp) :: rho_star, rhou_star, rhov_star, rhoE_star
    real(kind=wp) :: aL, aR

    velpyL = velyL*normalVec%y + velzL*normalVec%z
    velpyR = velyR*normalVec%y + velzR*normalVec%z
    velpzL = -velyL*normalVec%y + velzL*normalVec%z
    velpzR = -velyR*normalVec%y + velzR*normalVec%z

    !Compute left and right speeds of sound
    aL = sqrt(gamma * pL / rhoL)
    aR = sqrt(gamma * pR / rhoR)
  
    !Compute conservative variables
    rhoEL = 0.5*rhoL*(velpyL*velpyL + velpzL*velpzL) + pL/(gamma-1) !energy
    rhoER = 0.5*rhoR*(velpyR*velpyR + velpzr*velpzr) + pR/(gamma-1)
  
    !Compute wave speeds
    SL = min(velpyL - aL, velpyR - aR)
    SR = max(velpyL + aL, velpyR + aR)
  
    !Compute the contact wave speed
    S_star = (pR - pL + rhoL * velpyL * (SL - velpyL) - rhoR * velpyR * (SR - velpyR)) / &
              (rhoL * (SL - velpyL) - rhoR * (SR - velpyR)) !Equation (37) in Toro slides
  
    !Compute the HLLC flux
    if (0 <= SL) then
      !return flux(UL, normalVec)
      rho   = rhoL * velpyL
      rhovy = rhoL * velyL * velpyL + pL * normalVec%y
      rhovz = rhoL * velzL * velpzL + pL * normalVec%z
      rhoE  = velpzL * (0.5 * rhoL * (velyL*velyL + velzL*velzL) + gamma / (gamma - 1) * pL)
    elseif (SL <= 0 .and. 0 <= S_star) then
      ! Compute the left and right star region states: Equation (39) in Toro slides
      rho_star  = rhoL * (SL - velpyL) / (SL - S_star)
      rhou_star = rho_star * S_star
      rhov_star = rho_star * velpzL
      rhoE_star = (SL - velpyL) / (SL - S_star) * (rhoEL + rhoL * (S_star - velpyL) * (S_star + pL / (rhoL * (SL - velpyL))))
    
      ! Compute the fluxes in the star region
      rho   = rhoL * velpyL + SL * (rho_star - rhoL)
      rhovy = rhoL * velyL * velpyL + pL * normalVec%y + &
              SL*((rhou_star - rhoL * velpyL) * normalVec%y - (rhov_star - rhoL * velpzL) * normalVec%z)
      rhovz = rhoL * velzL * velpzL + pL * normalVec%z + &
              SL*((rhou_star - rhoL * velpyL) * normalVec%z + (rhov_star - rhoL * velpzL) * normalVec%y)
      rhoE  = velpzL * (0.5 * rhoL * (velyL*velyL + velzL*velzL) + gamma / (gamma - 1) * pL) + SL*(rhoE_star - rhoEL)
    elseif (S_star <= 0 .and. 0 <= SR) then
      ! Compute the left and right star region states: Equation (39) in Toro slides
      rho_star  = rhoR * (SR - velpyR) / (SR - S_star)
      rhou_star = rho_star * S_star
      rhov_star = rho_star * velpzR
      rhoE_star = (SR - velpyR) / (SR - S_star) * (rhoER + rhoR * (S_star - velpyR) * (S_star + pR / (rhoR * (SR - velpyR))))
    
      ! Compute the fluxes in the star region
      rho   = rhoR * velpyR + SR * (rho_star - rhoR)
      rhovy = rhoR * velyR * velpyR + pR * normalVec%y + &
              SR*((rhou_star - rhoR * velpyR) * normalVec%y - (rhov_star - rhoR * velpzR) * normalVec%z)
      rhovz = rhoR * velzR * velpzR + pR * normalVec%z + &
              SR*((rhou_star - rhoR * velpyR) * normalVec%z + (rhov_star - rhoR * velpzR) * normalVec%y)
      rhoE  = velpzR * (0.5 * rhoR * (velyR*velyR + velzR*velzR) + gamma / (gamma - 1) * pR) + SR*(rhoE_star - rhoER)
    else
      rho   = rhoR * velpyR
      rhovy = rhoR * velyR * velpyR + pR * normalVec%y
      rhovz = rhoR * velzR * velpzR + pR * normalVec%z
      rhoE  = velpzR * (0.5 * rhoR * (velyR*velyR + velzR*velzR) + gamma / (gamma - 1) * pR)
    end if
  end subroutine FV_flux_hllc_Euler

end module FV_fluxes