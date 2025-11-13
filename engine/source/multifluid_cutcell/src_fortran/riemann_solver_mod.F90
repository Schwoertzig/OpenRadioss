#include "my_real.inc"
module riemann_solver_mod
contains
  function fl(pL, aL, gamma, p_star)
    implicit none
    my_real :: pL, aL, gamma, p_star
    my_real :: fl
    if (p_star > pL) then
      ! Left shock
      fl = 2 * aL / sqrt(2 * gamma * (gamma - 1)) * (1 - p_star / pL) / sqrt(1 + (gamma + 1) / (gamma - 1) * p_star / pL) !Leveque (14.52)
    else
      ! Left rarefaction
      fl = 2 * aL / (gamma - 1) * (1 - (p_star / pL)**((gamma - 1) / (2 * gamma))) !Leveque (14.51)
    end if
  end function

  ! Gradient of function fl
  function Dfl(pL, aL, gamma, p_star)
    implicit none
    my_real :: pL, aL, gamma, p_star
    my_real :: Dfl
    if (p_star > pL) then
      ! Left shock
      Dfl = -2 * aL / sqrt(2 * gamma * (gamma - 1)) / sqrt(1 + (gamma + 1) / (gamma - 1) * p_star / pL) &
            * (1.0 + 0.5 * (gamma + 1) / (gamma - 1) / (1 + (gamma + 1) / (gamma - 1) * p_star / pL)) / pL
    else
      ! Left rarefaction
      Dfl = -aL / (gamma * pL) * (p_star / pL)**(-(gamma + 1) / (2 * gamma))
    end if
  end function Dfl

  function fr(pR, aR, gamma, p_star)
    implicit none
    my_real :: pR, aR, gamma, p_star
    my_real ::  fr
    if (p_star > pR) then
      ! Right shock
      fr = 2 * aR / sqrt(2 * gamma * (gamma - 1)) * (1 - p_star / pR) / sqrt(1 + (gamma + 1) / (gamma - 1) * p_star / pR) !Leveque (14.55)
    else
      ! Right rarefaction
      fr = 2 * aR / (gamma - 1) * (1 - (p_star / pR)**((gamma - 1) / (2 * gamma))) !Leveque (14.54)
    end if
  end function fr

  ! Gradient of function fr
  function Dfr(pR, aR, gamma, p_star)
    implicit none
    my_real :: pR, aR, gamma, p_star
    my_real ::  Dfr
    if (p_star > pR) then
      ! Left shock
      Dfr = -2 * aR / sqrt(2 * gamma * (gamma - 1)) / sqrt(1 + (gamma + 1) / (gamma - 1) * p_star / pR) &
              * (1.0 + 0.5 * (gamma + 1) / (gamma - 1) / (1 + (gamma + 1) / (gamma - 1) * p_star / pR)) / pR
    else
      ! Left rarefaction
      Dfr = -aR / (gamma * pR) * (p_star / pR)**(-(gamma + 1) / (2 * gamma))
    end if 
  end function Dfr

  ! Function to solve for pressure in the star region
  function phi(pL, aL, uL, pR, aR, uR, gamma, p_star)
    implicit none 
    my_real :: pL, aL, uL 
    my_real :: pR, aR, uR 
    my_real :: gamma, p_star
    my_real :: phi

    phi = fl(pL, aL, gamma, p_star) + fr(pR, aR, gamma, p_star) + uL - uR
  end function phi

  ! Gradient of function phi
  function Dphi(pL, aL, uL, pR, aR, uR, gamma, p_star)
    implicit none 
    my_real :: pL, aL, uL 
    my_real :: pR, aR, uR 
    my_real :: gamma, p_star
    my_real :: Dphi
   
    Dphi = Dfl(pL, aL, gamma, p_star) + Dfr(pR, aR, gamma, p_star)
  end function Dphi


  ! Use a Newton method to find the root
  function find_root(pL, aL, uL, pR, aR, uR, gamma, eps) result(p)
    implicit none 
    my_real :: pL, aL, uL 
    my_real :: pR, aR, uR 
    my_real :: gamma, eps
    my_real :: p 

    my_real :: p_prev, res

    p = pL
    p_prev = pR
    do while (abs(p - p_prev) > eps)
      res = phi(pL, aL, uL, pR, aR, uR, gamma, p)
      p_prev = p
      p = p - res / Dphi(pL, aL, uL, pR, aR, uR, gamma, p)
    end do
  end function find_root

  ! Function to solve the Riemann problem for the Euler equations and return wave velocity and flux at the interface
  ! wave-types: 1=left shock, 2=contact discontinuity, 3=right shock
  subroutine solve_riemann_problem(gamma, rhoL, rhoR, velyL, velyR, velzL, velzR, pL, pR, wave_type, &
                                  normalVecy, normalVecz, us, vsL, vsR, ps)
    implicit none

    my_real :: gamma, rhoL, rhoR, velyL, velyR, velzL, velzR, pL, pR
    integer :: wave_type
    my_real :: normalVecy, normalVecz
    my_real :: us, vsL, vsR, ps

    my_real, parameter :: prec_root_find = 1e-10

    my_real :: uL, uR, vL, vR, aL, aR
    my_real :: SL, SR
    my_real :: p_star, u_star
    my_real :: rho_starL, rho_starR, v_starL, v_starR

    uL = velyL*normalVecy + velzL*normalVecz
    uR = velyR*normalVecy + velzR*normalVecz

    vL = -velyL*normalVecz + velzL*normalVecy
    vR = -velyR*normalVecz + velzR*normalVecy

    !gamma = UL.gamma #TODO: change this to accomodate differents gammas for different materials

    ! Compute primitive variables
    aL = sqrt(gamma * pL / rhoL)
    aR = sqrt(gamma * pR / rhoR)

    ! Find the pressure in the star region
    p_star = find_root(pL, aL, uL, pR, aR, uR, gamma, prec_root_find)

    ! Compute the velocity in the star region
    u_star = uL + fl(pL, aL, gamma, p_star)


    !Determine left and right shock velocities
    SL = uL - aL * sqrt((gamma + 1) / (2 * gamma) * p_star / pL + (gamma - 1) / (2 * gamma))
    SR = uR + aR * sqrt((gamma + 1) / (2 * gamma) * p_star / pR + (gamma - 1) / (2 * gamma))

    ! Compute densities in the star region
    if (p_star > pL) then
      ! Left shock
      rho_starL = rhoL * ((p_star / pL) + (gamma - 1) / (gamma + 1)) / (((gamma - 1) / (gamma + 1)) * (p_star / pL) + 1)
      v_starL = SL + rhoL / rho_starL * (vL - SL)
    else
      ! Left rarefaction
      rho_starL = rhoL * (p_star / pL)**(1 / gamma)
      v_starL = vL
    end if

    if (p_star > pR) then
      !Right shock
      rho_starR = rhoR * ((p_star / pR) + (gamma - 1) / (gamma + 1)) / (((gamma - 1) / (gamma + 1)) * (p_star / pR) + 1)
      v_starR = SR + rhoR / rho_starR * (vR - SR)
    else
      ! Right rarefaction
      rho_starR = rhoR * (p_star / pR)**(1 / gamma)
      v_starR = vR
    end if

    if (wave_type == 2) then
      us = u_star
      vsL = v_starL
      vsR = v_starR
      ps = p_star
    elseif (wave_type == 1) then
      us = SL
      vsL = vL
      vsR = v_starL
      ps = pL
    elseif (wave_type == 3) then
      us = SR
      vsL = v_starR
      vsR = vR
      ps = pR
    end if
  end subroutine solve_riemann_problem
end module riemann_solver_mod
