#include <define.h>

MODULE MOD_Lake

!-----------------------------------------------------------------------
! DESCRIPTION:
! Simulating energy balance processes of land water body
!
! REFERENCE:
! Dai et al, 2018, The lake scheme of the common land model and its performance evaluation.
! Chinese Science Bulletin, 63(28-29), 3002–3021, https://doi.org/10.1360/N972018-00609
!
! Original author: Yongjiu Dai 04/2014/
!
! Revisions:
! Nan Wei,  01/2018: interaction btw prec and lake surface including phase change of prec and water body
! Nan Wei,  06/2018: update heat conductivity of water body and soil below and snow hydrology
! Hua Yuan, 01/2023: added snow layer absorption in melting calculation
!-----------------------------------------------------------------------

 use MOD_Precision
 USE MOD_Vars_Global, only: maxsnl
 USE MOD_SPMD_Task
 IMPLICIT NONE
!  SAVE
 real(R8),parameter :: SHR_CONST_PI     = 3.14159265358979323846_R8
 real(R8),parameter :: SHR_CONST_RHOICE = 0.917e3_R8           !  density of ice (kg/m^3)

 integer, parameter :: iulog = 6                               ! "stdout" log file unit number, default is 6
 integer, parameter :: numrad  = 2                             !  number of solar radiation bands: vis, nir

 logical, parameter :: use_extrasnowlayers = .false.
 character(len=256), parameter :: snow_shape = 'sphere'  
 logical, parameter :: use_dust_snow_internal_mixing = .false.
 character(len=256), parameter :: snicar_atm_type = 'mid-latitude_winter' 
! PUBLIC MEMBER FUNCTIONS:
  public :: newsnow_lake
  public :: laketem
  public :: snowwater_lake
  public :: LakeOptics_init

! PRIVATE MEMBER FUNCTIONS:
  private :: roughness_lake
  private :: hConductivity_lake

  ! !PUBLIC DATA MEMBERS:
  integer,  public, parameter :: sno_nbr_aer = 8         ! number of aerosol species in snowpack
                                                         ! (indices described above) [nbr]
  logical,  public, parameter :: DO_SNO_OC   = .false.   ! parameter to include organic carbon (OC)
                                                         ! in snowpack radiative calculations
  logical,  public, parameter :: DO_SNO_AER  = .true.    ! parameter to include aerosols in snowpack radiative calculations
  ! !PRIVATE DATA MEMBERS:
  integer,  parameter :: numrad_snw     = 5              ! number of spectral bands used in snow model [nbr]
  integer,  parameter :: nir_bnd_bgn    = 2              ! first band index in near-IR spectrum [idx]
  integer,  parameter :: nir_bnd_end    = 5              ! ending near-IR band index [idx]
  integer,  parameter :: idx_Mie_snw_mx = 1471           ! number of effective radius indices used in Mie lookup table [idx]
  integer,  parameter :: idx_T_max      = 11             ! maxiumum temperature index used in aging lookup table [idx]
  integer,  parameter :: idx_T_min      = 1              ! minimum temperature index used in aging lookup table [idx]
  integer,  parameter :: idx_Tgrd_max   = 31             ! maxiumum temperature gradient index used in aging lookup table [idx]
  integer,  parameter :: idx_Tgrd_min   = 1              ! minimum temperature gradient index used in aging lookup table [idx]
  integer,  parameter :: idx_rhos_max   = 8              ! maxiumum snow density index used in aging lookup table [idx]
  integer,  parameter :: idx_rhos_min   = 1              ! minimum snow density index used in aging lookup table [idx]

#ifdef MODAL_AER
  ! NOTE: right now the macro 'MODAL_AER' is not defined anywhere, i.e.,
  ! the below (modal aerosol scheme) is not available and can not be
  ! active either. It depends on the specific input aerosol deposition
  ! data which is suitable for modal scheme. [06/15/2023, Hua Yuan]
  !mgf++
  integer,  parameter :: idx_bc_nclrds_min     = 1       ! minimum index for BC particle size in optics lookup table
  integer,  parameter :: idx_bc_nclrds_max     = 10      ! maximum index for BC particle size in optics lookup table
  integer,  parameter :: idx_bcint_icerds_min  = 1       ! minimum index for snow grain size in optics lookup table for within-ice BC
  integer,  parameter :: idx_bcint_icerds_max  = 8       ! maximum index for snow grain size in optics lookup table for within-ice BC
  !mgf--
#endif

  integer,  parameter :: snw_rds_max_tbl = 1500          ! maximum effective radius defined in Mie lookup table [microns]
  integer,  parameter :: snw_rds_min_tbl = 30            ! minimium effective radius defined in Mie lookup table [microns]
  real(r8), parameter :: snw_rds_max     = 1500._r8      ! maximum allowed snow effective radius [microns]
  real(r8), parameter :: snw_rds_min     = 54.526_r8     ! minimum allowed snow effective radius (also "fresh snow" value) [microns
  real(r8), parameter :: snw_rds_refrz   = 1000._r8      ! effective radius of re-frozen snow [microns]
  real(r8), parameter :: min_snw = 1.0E-30_r8            ! minimum snow mass required for SNICAR RT calculation [kg m-2]
  !real(r8), parameter :: C1_liq_Brun89  = 1.28E-17_r8   ! constant for liquid water grain growth [m3 s-1],
                                                         ! from Brun89
  real(r8), parameter :: C1_liq_Brun89   = 0._r8         ! constant for liquid water grain growth [m3 s-1],
                                                         ! from Brun89: zeroed to accomodate dry snow aging
  real(r8), parameter :: C2_liq_Brun89   = 4.22E-13_r8   ! constant for liquid water grain growth [m3 s-1],
                                                         ! from Brun89: corrected for LWC in units of percent

  real(r8), parameter :: tim_cns_bc_rmv  = 2.2E-8_r8     ! time constant for removal of BC in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)
  real(r8), parameter :: tim_cns_oc_rmv  = 2.2E-8_r8     ! time constant for removal of OC in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)
  real(r8), parameter :: tim_cns_dst_rmv = 2.2E-8_r8     ! time constant for removal of dust in snow on sea-ice
                                                         ! [s-1] (50% mass removal/year)
  !$acc declare copyin(C1_liq_Brun89, C2_liq_Brun89, &
  !$acc tim_cns_bc_rmv, tim_cns_oc_rmv, tim_cns_dst_rmv)

  ! scaling of the snow aging rate (tuning option):
  logical :: flg_snoage_scl    = .false.                 ! flag for scaling the snow aging rate by some arbitrary factor
  real(r8), parameter :: xdrdt = 1.0_r8                  ! arbitrary factor applied to snow aging rate
  ! snow and aerosol Mie parameters:
  ! (arrays declared here, but are set in iniTimeConst)
  ! (idx_Mie_snw_mx is number of snow radii with defined parameters (i.e. from 30um to 1500um))

  ! direct-beam weighted ice optical properties
  real(r8), allocatable :: refindx_im_clr     (:,:) ! (numrad_snw,90);
  real(r8), allocatable :: refindx_re_clr     (:,:) ! (numrad_snw,90);
  real(r8), allocatable :: refindxwat_im_clr     (:,:) ! (numrad_snw,90);
  real(r8), allocatable :: refindxwat_re_clr     (:,:) ! (numrad_snw,90);
  real(r8), allocatable :: ss_alb_snw_avg_clr     (:,:,:) ! (numrad_snw,90);
  real(r8), allocatable :: asm_prm_snw_avg_clr     (:,:,:) ! (numrad_snw,90);
  real(r8), allocatable :: ext_cff_mss_snw_avg_clr     (:,:,:) ! (numrad_snw,90);
  real(r8), allocatable :: sca_cff_vlm_avg_clr     (:,:,:) ! (numrad_snw,90);
  real(r8), allocatable :: asm_prm_ice_avg_clr     (:,:,:) ! (numrad_snw,90);
  real(r8), allocatable :: abs_cff_mss_avg_clr     (:,:,:) ! (numrad_snw,90);
  real(r8), allocatable :: ss_alb_wtr_avg_clr     (:,:) ! (numrad_snw,90);
  real(r8), allocatable :: ext_cff_mss_wtr_avg_clr     (:,:) ! (numrad_snw,90);
  real(r8), allocatable :: refindx_im_cld     (:) ! (numrad_snw);
  real(r8), allocatable :: refindx_re_cld     (:) ! (numrad_snw);
  real(r8), allocatable :: refindxwat_im_cld     (:) ! (numrad_snw);
  real(r8), allocatable :: refindxwat_re_cld    (:) ! (numrad_snw);
  real(r8), allocatable :: ss_alb_snw_avg     (:,:) ! (numrad_snw);
  real(r8), allocatable :: asm_prm_snw_avg     (:,:) ! (numrad_snw);
  real(r8), allocatable :: ext_cff_mss_snw_avg     (:,:) ! (numrad_snw);
  real(r8), allocatable :: sca_cff_vlm_avg    (:,:) ! (numrad_snw);
  real(r8), allocatable :: asm_prm_ice_avg     (:,:) ! (numrad_snw);
  real(r8), allocatable :: abs_cff_mss_avg     (:,:) ! (numrad_snw);

  ! direct & diffuse flux
  real(r8), allocatable :: flx_wgt_dir (:,:,:) ! (6, 90, numrad_snw) ! direct flux, six atmospheric types, 0-89 SZA
  real(r8), allocatable :: flx_wgt_dif (:,:)   ! (6, numrad_snw)     ! diffuse flux, six atmospheric types

  ! snow grain shape
  integer, parameter :: snow_shape_sphere            = 1
  integer, parameter :: snow_shape_spheroid          = 2
  integer, parameter :: snow_shape_hexagonal_plate   = 3
  integer, parameter :: snow_shape_koch_snowflake    = 4

  ! atmospheric condition for SNICAR-AD
  integer, parameter :: atm_type_default             = 0
  integer, parameter :: atm_type_mid_latitude_winter = 1
  integer, parameter :: atm_type_mid_latitude_summer = 2
  integer, parameter :: atm_type_sub_Arctic_winter   = 3
  integer, parameter :: atm_type_sub_Arctic_summer   = 4
  integer, parameter :: atm_type_summit_Greenland    = 5
  integer, parameter :: atm_type_high_mountain       = 6

#ifdef MODAL_AER
  !mgf++
  ! Size-dependent BC optical properties. Currently a fixed BC size is
  ! assumed, but this framework enables optical properties to be
  ! assigned based on the BC effective radius, should this be
  ! implemented in the future.
  !
  ! within-ice BC (i.e., BC that was deposited within hydrometeors)
  real(r8), allocatable :: ss_alb_bc1     (:,:)  ! (numrad_snw,idx_bc_nclrds_max);
  real(r8), allocatable :: asm_prm_bc1    (:,:)  ! (numrad_snw,idx_bc_nclrds_max);
  real(r8), allocatable :: ext_cff_mss_bc1(:,:)  ! (numrad_snw,idx_bc_nclrds_max);

  ! external BC
  real(r8), allocatable :: ss_alb_bc2     (:,:)  ! (numrad_snw,idx_bc_nclrds_max);
  real(r8), allocatable :: asm_prm_bc2    (:,:)  ! (numrad_snw,idx_bc_nclrds_max);
  real(r8), allocatable :: ext_cff_mss_bc2(:,:)  ! (numrad_snw,idx_bc_nclrds_max);
  !mgf--
#else
  ! hydrophiliic BC
  real(r8), allocatable :: ss_alb_bc1     (:) ! (numrad_snw);
  real(r8), allocatable :: asm_prm_bc1    (:) ! (numrad_snw);
  real(r8), allocatable :: ext_cff_mss_bc1(:) ! (numrad_snw);

  ! hydrophobic BC
  real(r8), allocatable :: ss_alb_bc2     (:) ! (numrad_snw);
  real(r8), allocatable :: asm_prm_bc2    (:) ! (numrad_snw);
  real(r8), allocatable :: ext_cff_mss_bc2(:) ! (numrad_snw);
#endif

  ! hydrophobic OC
  real(r8), allocatable :: ss_alb_oc1     (:) ! (numrad_snw);
  real(r8), allocatable :: asm_prm_oc1    (:) ! (numrad_snw);
  real(r8), allocatable :: ext_cff_mss_oc1(:) ! (numrad_snw);

  ! hydrophilic OC
  real(r8), allocatable :: ss_alb_oc2     (:) ! (numrad_snw);
  real(r8), allocatable :: asm_prm_oc2    (:) ! (numrad_snw);
  real(r8), allocatable :: ext_cff_mss_oc2(:) ! (numrad_snw);

  ! dust species 1:
  real(r8), allocatable :: ss_alb_dst1     (:) ! (numrad_snw);
  real(r8), allocatable :: asm_prm_dst1    (:) ! (numrad_snw);
  real(r8), allocatable :: ext_cff_mss_dst1(:) ! (numrad_snw);

  ! dust species 2:
  real(r8), allocatable :: ss_alb_dst2     (:) ! (numrad_snw);
  real(r8), allocatable :: asm_prm_dst2    (:) ! (numrad_snw);
  real(r8), allocatable :: ext_cff_mss_dst2(:) ! (numrad_snw);

  ! dust species 3:
  real(r8), allocatable :: ss_alb_dst3     (:) ! (numrad_snw);
  real(r8), allocatable :: asm_prm_dst3    (:) ! (numrad_snw);
  real(r8), allocatable :: ext_cff_mss_dst3(:) ! (numrad_snw);

  ! dust species 4:
  real(r8), allocatable :: ss_alb_dst4     (:) ! (numrad_snw);
  real(r8), allocatable :: asm_prm_dst4    (:) ! (numrad_snw);
  real(r8), allocatable :: ext_cff_mss_dst4(:) ! (numrad_snw);

#ifdef MODAL_AER
  ! Absorption enhancement factors for within-ice BC
  real(r8), allocatable :: bcenh (:,:,:) ! (numrad_snw,idx_bc_nclrds_max,idx_bcint_icerds_max);
#endif
!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------



  subroutine newsnow_lake ( &
            ! "in" arguments
            ! ---------------
            maxsnl    , nl_lake   , deltim      , dz_lake ,&
            pg_rain   , pg_snow   , t_precip    , bifall  ,&

            ! "inout" arguments
            ! ------------------
            t_lake    , zi_soisno , z_soisno ,&
            dz_soisno , t_soisno  , wliq_soisno , wice_soisno ,&
            fiold     , snl       , sag         , scv         ,&
            snowdp    , lake_icefrac )

!-----------------------------------------------------------------------
! DESCRIPTION:
! Add new snow nodes and interaction btw prec and lake surface including phase change of prec and water body
!
! Original author : Yongjiu Dai, 04/2014
!
! Revisions:
! Nan Wei,  01/2018: update interaction btw prec and lake surface
!-----------------------------------------------------------------------

  use MOD_Precision
  use MOD_Const_Physical, only : tfrz, denh2o, cpliq, cpice, hfus
  implicit none
! ------------------------ Dummy Argument ------------------------------
  integer, INTENT(in) :: maxsnl    ! maximum number of snow layers
  integer, INTENT(in) :: nl_lake   ! number of soil layers
  real(r8), INTENT(in) :: deltim   ! seconds in a time step [second]
  real(r8), INTENT(inout) :: pg_rain  ! liquid water onto ground [kg/(m2 s)]
  real(r8), INTENT(inout) :: pg_snow  ! ice onto ground [kg/(m2 s)]
  real(r8), INTENT(in) :: t_precip ! snowfall/rainfall temperature [kelvin]
  real(r8), INTENT(in) :: bifall   ! bulk density of newly fallen dry snow [kg/m3]

  real(r8), INTENT(in) :: dz_lake(1:nl_lake) ! lake layer thickness (m)
  real(r8), INTENT(inout) ::   zi_soisno(maxsnl:0)   ! interface level below a "z" level (m)
  real(r8), INTENT(inout) ::    z_soisno(maxsnl+1:0) ! snow layer depth (m)
  real(r8), INTENT(inout) ::   dz_soisno(maxsnl+1:0) ! snow layer thickness (m)
  real(r8), INTENT(inout) ::    t_soisno(maxsnl+1:0) ! snow layer temperature [K]
  real(r8), INTENT(inout) :: wliq_soisno(maxsnl+1:0) ! snow layer liquid water (kg/m2)
  real(r8), INTENT(inout) :: wice_soisno(maxsnl+1:0) ! snow layer ice lens (kg/m2)
  real(r8), INTENT(inout) ::       fiold(maxsnl+1:0) ! fraction of ice relative to the total water
   integer, INTENT(inout) :: snl    ! number of snow layers
  real(r8), INTENT(inout) :: sag    ! non dimensional snow age [-]
  real(r8), INTENT(inout) :: scv    ! snow mass (kg/m2)
  real(r8), INTENT(inout) :: snowdp ! snow depth (m)
  real(r8), INTENT(inout) :: lake_icefrac(1:nl_lake) ! mass fraction of lake layer that is frozen
  real(r8), INTENT(inout) :: t_lake(1:nl_lake)       ! lake layer temperature (m)

! ----------------------- Local  Variables -----------------------------

  integer lb
  integer newnode    ! signification when new snow node is set, (1=yes, 0=non)
  real(r8) dz_snowf  ! layer thickness rate change due to precipitation [m/s]
  real(r8) a, b, c, d, e, f, g, h
  real(r8) wice_lake(1:nl_lake), wliq_lake(1:nl_lake), tw

!-----------------------------------------------------------------------

      newnode = 0
      dz_snowf = pg_snow/bifall
      snowdp = snowdp + dz_snowf*deltim
      scv = scv + pg_snow*deltim      ! snow water equivalent (mm)


      zi_soisno(0) = 0.

      IF (snl==0 .and. snowdp < 0.01) then       ! no snow layer, energy exchange between prec and lake surface

          a = cpliq*pg_rain*deltim*(t_precip-tfrz)                          !cool down rainfall to tfrz
          b = pg_rain*deltim*hfus                                           !all rainfall frozen
          c = cpice*denh2o*dz_lake(1)*lake_icefrac(1)*(tfrz-t_lake(1))      !warm up lake surface ice to tfrz
          d = denh2o*dz_lake(1)*lake_icefrac(1)*hfus                        !all lake surface ice melt
          e = cpice*pg_snow*deltim*(tfrz-t_precip)                          !warm up snowfall to tfrz
          f = pg_snow*deltim*hfus                                           !all snowfall melt
          g = cpliq*denh2o*dz_lake(1)*(1-lake_icefrac(1))*(t_lake(1)-tfrz)  !cool down lake surface water to tfrz
          h = denh2o*dz_lake(1)*(1-lake_icefrac(1))*hfus                    !all lake surface water frozen
          sag = 0.0

          if (lake_icefrac(1) > 0.999) then
              ! all rainfall frozen, release heat to warm up frozen lake surface
              if (a+b<=c) then
                  tw=min(tfrz,t_precip)
                  t_lake(1)=(a+b+cpice*(pg_rain+pg_snow)*deltim*tw+cpice*denh2o*dz_lake(1)*t_lake(1)*lake_icefrac(1))/&
                            (cpice*denh2o*dz_lake(1)*lake_icefrac(1)+cpice*(pg_rain+pg_snow)*deltim)
                  scv = scv+pg_rain*deltim
                  snowdp = snowdp + pg_rain*deltim/bifall
                  pg_snow = pg_snow+pg_rain
                  pg_rain = 0.0
              ! prec tem at tfrz, partial rainfall frozen ->release heat -> warm up lake surface to tfrz (no latent heat)
              else if (a<=c) then
                  t_lake(1)=tfrz
                  scv = scv + (c-a)/hfus
                  snowdp = snowdp + (c-a)/(hfus*bifall)
                  pg_snow = pg_snow + min(pg_rain,(c-a)/(hfus*deltim))
                  pg_rain = max(0.0,pg_rain - (c-a)/(hfus*deltim))
              ! lake surface tem at tfrz, partial lake surface melt -> absorb heat -> cool down rainfall to tfrz (no latent heat)
              else if (a<=c+d) then
                  t_lake(1)=tfrz
                  wice_lake(1) = denh2o*dz_lake(1) - (a-c)/hfus
                  wliq_lake(1) = (a-c)/hfus
                  lake_icefrac(1) = wice_lake(1)/(wice_lake(1) + wliq_lake(1))
              ! all lake surface melt, absorb heat to cool down rainfall
              else  !(a>c+d)
                  t_lake(1)=(cpliq*pg_rain*deltim*t_precip+cpliq*denh2o*dz_lake(1)*tfrz-c-d)/&
                            (cpliq*denh2o*dz_lake(1)+cpliq*pg_rain*deltim)
                  lake_icefrac(1) = 0.0
              end if

              if (snowdp>=0.01) then  !frozen rain may make new snow layer
                  snl = -1
                  newnode = 1
                  dz_soisno(0)  = snowdp             ! meter
                  z_soisno (0)  = -0.5*dz_soisno(0)
                  zi_soisno(-1) = -dz_soisno(0)
                  sag = 0.                           ! snow age

                  t_soisno (0) = t_lake(1)           ! K
                  wice_soisno(0) = scv               ! kg/m2
                  wliq_soisno(0) = 0.                ! kg/m2
                  fiold(0) = 1.
              end if

          else if (lake_icefrac(1) >= 0.001) then
              if (pg_rain > 0.0 .and. pg_snow > 0.0) then
                  t_lake(1)=tfrz
              else if (pg_rain > 0.0) then
                  if (a>=d) then
                      t_lake(1)=(cpliq*pg_rain*deltim*t_precip+cpliq*denh2o*dz_lake(1)*tfrz-d)/&
                                (cpliq*denh2o*dz_lake(1)+cpliq*pg_rain*deltim)
                      lake_icefrac(1) = 0.0
                  else
                      t_lake(1)=tfrz
                      wice_lake(1) = denh2o*dz_lake(1)*lake_icefrac(1) - a/hfus
                      wliq_lake(1) = denh2o*dz_lake(1)*(1-lake_icefrac(1)) + a/hfus
                      lake_icefrac(1) = wice_lake(1)/(wice_lake(1) + wliq_lake(1))
                  end if
              else if (pg_snow > 0.0) then
                  if (e>=h) then
                      t_lake(1)=(h+cpice*denh2o*dz_lake(1)*tfrz+cpice*pg_snow*deltim*t_precip)/&
                                (cpice*pg_snow*deltim+cpice*denh2o*dz_lake(1))
                      lake_icefrac(1) = 1.0
                  else
                      t_lake(1)=tfrz
                      wice_lake(1) = denh2o*dz_lake(1)*lake_icefrac(1) + e/hfus
                      wliq_lake(1) = denh2o*dz_lake(1)*(1-lake_icefrac(1)) - e/hfus
                      lake_icefrac(1) = wice_lake(1)/(wice_lake(1) + wliq_lake(1))
                  end if
              end if

          else
              ! all snowfall melt, absorb heat to cool down lake surface water
              if (e+f<=g) then
                  tw=max(tfrz,t_precip)
                  t_lake(1)=(cpliq*denh2o*dz_lake(1)*t_lake(1)*(1-lake_icefrac(1))+cpliq*(pg_rain+pg_snow)*deltim*tw-e-f)/&
                            (cpliq*(pg_rain+pg_snow)*deltim+cpliq*denh2o*dz_lake(1)*(1-lake_icefrac(1)))
                  scv = scv - pg_snow*deltim
                  snowdp = snowdp - dz_snowf*deltim
                  pg_rain = pg_rain + pg_snow
                  pg_snow = 0.0
              ! prec tem at tfrz, partial snowfall melt ->absorb heat -> cool down lake surface to tfrz (no latent heat)
              else if (e<=g) then
                  t_lake(1) = tfrz
                  scv = scv - (g-e)/hfus
                  snowdp = snowdp - (g-e)/(hfus*bifall)
                  pg_rain = pg_rain + min(pg_snow, (g-e)/(hfus*deltim))
                  pg_snow = max(0.0, pg_snow - (g-e)/(hfus*deltim))
              ! lake surface tem at tfrz, partial lake surface frozen -> release heat -> warm up snowfall to tfrz (no latent heat)
              else if (e<=g+h) then
                  t_lake(1) = tfrz
                  wice_lake(1) = (e-g)/hfus
                  wliq_lake(1) = denh2o*dz_lake(1) - (e-g)/hfus
                  lake_icefrac(1) = wice_lake(1)/(wice_lake(1) + wliq_lake(1))
              ! all lake surface frozen, release heat to warm up snowfall
              else       !(e>g+h)
                  t_lake(1) = (g+h+cpice*denh2o*dz_lake(1)*tfrz+cpice*pg_snow*deltim*t_precip)/&
                              (cpice*pg_snow*deltim+cpice*denh2o*dz_lake(1))
                  lake_icefrac(1) = 1.0
              end if
          end if

      ELSE IF (snl==0 .and. snowdp >= 0.01) then

          ! only ice part of snowfall is added here, the liquid part will be added later
          snl = -1
          newnode = 1
          dz_soisno(0)  = snowdp             ! meter
          z_soisno (0)  = -0.5*dz_soisno(0)
          zi_soisno(-1) = -dz_soisno(0)
          sag = 0.                           ! snow age

          t_soisno (0) = min(tfrz, t_precip) ! K
          wice_soisno(0) = scv               ! kg/m2
          wliq_soisno(0) = 0.                ! kg/m2
          fiold(0) = 1.

      ELSE                                   ! ( snl<0 .and. newnode ==0 )

          lb = snl + 1
          t_soisno(lb) = ( (wice_soisno(lb)*cpice+wliq_soisno(lb)*cpliq)*t_soisno(lb) &
                       +   (pg_rain*cpliq + pg_snow*cpice)*deltim*t_precip ) &
                       / ( wice_soisno(lb)*cpice + wliq_soisno(lb)*cpliq &
                       +   pg_rain*deltim*cpliq + pg_snow*deltim*cpice )

          t_soisno(lb) = min(tfrz, t_soisno(lb))
          wice_soisno(lb) = wice_soisno(lb)+deltim*pg_snow
          dz_soisno(lb) = dz_soisno(lb)+dz_snowf*deltim
          z_soisno(lb) = zi_soisno(lb) - 0.5*dz_soisno(lb)
          zi_soisno(lb-1) = zi_soisno(lb) - dz_soisno(lb)

      END IF

  end subroutine newsnow_lake



  subroutine laketem (&
           ! "in" arguments
           ! -------------------
           patchtype    , maxsnl      , nl_soil      , nl_lake   ,&
           dlat         , deltim      , forc_hgt_u   , forc_hgt_t,&
           forc_hgt_q   , forc_us     , forc_vs      , forc_t    ,&
           forc_q       , forc_rhoair , forc_psrf    , forc_sols ,&
           forc_soll    , forc_solsd  , forc_solld   , sabg      ,&
           forc_frl     , dz_soisno   , z_soisno     , zi_soisno ,&
           dz_lake      , lakedepth   , vf_quartz    , vf_gravels,&
           vf_om        , vf_sand     , wf_gravels   , wf_sand   ,&
           porsl        , csol        , k_solids     , &
           dksatu       , dksatf      , dkdry        , &
           BA_alpha     , BA_beta     , hpbl         , &
           dlon         , idate       , sabgv        , &
           sabgvd      , sabgvda      , sabgvdu      , sabgn     ,&
           sabgnd       ,&
           h2osno, frac_sno, h2osno_liq, h2osno_ice, snw_rds,&
           albsfc, &
           mss_bcpho    ,mss_bcphi    ,mss_ocpho       ,mss_ocphi       ,&
           mss_dst1     ,mss_dst2     ,mss_dst3        ,mss_dst4,         &   

           ! "inout" arguments
           ! -------------------
           t_grnd       , scv         , snowdp       , t_soisno  ,&
           wliq_soisno  , wice_soisno , imelt_soisno , t_lake    ,&
           lake_icefrac , savedtke1, &

! SNICAR model variables
           snofrz       ,sabg_lyr     ,&
! END SNICAR model variables

           ! "out" arguments
           ! -------------------
           taux         , tauy        , fsena                    ,&
           fevpa        , lfevpa      , fseng        , fevpg     ,&
           qseva        , qsubl       , qsdew        , qfros     ,&
           olrg         , fgrnd       , tref         , qref      ,&
           trad         , emis        , z0m          , zol       ,&
           rib          , ustar       , qstar        , tstar     ,&
           fm           , fh          , fq           , sm        ,&
           urban_call)

! ------------------------ code history ---------------------------
! purpose: lake temperature and snow on frozen lake
! initial  Yongjiu Dai, 2000
!          Zack Subin, 2009
!          Yongjiu Dai, /12/2012/, /04/2014/, 06/2018
!          Nan Wei, /06/2018/
!
! ------------------------ notes ----------------------------------
! Lakes have variable depth, possible snow layers above, freezing & thawing of lake water,
! and soil layers with active temperature and gas diffusion below.
!
! Calculates temperatures in the 25-30 layer column of (possible) snow,
! lake water, soil, and bedrock beneath lake.
! Snow and soil temperatures are determined as in SoilTemperature, except
! for appropriate boundary conditions at the top of the snow (the flux is fixed
! to be the ground heat flux), the bottom of the snow (adjacent to top lake layer),
! and the top of the soil (adjacent to the bottom lake layer).
! Also, the soil is kept fully saturated.
! The whole column is solved simultaneously as one tridiagonal matrix.
!
! calculate lake temperatures from one-dimensional thermal
! stratification model based on eddy diffusion concepts to
! represent vertical mixing of heat
!
! d ts    d            d ts     1 ds
! ---- = -- [(km + ke) ----] + -- --
!  dt    dz             dz     cw dz
! where: ts = temperature (kelvin)
!         t = time (s)
!         z = depth (m)
!        km = molecular diffusion coefficient (m**2/s)
!        ke = eddy diffusion coefficient (m**2/s)
!        cw = heat capacity (j/m**3/kelvin)
!         s = heat source term (w/m**2)
!
! use crank-nicholson method to set up tridiagonal system of equations to
! solve for ts at time n+1, where the temperature equation for layer i is
! r_i = a_i [ts_i-1] n+1 + b_i [ts_i] n+1 + c_i [ts_i+1] n+1
! the solution conserves energy as
! cw*([ts(  1)] n+1 - [ts(  1)] n)*dz(  1)/dt + ... +
! cw*([ts(nl_lake)] n+1 - [ts(nl_lake)] n)*dz(nl_lake)/dt = fin
! where
! [ts] n   = old temperature (kelvin)
! [ts] n+1 = new temperature (kelvin)
! fin      = heat flux into lake (w/m**2)
!          = beta*sabg_lyr(1)+forc_frl-olrg-fsena-lfevpa-hm + phi(1) + ... + phi(nl_lake)
!
! REVISIONS:
! Yongjiu Dai and Hua Yuan, 01/2023: added SNICAR for layer solar absorption, ground heat
!                                    flux, temperature and freezing mass calculations
! Shaofeng Liu, 05/2023: add option to call moninobuk_leddy, the LargeEddy
!                        surface turbulence scheme (LZD2022);
!                        make a proper update of um.
!
! -----------------------------------------------------------------
  use MOD_Precision
  use MOD_Const_Physical, only : tfrz,hvap,hfus,hsub,tkwat,tkice,tkair,stefnc,&
                                  vonkar,grav,cpliq,cpice,cpair,denh2o,denice,rgas
  use MOD_FrictionVelocity
  USE MOD_Namelist, only: DEF_USE_CBL_HEIGHT, DEF_USE_SNICAR, DEF_USE_ORIGINAL, DEF_USE_ONEBANDICE, DEF_USE_TWOBANDICE, DEF_USE_TWOSTREAM
  USE MOD_TurbulenceLEddy
  USE MOD_Qsadv
  USE MOD_SoilThermalParameters
  USE MOD_Utils
  USE MOD_OrbCoszen
  use MOD_TimeManager
!   use MOD_Vars_Global, only: maxsnl
  use MOD_NetCDFSerial
  USE MOD_Aerosol, only: AerosolMasses
  IMPLICIT NONE
! ------------------------ input/output variables -----------------

  integer, INTENT(in) :: patchtype    ! land patch type (4=deep lake, 5=shallow lake)
  integer, INTENT(in) :: maxsnl       ! maximum number of snow layers
  integer, INTENT(in) :: nl_soil      ! number of soil layers
  integer, INTENT(in) :: nl_lake      ! number of lake layers
  integer, INTENT(in) :: idate(3)      ! next time-step /year/julian day/second in a day/
  real(r8), INTENT(in) :: dlat        ! latitude (radians)
  real(r8), INTENT(in) :: dlon        ! longtitude (radians)
  real(r8), INTENT(in) :: deltim      ! seconds in a time step (s)
  real(r8), INTENT(in) :: forc_hgt_u  ! observational height of wind [m]
  real(r8), INTENT(in) :: forc_hgt_t  ! observational height of temperature [m]
  real(r8), INTENT(in) :: forc_hgt_q  ! observational height of humidity [m]
  real(r8), INTENT(in) :: forc_us     ! wind component in eastward direction [m/s]
  real(r8), INTENT(in) :: forc_vs     ! wind component in northward direction [m/s]
  real(r8), INTENT(in) :: forc_t      ! temperature at agcm reference height [kelvin]
  real(r8), INTENT(in) :: forc_q      ! specific humidity at agcm reference height [kg/kg]
  real(r8), INTENT(in) :: forc_rhoair ! density air [kg/m3]
  real(r8), INTENT(in) :: forc_psrf   ! atmosphere pressure at the surface [pa]
  real(r8), INTENT(in) :: forc_sols   ! atm vis direct beam solar rad onto srf [W/m2]
  real(r8), INTENT(in) :: forc_soll   ! atm nir direct beam solar rad onto srf [W/m2]
  real(r8), INTENT(in) :: forc_solsd  ! atm vis diffuse solar rad onto srf [W/m2]
  real(r8), INTENT(in) :: forc_solld  ! atm nir diffuse solar rad onto srf [W/m2]
  real(r8), INTENT(in) :: forc_frl    ! atmospheric infrared (longwave) radiation [W/m2]
  real(r8), INTENT(in) :: sabg        ! solar radiation absorbed by ground [W/m2]
  real(r8), INTENT(in) :: sabgv       ! direct beam vis solar absorbed by ground  [W/m2] 
  real(r8), INTENT(in) :: sabgvd      ! diffuse beam vis solar absorbed by ground  [W/m2]
  real(r8), INTENT(in) :: sabgvdu     ! diffuse beam vis solar reflected by ground  [W/m2]
  real(r8), INTENT(in) :: sabgvda     ! diffuse beam vis solar reflected and absorbed by ground  [W/m2]
  real(r8), INTENT(in) :: sabgn       ! direct beam nir solar absorbed by ground  [W/m2]
  real(r8), INTENT(in) :: sabgnd      ! diffuse beam nir solar absorbed by ground  [W/m2] 

  ! ------------------------SNICAR-----------------------------------------
  integer   :: flg_snw_ice                  ! flag: =1 when called from CLM, =2 when called from CSIM
  integer  :: flg_slr_in                   ! flag: =1 for direct-beam incident flux,=2 for diffuse incident flux
!   integer, INTENT(in) :: maxsnl       ! maximum number of snow layers
  real(r8) , intent(in)  :: h2osno          !  snow liquid water equivalent (col) [kg/m2]
  real(r8) , intent(in)  :: frac_sno        !  fraction of ground covered by snow (0 to 1)
  real(r8) , intent(in)  :: h2osno_liq     ( maxsnl+1:0 )   ! liquid water content (col,lyr) [kg/m2]
  real(r8) , intent(in)  :: h2osno_ice     ( maxsnl+1:0 )   ! ice content (col,lyr) [kg/m2]
  real(r8) , intent(inout)  :: snw_rds        ( maxsnl+1: 0)   ! snow effective radius (col,lyr) [microns, m^-6]
  real(r8) , intent(in)  :: albsfc         ( 1:numrad )     ! albedo of surface underlying snow (col,bnd) [frc]
  real(r8) :: flx_abs        ( maxsnl+1:12 , 1:numrad)! absorbed flux in each layer per unit flux incident (col, lyr, bnd)
! -----------------------------------------------------------------------  

  real(r8), INTENT(in) :: dz_soisno(maxsnl+1:nl_soil) ! soil/snow layer thickness (m)
  real(r8), INTENT(in) :: z_soisno(maxsnl+1:nl_soil)  ! soil/snow node depth [m]
  real(r8), INTENT(in) :: zi_soisno(maxsnl:nl_soil)   ! soil/snow depth of layer interface [m]

  real(r8), INTENT(in) :: dz_lake(nl_lake)      ! lake layer thickness (m)
  real(r8), INTENT(in) :: lakedepth             ! column lake depth (m)

  real(r8), INTENT(in) :: vf_quartz (1:nl_soil) ! volumetric fraction of quartz within mineral soil
  real(r8), INTENT(in) :: vf_gravels(1:nl_soil) ! volumetric fraction of gravels
  real(r8), INTENT(in) :: vf_om     (1:nl_soil) ! volumetric fraction of organic matter
  real(r8), INTENT(in) :: vf_sand   (1:nl_soil) ! volumetric fraction of sand
  real(r8), INTENT(in) :: wf_gravels(1:nl_soil) ! gravimetric fraction of gravels
  real(r8), INTENT(in) :: wf_sand   (1:nl_soil) ! gravimetric fraction of sand
  real(r8), INTENT(in) :: porsl(1:nl_soil)      ! soil porosity [-]

  real(r8), INTENT(in) :: csol(1:nl_soil)       ! heat capacity of soil solids [J/(m3 K)]
  real(r8), INTENT(in) :: k_solids(1:nl_soil)   ! thermal conductivity of mineralssoil [W/m-K]
  real(r8), INTENT(in) :: dksatu(1:nl_soil)     ! thermal conductivity of saturated unfrozen soil [W/m-K]
  real(r8), INTENT(in) :: dksatf(1:nl_soil)     ! thermal conductivity of saturated frozen soil [W/m-K]
  real(r8), INTENT(in) :: dkdry(1:nl_soil)      ! thermal conductivity of dry soil [W/m-K]
  real(r8), INTENT(in) :: BA_alpha(1:nl_soil)   ! alpha in Balland and Arp(2005) thermal conductivity scheme
  real(r8), INTENT(in) :: BA_beta(1:nl_soil)    ! beta in Balland and Arp(2005) thermal conductivity scheme
  real(r8), INTENT(in) :: hpbl                  ! atmospheric boundary layer height [m]

  real(r8), INTENT(inout) :: t_grnd             ! surface temperature (kelvin)
  real(r8), INTENT(inout) :: scv                ! snow water equivalent [mm]
  real(r8), INTENT(inout) :: snowdp             ! snow depth [mm]

  real(r8), INTENT(inout) :: t_soisno    (maxsnl+1:nl_soil) ! soil/snow temperature [K]
  real(r8), INTENT(inout) :: wliq_soisno (maxsnl+1:nl_soil) ! soil/snow liquid water (kg/m2)
  real(r8), INTENT(inout) :: wice_soisno (maxsnl+1:nl_soil) ! soil/snow ice lens (kg/m2)
  integer,  INTENT(inout) :: imelt_soisno(maxsnl+1:nl_soil) ! soil/snow flag for melting (=1), freezing (=2), Not=0 (new)

  real(r8), INTENT(inout) :: t_lake(nl_lake)       ! lake temperature (kelvin)
  real(r8), INTENT(inout) :: lake_icefrac(nl_lake) ! lake mass fraction of lake layer that is frozen
  real(r8), INTENT(inout) :: savedtke1             ! top level eddy conductivity (W/m K)

! SNICAR model variables
  REAL(r8), intent(out) :: snofrz   (maxsnl+1:0)   ! snow freezing rate (col,lyr) [kg m-2 s-1]
  REAL(r8), intent(in)  :: sabg_lyr (maxsnl+1:1)   ! solar radiation absorbed by ground [W/m2]
! END SNICAR model variables

  real(r8), INTENT(out) :: taux   ! wind stress: E-W [kg/m/s**2]
  real(r8), INTENT(out) :: tauy   ! wind stress: N-S [kg/m/s**2]
  real(r8), INTENT(out) :: fsena  ! sensible heat from canopy height to atmosphere [W/m2]
  real(r8), INTENT(out) :: fevpa  ! evapotranspiration from canopy height to atmosphere [mm/s]
  real(r8), INTENT(out) :: lfevpa ! latent heat flux from canopy height to atmosphere [W/m2]

  real(r8), INTENT(out) :: fseng  ! sensible heat flux from ground [W/m2]
  real(r8), INTENT(out) :: fevpg  ! evaporation heat flux from ground [mm/s]

  real(r8), INTENT(out) :: qseva  ! ground surface evaporation rate (mm h2o/s)
  real(r8), INTENT(out) :: qsubl  ! sublimation rate from snow pack (mm H2O /s) [+]
  real(r8), INTENT(out) :: qsdew  ! surface dew added to snow pack (mm H2O /s) [+]
  real(r8), INTENT(out) :: qfros  ! ground surface frosting formation (mm H2O /s) [+]

  real(r8), INTENT(out) :: olrg   ! outgoing long-wave radiation from ground+canopy
  real(r8), INTENT(out) :: fgrnd  ! ground heat flux [W/m2]

  real(r8), INTENT(out) :: tref   ! 2 m height air temperature [kelvin]
  real(r8), INTENT(out) :: qref   ! 2 m height air specific humidity
  real(r8), INTENT(out) :: trad   ! radiative temperature [K]

  real(r8), INTENT(out) :: emis   ! averaged bulk surface emissivity
  real(r8), INTENT(out) :: z0m    ! effective roughness [m]
  real(r8), INTENT(out) :: zol    ! dimensionless height (z/L) used in Monin-Obukhov theory
  real(r8), INTENT(out) :: rib    ! bulk Richardson number in surface layer
  real(r8), INTENT(out) :: ustar  ! u* in similarity theory [m/s]
  real(r8), INTENT(out) :: qstar  ! q* in similarity theory [kg/kg]
  real(r8), INTENT(out) :: tstar  ! t* in similarity theory [K]
  real(r8), INTENT(out) :: fm     ! integral of profile function for momentum
  real(r8), INTENT(out) :: fh     ! integral of profile function for heat
  real(r8), INTENT(out) :: fq     ! integral of profile function for moisture
  real(r8), INTENT(out) :: sm     ! rate of snowmelt [mm/s, kg/(m2 s)]
  logical, optional, intent(in) :: urban_call   ! whether it is a urban CALL

! ---------------- local variables in surface temp and fluxes calculation -----------------
  integer idlak     ! index of lake, 1 = deep lake, 2 = shallow lake
  real(r8) z_lake (nl_lake)  ! lake node depth (middle point of layer) (m)

  real(r8) ax       ! used in iteration loop for calculating t_grnd (numerator of NR solution)
  real(r8) bx       ! used in iteration loop for calculating t_grnd (denomin. of NR solution)
  real(r8) beta1    ! coefficient of conective velocity [-]
  real(r8) degdT    ! d(eg)/dT
  real(r8) displax  ! zero- displacement height [m]
  real(r8) dqh      ! diff of humidity between ref. height and surface
  real(r8) dth      ! diff of virtual temp. between ref. height and surface
  real(r8) dthv     ! diff of vir. poten. temp. between ref. height and surface
  real(r8) dzsur    ! 1/2 the top layer thickness (m)
  real(r8) tsur     ! top layer temperature
  real(r8) rhosnow  ! partitial density of water (ice + liquid)
  real(r8) eg       ! water vapor pressure at temperature T [pa]
  real(r8) emg      ! ground emissivity (0.97 for snow,
  real(r8) errore   ! lake temperature energy conservation error (w/m**2)
  real(r8) hm       ! energy residual [W/m2]
  real(r8) htvp     ! latent heat of vapor of water (or sublimation) [j/kg]
  real(r8) obu      ! monin-obukhov length (m)
  real(r8) obuold   ! monin-obukhov length of previous iteration
  real(r8) qsatg    ! saturated humidity [kg/kg]
  real(r8) qsatgdT  ! d(qsatg)/dT

  real(r8) ram      ! aerodynamical resistance [s/m]
  real(r8) rah      ! thermal resistance [s/m]
  real(r8) raw      ! moisture resistance [s/m]
  real(r8) stftg3   ! emg*sb*t_grnd*t_grnd*t_grnd
  real(r8) fh2m     ! relation for temperature at 2m
  real(r8) fq2m     ! relation for specific humidity at 2m
  real(r8) fm10m    ! integral of profile function for momentum at 10m
  real(r8) t_grnd_bef0   ! initial ground temperature
  real(r8) t_grnd_bef    ! initial ground temperature
  real(r8) thm      ! intermediate variable (forc_t+0.0098*forc_hgt_t)
  real(r8) th       ! potential temperature (kelvin)
  real(r8) thv      ! virtual potential temperature (kelvin)
  real(r8) thvstar  ! virtual potential temperature scaling parameter
  real(r8) tksur    ! thermal conductivity of snow/soil (w/m/kelvin)
  real(r8) um       ! wind speed including the stablity effect [m/s]
  real(r8) ur       ! wind speed at reference height [m/s]
  real(r8) visa     ! kinematic viscosity of dry air [m2/s]
  real(r8) wc       ! convective velocity [m/s]
  real(r8) wc2      ! wc*wc
  real(r8) zeta     ! dimensionless height used in Monin-Obukhov theory
  real(r8) zii      ! convective boundary height [m]
  real(r8) zldis    ! reference height "minus" zero displacement heght [m]
  real(r8) z0mg     ! roughness length over ground, momentum [m]
  real(r8) z0hg     ! roughness length over ground, sensible heat [m]
  real(r8) z0qg     ! roughness length over ground, latent heat [m]

  real(r8) wliq_lake(nl_lake)  ! lake liquid water (kg/m2)
  real(r8) wice_lake(nl_lake)  ! lake ice lens (kg/m2)
  real(r8) vf_water(1:nl_soil) ! volumetric fraction liquid water within underlying soil
  real(r8) vf_ice(1:nl_soil)   ! volumetric fraction ice len within underlying soil

  real(r8) fgrnd1  ! ground heat flux into the first snow/lake layer [W/m2]

! ---------------- local variables in lake/snow/soil temperature calculation --------------
  real(r8), parameter :: cur0 = 0.01     ! min. Charnock parameter
  real(r8), parameter :: curm = 0.1      ! maximum Charnock parameter
  real(r8), parameter :: fcrit = 22.     ! critical dimensionless fetch for Charnock parameter (Vickers & Mahrt 1997)
                                         ! but converted to use u instead of u* (Subin et al. 2011)
  real(r8), parameter :: mixfact = 5.    ! Mixing enhancement factor.
  real(r8), parameter :: depthcrit = 25. ! (m) Depth beneath which to enhance mixing
  real(r8), parameter :: fangmult = 5.   ! Multiplier for unfrozen diffusivity
  real(r8), parameter :: minmultdepth = 20. ! (m) Minimum depth for imposing fangmult
  real(r8), parameter :: cnfac  = 0.5    ! Crank Nicholson factor between 0 and 1

  !--------------------
  real(r8) fetch     ! lake fetch (m)
  real(r8) cur       ! Charnock parameter (-)
  real(r8) betavis   !
  real(r8) betaprime ! Effective beta
  real(r8) tdmax     ! temperature of maximum water density
  real(r8) cfus      ! effective heat of fusion per unit volume
  real(r8) tkice_eff ! effective conductivity since layer depth is constant
  real(r8) cice_eff  ! effective heat capacity of ice (using density of
                     ! water because layer depth is not adjusted when freezing
  real(r8) cwat      ! specific heat capacity of water (j/m**3/kelvin)

  !--------------------
  real(r8) rhow(nl_lake) ! density of water (kg/m**3)
  real(r8) fin           ! heat flux into lake - flux out of lake (w/m**2)
  real(r8) phi(maxsnl+1:nl_lake)  ! solar radiation absorbed by layer (w/m**2)
  real(r8) phi_soil      ! solar radiation into top soil layer (W/m^2)
  real(r8) phidum        ! temporary value of phi

  integer  imelt_lake(1:nl_lake)       ! lake flag for melting or freezing snow and soil layer [-]
  real(r8) cv_lake(1:nl_lake)          ! heat capacity [J/(m2 K)]
  real(r8) tk_lake(1:nl_lake)          ! thermal conductivity at layer node [W/(m K)]
  real(r8) cv_soisno(maxsnl+1:nl_soil) ! heat capacity of soil/snow [J/(m2 K)]
  real(r8) tk_soisno(maxsnl+1:nl_soil) ! thermal conductivity of soil/snow [W/(m K)] (at interface below, except for j=0)
  real(r8) hcap(1:nl_soil)           ! J/(m3 K)
  real(r8) thk(maxsnl+1:nl_soil)     ! W/(m K)
  real(r8) tktopsoil                   ! thermal conductivity of the top soil layer [W/(m K)]

  real(r8) t_soisno_bef(maxsnl+1:nl_soil) ! beginning soil/snow temp for E cons. check [K]
  real(r8) t_lake_bef(1:nl_lake)          ! beginning lake temp for energy conservation check [K]
  real(r8) wice_soisno_bef(maxsnl+1:0)    ! ice lens [kg/m2]

  real(r8) cvx    (maxsnl+1:nl_lake+nl_soil) ! heat capacity for whole column [J/(m2 K)]
  real(r8) tkix   (maxsnl+1:nl_lake+nl_soil) ! thermal conductivity at layer interfaces for whole column [W/(m K)]
  real(r8) phix   (maxsnl+1:nl_lake+nl_soil) ! solar source term for whole column [W/m**2]
  real(r8) zx     (maxsnl+1:nl_lake+nl_soil) ! interface depth (+ below surface) for whole column [m]
  real(r8) tx     (maxsnl+1:nl_lake+nl_soil) ! temperature of whole column [K]
  real(r8) tx_bef (maxsnl+1:nl_lake+nl_soil) ! beginning lake/snow/soil temp for energy conservation check [K]
  real(r8) factx  (maxsnl+1:nl_lake+nl_soil) ! coefficient used in computing tridiagonal matrix
  real(r8) fnx    (maxsnl+1:nl_lake+nl_soil) ! heat diffusion through the layer interface below [W/m2]
  real(r8) a      (maxsnl+1:nl_lake+nl_soil) ! "a" vector for tridiagonal matrix
  real(r8) b      (maxsnl+1:nl_lake+nl_soil) ! "b" vector for tridiagonal matrix
  real(r8) c      (maxsnl+1:nl_lake+nl_soil) ! "c" vector for tridiagonal matrix
  real(r8) r      (maxsnl+1:nl_lake+nl_soil) ! "r" vector for tridiagonal solution
  real(r8) fn1    (maxsnl+1:nl_lake+nl_soil) ! heat diffusion through the layer interface below [W/m2]
  real(r8) brr    (maxsnl+1:nl_lake+nl_soil) !
  integer  imelt_x(maxsnl+1:nl_lake+nl_soil) ! flag for melting (=1), freezing (=2), Not=0 (new)

  real(r8) dzm       ! used in computing tridiagonal matrix [m]
  real(r8) dzp       ! used in computing tridiagonal matrix [m]
  real(r8) zin       ! depth at top of layer (m)
  real(r8) zout      ! depth at bottom of layer (m)
  real(r8) rsfinv    ! relative flux of ice visible solar radiation into layer
  real(r8) rsfoutv   ! relative flux of ice visible solar radiation out of layer
  real(r8) rsfinn    ! relative flux of ice near-infrared solar radiation into layer
  real(r8) rsfoutn   ! relative flux of ice near-infrared solar radiation out of layer
  real(r8) rsfin     ! relative flux of solar radiation into layer
  real(r8) rsfout    ! relative flux of solar radiation out of layer
  real(r8) eta       ! light extinction coefficient (/m): depends on lake type
  real(r8) etaiv     ! light extinction coefficient (/m): for ice,visible
  real(r8) etain     ! light extinction coefficient (/m): for ice,near-infrared
  real(r8) etav 
  real(r8) etan 
  real(r8) rsfinvw 
  real(r8) rsfinnw   
  real(r8) rsfoutvw 
  real(r8) rsfoutnw 
  real(r8) etatsa(nl_lake)       ! light extinction coefficient (/m): depends on lake type
  real(r8) topc 
  real(r8) daT  
  real(r8) rou 
  real(r8) D 
  real(r8) w0
  real(r8) uavg 
!   real(r8) pi 
  real(r8) sum
  real(r8) dtao(nl_lake)
  real(r8) tao(nl_lake+1)
  real(r8) eup(nl_lake+1)
  real(r8) eupd(nl_lake+1)
  real(r8) eupz(nl_lake+1)
  real(r8) ssk
  real(r8) ssh
  real(r8) edown(nl_lake)
!   real(r8) g         !asymmetry factor
  real(r8) f
  real(r8) s1
  real(r8) s2 
  real(r8) k 
  real(r8) ssa 
  real(r8) z1
  real(r8) z2
  real(r8) saa
  real(r8) sbe
  real(r8) se
  real(r8) sr
  real(r8) su
  real(r8) sv
  real(r8) se2
  real(r8) sr2
  real(r8) su2
  real(r8) sv2
  real(r8) k2
  real(r8) sk2
  real(r8) sh2
!   real(r8) r1
!   real(r8) r2
  real(r8) r3
  real(r8) tu
  real(r8) ttu
  real(r8) coszen
  real(r8) czen
  real(r8) calday
  real(r8) sum2
  real(r8) za(2)     ! base of surface absorption layer (m): depends on lake type
  !-----------------SNICAR--------------------------------------------
  integer :: snl_lcl                                   ! negative number of snow layers [nbr]
  integer :: snw_rds_lcl(maxsnl+1:11)                   ! snow effective radius [m^-6]
  integer :: snw_rds_idx(maxsnl+1:11)
  real(r8):: dif_a(1:numrad_snw)                      ! +CAW
  real(r8):: dif_b(1:numrad_snw)                      ! +CAW
  real(r8):: flx_abs_btm(1:numrad_snw)
  real(r8):: phifsum
  real(r8):: phi_soilfsum
  real(r8):: srf_lyr
  real(r8):: phi_srf
  real(r8):: srf_prop
  real(r8) :: dziw(1:11)     !(1:nbr_lyr) The thickness of the ice and water layers
  real(r8) :: ziw (1:11)
  real(r8) :: rho_snw(1:11)
  real(r8) :: phi_up
  real(r8):: phif(2,maxsnl+1:nl_lake)
  real(r8):: phiup(maxsnl+1:0)
  real(r8):: phidn(2,1:numrad_snw)
  real(r8):: flx_abs_dn_lcl(1:numrad_snw)
  real(r8):: flx_abs_up_lcl(1:numrad_snw)
  real(r8):: flx_abs_dn(1:numrad)
  real(r8):: flx_abs_up(1:numrad)
  real(r8):: flx_abs_dnd(1:numrad)
  real(r8):: flx_abs_upd(1:numrad)
  real(r8):: flx_abs_dni(1:numrad)
  real(r8):: flx_abs_upi(1:numrad)
  real(r8):: flx_absi(maxsnl+1:12,1:numrad)
  real(r8):: flx_absd(maxsnl+1:12,1:numrad)
  real(r8):: phi_soilf(2)
  real(r8):: flx_slrd_lcl(1:numrad_snw)                ! direct beam incident irradiance [W/m2] (set to 1)
  real(r8):: flx_slri_lcl(1:numrad_snw)                ! diffuse incident irradiance [W/m2] (set to 1)
  real(r8):: mss_cnc_aer_lcl(maxsnl+1:11,1:sno_nbr_aer) ! aerosol mass concentration (lyr,aer_nbr) [kg/kg]
  real(r8):: h2osno_lcl                                ! total column snow mass [kg/m2]
  real(r8):: h2osno_liq_lcl(maxsnl+1:11)                ! liquid water mass [kg/m2]
  real(r8):: h2osno_ice_lcl(maxsnl+1:11)                ! ice mass [kg/m2]
  real(r8):: albsfc_lcl(1:numrad_snw)                  ! albedo of underlying surface [frc]
  real(r8):: ss_alb_snw_lcl(maxsnl+1:11)                ! single-scatter albedo of ice grains (lyr) [frc]
  real(r8):: asm_prm_snw_lcl(maxsnl+1:11)               ! asymmetry parameter of ice grains (lyr) [frc]
  real(r8):: ext_cff_mss_snw_lcl(maxsnl+1:11)           ! mass extinction coefficient of ice grains (lyr) [m2/kg]
  real(r8):: ss_alb_aer_lcl(sno_nbr_aer)               ! single-scatter albedo of aerosol species (aer_nbr) [frc]
  real(r8):: asm_prm_aer_lcl(sno_nbr_aer)              ! asymmetry parameter of aerosol species (aer_nbr) [frc]
  real(r8):: ext_cff_mss_aer_lcl(sno_nbr_aer)          ! mass extinction coefficient of aerosol species (aer_nbr) [m2/kg]
  real(r8) :: rds_bcint_lcl(maxsnl+1:11)       ! effective radius of within-ice BC [nm]
  real(r8) :: rds_bcext_lcl(maxsnl+1:11)       ! effective radius of external BC [nm]
  ! Other local variables
  integer :: DELTA                               ! flag to use Delta approximation (Joseph, 1976)
                                                 ! (1= use, 0= don't use)
  real(r8):: flx_wgt(1:numrad_snw)               ! weights applied to spectral bands,
                                                 ! specific to direct and diffuse cases (bnd) [frc]
  integer :: flg_nosnl                           ! flag: =1 if there is snow, but zero snow layers,
                                                 ! =0 if at least 1 snow layer [flg]
                                                 ! integer :: trip                                                                          ! flag: =1 to redo RT calculation if result is unrealistic
                                                 ! integer :: flg_dover  
  real(r8):: albedo                              ! temporary snow albedo [frc]
  real(r8):: albout         ( 1:numrad ) 
  real(r8):: flx_sum                             ! temporary summation variable for NIR weighting
  real(r8):: albout_lcl(numrad_snw)              ! snow albedo by band [frc]
  real(r8):: flx_abs_lcl(maxsnl+1:12,numrad_snw)  ! absorbed flux per unit incident flux at top of snowpack (lyr,bnd) [frc]
  real(r8):: flx_sum2(maxsnl+1:12)
  real(r8):: L_snw(maxsnl+1:11)                   ! h2o mass (liquid+solid) in snow layer (lyr) [kg/m2]
  real(r8):: tau_snw(maxsnl+1:11)                 ! snow optical depth (lyr) [unitless]
  real(r8):: L_aer(maxsnl+1:11,sno_nbr_aer)       ! aerosol mass in snow layer (lyr,nbr_aer) [kg/m2]
  real(r8):: tau_aer(maxsnl+1:11,sno_nbr_aer)     ! aerosol optical depth (lyr,nbr_aer) [unitless]
  real(r8):: tau_sum                              ! cumulative (snow+aerosol) optical depth [unitless]
  real(r8):: tau_elm(maxsnl+1:11)                 ! column optical depth from layer bottom to snowpack top (lyr) [unitless]
  real(r8):: omega_sum                            ! temporary summation of single-scatter albedo of all aerosols [frc]
  real(r8):: g_sum                                ! temporary summation of asymmetry parameter of all aerosols [frc]

  real(r8):: tau(maxsnl+1:11)                     ! weighted optical depth of snow+aerosol layer (lyr) [unitless]
  real(r8):: omega(maxsnl+1:11)                   ! weighted single-scatter albedo of snow+aerosol layer (lyr) [frc]
  real(r8):: g(maxsnl+1:11)                       ! weighted asymmetry parameter of snow+aerosol layer (lyr) [frc]
  real(r8):: tau_star(maxsnl+1:11)                ! transformed (i.e. Delta-Eddington) optical depth of snow+aerosol layer
                                                  ! (lyr) [unitless]
  real(r8):: omega_star(maxsnl+1:11)              ! transformed (i.e. Delta-Eddington) SSA of snow+aerosol layer (lyr) [frc]
  real(r8):: g_star(maxsnl+1:11)                  ! transformed (i.e. Delta-Eddington) asymmetry paramater of snow+aerosol layer
                                                  ! (lyr) [frc]    
  integer :: bnd_idx                             ! spectral band index (1 <= bnd_idx <= numrad_snw) [idx]
  integer :: rds_idx                             ! snow effective radius index for retrieving
  integer :: snl_btm                             ! index of bottom snow layer (0) [idx]
  integer :: snl_top                             ! index of top snow layer (-4 to 0) [idx]
  integer :: fc                                  ! column filter index
  integer :: m                                   ! secondary layer index [idx]
  integer :: nint_snw_rds_min                    ! nearest integer value of snw_rds_min

  real(r8):: F_abs(maxsnl+1:11)                   ! net absorbed radiative energy (lyr) [W/m^2]
  real(r8):: F_abs_sum                           ! total absorbed energy in column [W/m^2]
  real(r8):: F_sfc_pls                           ! upward radiative flux at snowpack top [W/m^2]
  real(r8):: F_btm_net                           ! net flux at bottom of snowpack [W/m^2]
  real(r8):: energy_sum                          ! sum of all energy terms; should be 0.0 [W/m^2]
  real(r8):: mu_not                              ! cosine of solar zenith angle (used locally) [frc]

  integer :: err_idx                             ! counter for number of times through error loop [nbr]
  real(r8):: pi                                  ! 3.1415...
  integer :: snw_shp_lcl(maxsnl+1:0)             ! Snow grain shape option:
                                                 ! 1=sphere; 2=spheroid; 3=hexagonal plate; 4=koch snowflake
  real(r8):: snw_fs_lcl(maxsnl+1:0)              ! Shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
                                                 ! 0=use recommended default value
                                                 ! others(0<fs<1)= use user-specified value
                                                 ! only activated when sno_shp > 1 (i.e. nonspherical)
  real(r8):: snw_ar_lcl(maxsnl+1:0)              ! % Aspect ratio: ratio of grain width to length
                                                 ! 0=use recommended default value
                                                 ! others(0.1<fs<20)= use user-specified value
                                                 ! only activated when sno_shp > 1 (i.e. nonspherical) 

     real(r8):: &
         diam_ice           , & ! effective snow grain diameter
         fs_sphd            , & ! shape factor for spheroid
         fs_hex0            , & ! shape factor for hexagonal plate
         fs_hex             , & ! shape factor for hexagonal plate (reference)
         fs_koch            , & ! shape factor for koch snowflake
         AR_tmp             , & ! aspect ratio for spheroid
         g_ice_Cg_tmp(7)    , & ! temporary for calculation of asymetry factor
         gg_ice_F07_tmp(7)  , & ! temporary for calculation of asymetry factor
         g_ice_F07          , & ! temporary for calculation of asymetry factor
         g_ice              , & ! asymmetry factor
         gg_F07_intp        , & ! temporary for calculation of asymetry factor (interpolated)
         g_Cg_intp          , & ! temporary for calculation of asymetry factor  (interpolated)
         R_1_omega_tmp      , & ! temporary for dust-snow mixing calculation
         C_dust_total           ! dust concentration

     integer :: atm_type_index  ! index for atmospheric type
     integer :: slr_zen         ! integer value of solar zenith angle

     ! SNICAR_AD new variables, follow sea-ice shortwave conventions
     real(r8):: &
        trndir(maxsnl+1:12)  , & ! solar beam down transmission from top
        trntdr(maxsnl+1:12)  , & ! total transmission to direct beam for layers above
        trndif(maxsnl+1:12)  , & ! diffuse transmission to diffuse beam for layers above
        rupdir(maxsnl+1:12)  , & ! reflectivity to direct radiation for layers below
        rupdif(maxsnl+1:12)  , & ! reflectivity to diffuse radiation for layers below
        rdndif(maxsnl+1:12)  , & ! reflectivity to diffuse radiation for layers above
        dfdir(maxsnl+1:12)   , & ! down-up flux at interface due to direct beam at top surface
        dfdif(maxsnl+1:12)   , & ! down-up flux at interface due to diffuse beam at top surface
        fdirdn(maxsnl+1:12)   , & 
        fdirup(maxsnl+1:12)   , &
        fdifdn(maxsnl+1:12)   , &
        fdifup(maxsnl+1:12)   , &
        F_up(maxsnl+1:12,numrad_snw)   , &
        F_dwn(maxsnl+1:12,numrad_snw)   , &
        F_net(maxsnl+1:12,numrad_snw)   , &
        dftmp(maxsnl+1:12)       ! temporary variable for down-up flux at interface

     real(r8):: &
        rdir(maxsnl+1:11)    , & ! layer reflectivity to direct radiation
        rdif_a(maxsnl+1:11)  , & ! layer reflectivity to diffuse radiation from above
        rdif_b(maxsnl+1:11)  , & ! layer reflectivity to diffuse radiation from below
        tdir(maxsnl+1:11)    , & ! layer transmission to direct radiation (solar beam + diffuse)
        tdif_a(maxsnl+1:11)  , & ! layer transmission to diffuse radiation from above
        tdif_b(maxsnl+1:11)  , & ! layer transmission to diffuse radiation from below
        trnlay(maxsnl+1:11)      ! solar beam transm for layer (direct beam only)


     real(r8):: &
         ts       , & ! layer delta-scaled extinction optical depth
         ws       , & ! layer delta-scaled single scattering albedo
         gs       , & ! layer delta-scaled asymmetry parameter
         extins   , & ! extinction
         alp      , & ! temporary for alpha
         gam      , & ! temporary for agamm
         amg      , & ! alp - gam
         apg      , & ! alp + gam
         ue       , & ! temporary for u
         refk     , & ! interface multiple scattering
         refkp1   , & ! interface multiple scattering for k+1
         refkm1   , & ! interface multiple scattering for k-1
         tdrrdir  , & ! direct tran times layer direct ref
         tdndif       ! total down diffuse = tot tran - direct tran

     real(r8) :: &
         alpha    , & ! term in direct reflectivity and transmissivity
         agamm    , & ! term in direct reflectivity and transmissivity
         el       , & ! term in alpha,agamm,n,u
         taus     , & ! scaled extinction optical depth
         omgs     , & ! scaled single particle scattering albedo
         asys     , & ! scaled asymmetry parameter
         u        , & ! term in diffuse reflectivity and transmissivity
         n        , & ! term in diffuse reflectivity and transmissivity
         lm       , & ! temporary for el
         mu       , & ! cosine solar zenith for either snow or water
         ne           ! temporary for n

     ! perpendicular and parallel relative to plane of incidence and scattering
     real(r8) :: &
         R1       , & ! perpendicular polarization reflection amplitude
         R2       , & ! parallel polarization reflection amplitude
         T1       , & ! perpendicular polarization transmission amplitude
         T2       , & ! parallel polarization transmission amplitude
         Rf_dir_a , & ! fresnel reflection to direct radiation
         Tf_dir_a , & ! fresnel transmission to direct radiation
         Rf_dif_a , & ! fresnel reflection to diff radiation from above
         Rf_dif_b , & ! fresnel reflection to diff radiation from below
         Tf_dif_a , & ! fresnel transmission to diff radiation from above
         Tf_dif_b     ! fresnel transmission to diff radiation from below

     real(r8) :: &
         gwt      , & ! gaussian weight
         swt      , & ! sum of weights
         trn      , & ! layer transmission
         rdr      , & ! rdir for gaussian integration
         tdr      , & ! tdir for gaussian integration
         smr      , & ! accumulator for rdif gaussian integration
         smt      , & ! accumulator for tdif gaussian integration
         exp_min      ! minimum exponential value

     integer :: &
         ng             , & ! gaussian integration index
         snl_btm_itf    , & ! index of bottom snow layer interfaces (1) [idx]
         ngmax = 8          ! gaussian integration index

     ! Gaussian integration angle and coefficients
     real(r8) :: difgauspt(1:8) , difgauswt(1:8)

     ! constants used in algorithm
     real(r8) :: &
         c0      = 0.0_r8     , &
         c1      = 1.0_r8     , &
         c2      = 2.0_r8     , &
         c3      = 3.0_r8     , &
         c4      = 4.0_r8     , &
         c6      = 6.0_r8     , &
         cp01    = 0.01_r8    , &
         cp5     = 0.5_r8     , &
         cp75    = 0.75_r8    , &
         c1p5    = 1.5_r8     , &
         trmin   = 0.001_r8   , &
         argmax  = 10.0_r8       ! maximum argument of exponential

     ! cconstant coefficients used for SZA parameterization
     real(r8) :: &
         sza_a0 =  0.085730_r8 , &
         sza_a1 = -0.630883_r8 , &
         sza_a2 =  1.303723_r8 , &
         sza_b0 =  1.467291_r8 , &
         sza_b1 = -3.338043_r8 , &
         sza_b2 =  6.807489_r8 , &
         puny   =  1.0e-11_r8  , &
         mu_75  =  0.2588_r8   , &  ! cosine of 75 degree
         SHR_CONST_PI     = 3.14159265358979323846_R8 , &
         snw_rds_min     = 54.526_r8         ,&
         min_snw = 1.0E-30_r8


     ! coefficients used for SZA parameterization
     real(r8) :: &
         sza_c1          , & ! coefficient, SZA parameteirzation
         sza_c0          , & ! coefficient, SZA parameterization
         sza_factor      , & ! factor used to adjust NIR direct albedo
         flx_sza_adjust  , & ! direct NIR flux adjustment from sza_factor
         mu0             , &   ! incident solar zenith angle
         ice_density_wgted

     !-----------------------NEW---------------------------------------------
     integer :: nbr_lyr      ! Total number of layers of ice and water
     integer :: nl_ice       ! number of layers of ice
     integer :: nl_wat       ! number of layers of water
     integer :: kfrsnl1      ! Fresnel layer between snow and ice
     integer :: kfrsnl2      ! Fresnel layer between ice and water
     integer :: kfrsnl3
     integer :: kfrsnl4
     integer :: ice_ri       ! ice refractive index dataset
   
     real(r8) :: rho_ice        ! Density of ice
   !   real(r8) :: rho_air        ! Density of air
     real(r8) :: ziwsum
     real(r8) :: refindx_re(numrad_snw)   !
     real(r8) :: refindx_im(numrad_snw)   ! ice refractive index
     real(r8) :: refindxwat_im(numrad_snw)   !
     real(r8) :: refindxwat_re(numrad_snw)   ! ice refractive index
     real(r8) :: ss_alb_wtr_avg_cld(numrad_snw)   ! ice refractive index
   !   real(r8) :: ss_alb_wtr_avg_clr(numrad_snw,90)   ! ice refractive index
   !   real(r8) :: ext_cff_mss_wtr_avg_clr(numrad_snw,90)   ! ice refractive index
     real(r8) :: ext_cff_mss_wtr_avg_cld(numrad_snw)   ! ice refractive index
     real(r8) :: temp1          ! transmissivity1
     real(r8) :: temp2          ! transmissivity2
     real(r8) :: nreal          ! Adjusted complex refractive index
     real(r8) :: nreal_wtr          ! Adjusted complex refractive index
     integer :: iwin           ! The number of layer at ice water interface
     real(r8) :: mu0n           ! 
     real(r8) :: rintfc


     !-----------------------------------------------------------------------
#ifdef MODAL_AER
         !mgf++
         integer :: idx_bcint_icerds  ! index of ice effective radius for optical properties lookup table
         integer :: idx_bcint_nclrds  ! index of within-ice BC effective radius for optical properties lookup table
         integer :: idx_bcext_nclrds  ! index of external BC effective radius for optical properties lookup table
         real(r8):: enh_fct           ! extinction/absorption enhancement factor for within-ice BC
         real(r8):: tmp1              ! temporary variable
         !mgf--
#endif

      ! Constants for non-spherical ice particles and dust-snow internal mixing
      real(r8) :: g_b2(7)
      real(r8) :: g_b1(7)
      real(r8) :: g_b0(7)
      real(r8) :: g_F07_c2(7)
      real(r8) :: g_F07_c1(7)
      real(r8) :: g_F07_c0(7)
      real(r8) :: g_F07_p2(7)
      real(r8) :: g_F07_p1(7)
      real(r8) :: g_F07_p0(7)
      real(r8) :: dust_clear_d0(3)
      real(r8) :: dust_clear_d1(3)
      real(r8) :: dust_clear_d2(3)
      real(r8) :: dust_cloudy_d0(3)
      real(r8) :: dust_cloudy_d1(3)
      real(r8) :: dust_cloudy_d2(3)

      real(r8), allocatable :: lyr_typ(:)  !(maxsnl+1:nbr_lyr) layer type
      real(r8) :: sca_cff_vlm_airbbl_lcl(nl_lake)
      real(r8) :: abs_cff_mss_ice_lcl(nl_lake)
      real(r8) :: vlm_frac_air(nl_lake)


      !!! factors for considering snow grain shape
      data g_b0(:) / 9.76029E-01_r8, 9.67798E-01_r8, 1.00111E+00_r8, 1.00224E+00_r8,&
                     9.64295E-01_r8, 9.97475E-01_r8, 9.97475E-01_r8/
      data g_b1(:) / 5.21042E-01_r8, 4.96181E-01_r8, 1.83711E-01_r8, 1.37082E-01_r8,&
                     5.50598E-02_r8, 8.48743E-02_r8, 8.48743E-02_r8/
      data g_b2(:) /-2.66792E-04_r8, 1.14088E-03_r8, 2.37011E-04_r8,-2.35905E-04_r8,&
                     8.40449E-04_r8,-4.71484E-04_r8,-4.71484E-04_r8/

      data g_F07_c2(:) / 1.349959E-1_r8, 1.115697E-1_r8, 9.853958E-2_r8, 5.557793E-2_r8,&
                        -1.233493E-1_r8, 0.0_r8, 0.0_r8/
      data g_F07_c1(:) /-3.987320E-1_r8,-3.723287E-1_r8,-3.924784E-1_r8,-3.259404E-1_r8,&
                         4.429054E-2_r8,-1.726586E-1_r8,-1.726586E-1_r8/
      data g_F07_c0(:) / 7.938904E-1_r8, 8.030084E-1_r8, 8.513932E-1_r8, 8.692241E-1_r8,&
                         7.085850E-1_r8, 6.412701E-1_r8, 6.412701E-1_r8/
      data g_F07_p2(:) / 3.165543E-3_r8, 2.014810E-3_r8, 1.780838E-3_r8, 6.987734E-4_r8,&
                        -1.882932E-2_r8,-2.277872E-2_r8,-2.277872E-2_r8/
      data g_F07_p1(:) / 1.140557E-1_r8, 1.143152E-1_r8, 1.143814E-1_r8, 1.071238E-1_r8,&
                         1.353873E-1_r8, 1.914431E-1_r8, 1.914431E-1_r8/
      data g_F07_p0(:) / 5.292852E-1_r8, 5.425909E-1_r8, 5.601598E-1_r8, 6.023407E-1_r8,&
                         6.473899E-1_r8, 4.634944E-1_r8, 4.634944E-1_r8/

      !!! factors for considring dust-snow internal mixing
      data dust_clear_d0(:) /1.0413E+00_r8,1.0168E+00_r8,1.0189E+00_r8/
      data dust_clear_d1(:) /1.0016E+00_r8,1.0070E+00_r8,1.0840E+00_r8/
      data dust_clear_d2(:) /2.4208E-01_r8,1.5300E-03_r8,1.1230E-04_r8/

      data dust_cloudy_d0(:) /1.0388E+00_r8,1.0167E+00_r8,1.0189E+00_r8/
      data dust_cloudy_d1(:) /1.0015E+00_r8,1.0061E+00_r8,1.0823E+00_r8/
      data dust_cloudy_d2(:) /2.5973E-01_r8,1.6200E-03_r8,1.1721E-04_r8/
  !--------------------

  real(r8) hs        ! net ground heat flux into the surface
  real(r8) dhsdT     ! temperature derivative of "hs"
  real(r8) heatavail ! available energy for melting or freezing (J/m^2)
  real(r8) heatrem   ! energy residual or loss after melting or freezing
  real(r8) melt      ! actual melting (+) or freezing (-) [kg/m2]
  real(r8) xmf       ! total per-column latent heat abs. from phase change  (J/m^2)
  !--------------------

  real(r8) ocvts     ! (cwat*(t_lake[n  ])*dz_lake
  real(r8) ncvts     ! (cwat*(t_lake[n+1])*dz_lake
  real(r8) esum1     ! temp for checking energy (J/m^2)
  real(r8) esum2     ! ""
  real(r8) zsum      ! temp for putting ice at the top during convection (m)
  real(r8) errsoi    ! soil/lake energy conservation error (W/m^2)

  real(r8) iceav     ! used in calc aver ice for convectively mixed layers
  real(r8) qav       ! used in calc aver heat content for conv. mixed layers
  real(r8) tav       ! used in aver temp for convectively mixed layers
  real(r8) tav_froz  ! used in aver temp for convectively mixed layers (C)
  real(r8) tav_unfr  ! "
  real(r8) nav       ! used in aver temp for convectively mixed layers

  real(r8) fevpg_lim ! temporary evap_soi limited by top snow layer content [mm/s]
  real(r8) scv_temp  ! temporary h2osno [kg/m^2]
  real(r8) tmp       !
  real(r8) h_fin     !
  real(r8) h_finDT   !
  real(r8) del_T_grnd   !
!  real(r8) savedtke1

  integer iter       ! iteration index
  integer convernum  ! number of time when del_T_grnd < 0.01
  integer nmozsgn    ! number of times moz changes sign

! assign iteration parameters
  integer, parameter :: itmax  = 40   ! maximum number of iteration
  integer, parameter :: itmin  = 6    ! minimum number of iteration
  real(r8),parameter :: delmax = 3.0  ! maximum change in lake temperature [K]
  real(r8),parameter :: dtmin  = 0.01 ! max limit for temperature convergence [K]
  real(r8),parameter :: dlemin = 0.1  ! max limit for energy flux convergence [w/m2]

  !--------------------

  integer nl_sls  ! abs(snl)+nl_lake+nl_soil
  integer snl     ! number of snow layers (minimum -5)
  integer lb      ! lower bound of arrays
  integer jprime  ! j - nl_lake

  integer i,j     ! do loop or array index
  integer aer

  real(r8), intent(inout) :: &
  mss_bcpho    ( maxsnl+1:0 ), &! mass of hydrophobic BC in snow  (col,lyr) [kg]
  mss_bcphi    ( maxsnl+1:0 ), &! mass of hydrophillic BC in snow (col,lyr) [kg]
  mss_ocpho    ( maxsnl+1:0 ), &! mass of hydrophobic OC in snow  (col,lyr) [kg]
  mss_ocphi    ( maxsnl+1:0 ), &! mass of hydrophillic OC in snow (col,lyr) [kg]
  mss_dst1     ( maxsnl+1:0 ), &! mass of dust species 1 in snow  (col,lyr) [kg]
  mss_dst2     ( maxsnl+1:0 ), &! mass of dust species 2 in snow  (col,lyr) [kg]
  mss_dst3     ( maxsnl+1:0 ), &! mass of dust species 3 in snow  (col,lyr) [kg]
  mss_dst4     ( maxsnl+1:0 )   ! mass of dust species 4 in snow  (col,lyr) [kg]

  logical do_capsnow      ! true => DO snow capping

  real(r8) snwcp_ice                        !excess precipitation due to snow capping [kg m-2 s-1]
  real(r8) mss_cnc_bcphi ( maxsnl+1:0 )     !mass concentration of hydrophilic BC (col,lyr) [kg/kg]
  real(r8) mss_cnc_bcpho ( maxsnl+1:0 )     !mass concentration of hydrophobic BC (col,lyr) [kg/kg]
  real(r8) mss_cnc_ocphi ( maxsnl+1:0 )     !mass concentration of hydrophilic OC (col,lyr) [kg/kg]
  real(r8) mss_cnc_ocpho ( maxsnl+1:0 )     !mass concentration of hydrophobic OC (col,lyr) [kg/kg]
  real(r8) mss_cnc_dst1  ( maxsnl+1:0 )     !mass concentration of dust aerosol species 1 (col,lyr) [kg/kg]
  real(r8) mss_cnc_dst2  ( maxsnl+1:0 )     !mass concentration of dust aerosol species 2 (col,lyr) [kg/kg]
  real(r8) mss_cnc_dst3  ( maxsnl+1:0 )     !mass concentration of dust aerosol species 3 (col,lyr) [kg/kg]
  real(r8) mss_cnc_dst4  ( maxsnl+1:0 )     !mass concentration of dust aerosol species 4 (col,lyr) [kg/kg]
  real(r8) mss_cnc_aer_in (maxsnl+1:0,sno_nbr_aer) ! mass concentration of aerosol species for forcing calculation (zero) (col,lyr,aer) [kg kg-1]
  real(r8) snw_rds_in        ( maxsnl+1: 0) 
! ======================================================================
!*[1] constants and model parameters
! ======================================================================

! constants for lake temperature model
  write(*,*)'twostream_first',DEF_USE_TWOSTREAM
      za = (/0.5, 0.6/)
      cwat = cpliq*denh2o     ! water heat capacity per unit volume
      cice_eff = cpice*denh2o ! use water density because layer depth is not adjusted for freezing
      cfus = hfus*denh2o      ! latent heat per unit volume
      tkice_eff = tkice * denice/denh2o ! effective conductivity since layer depth is constant
      emg = 0.97              ! surface emissivity

! define snow layer on ice lake
      snl = 0
      do j=maxsnl+1,0
         if(wliq_soisno(j)+wice_soisno(j)>0.) snl=snl-1
      enddo
      lb = snl + 1

! latent heat
      if (t_grnd > tfrz )then
         htvp = hvap
      else
         htvp = hsub
      end if

! define levels
      z_lake(1) = dz_lake(1) / 2.
      do j = 2, nl_lake
         z_lake(j) = z_lake(j-1) + (dz_lake(j-1) + dz_lake(j))/2.
      end do

! Base on lake depth, assuming that small lakes are likely to be shallower
! Estimate crudely based on lake depth
      if (z_lake(nl_lake) < 4.) then
          idlak = 1
          fetch = 100. ! shallow lake
      else
          idlak = 2
          fetch = 25.*z_lake(nl_lake) ! deep lake
      end if


! ======================================================================
!*[2] pre-processing for the calcilation of the surface temperature and fluxes
! ======================================================================

      ! IF (.not. DEF_USE_SNICAR .or. present(urban_call)) THEN
      !    IF (DEF_USE_TWOSTREAM) THEN
      !       if (snl == 0) then
      !          ! calculate the nir fraction of absorbed solar.
      !          betaprime = (sabgn+sabgnd)/max(1.e-5,sabgv+sabgn+sabgvd+sabgnd)
      !          betavis = 0. ! The fraction of the visible (e.g. vis not nir from atm) sunlight
      !                      ! absorbed in ~1 m of water (the surface layer za_lake).
      !                      ! This is roughly the fraction over 700 nm but may depend on the details
      !                      ! of atmospheric radiative transfer.
      !                      ! As long as NIR = 700 nm and up, this can be zero.
      !          betaprime = betaprime + (1.0-betaprime)*betavis
      !       else
      !          ! or frozen but no snow layers or
      !          ! currently ignor the transmission of solar in snow and ice layers
      !          ! to be updated in the future version
      !          betaprime = 1.0
      !       end if
      !    ELSE 
      !       if (snl == 0) then
      !          ! calculate the nir fraction of absorbed solar.
      !          betaprime = (forc_soll+forc_solld)/max(1.e-5,forc_sols+forc_soll+forc_solsd+forc_solld)
      !          betavis = 0. ! The fraction of the visible (e.g. vis not nir from atm) sunlight
      !                      ! absorbed in ~1 m of water (the surface layer za_lake).
      !                      ! This is roughly the fraction over 700 nm but may depend on the details
      !                      ! of atmospheric radiative transfer.
      !                      ! As long as NIR = 700 nm and up, this can be zero.
      !          betaprime = betaprime + (1.0-betaprime)*betavis
      !       else
      !          ! or frozen but no snow layers or
      !          ! currently ignor the transmission of solar in snow and ice layers
      !          ! to be updated in the future version
      !          betaprime = 1.0
      !       end if 
      !    END IF 

      ! ELSE
      !    ! calculate the nir fraction of absorbed solar.
      !    betaprime = (forc_soll+forc_solld)/max(1.e-5,forc_sols+forc_soll+forc_solsd+forc_solld)
      !    betavis = 0. ! The fraction of the visible (e.g. vis not nir from atm) sunlight
      !                 ! absorbed in ~1 m of water (the surface layer za_lake).
      !                 ! This is roughly the fraction over 700 nm but may depend on the details
      !                 ! of atmospheric radiative transfer.
      !                 ! As long as NIR = 700 nm and up, this can be zero.
      !    betaprime = betaprime + (1.0-betaprime)*betavis
      ! ENDIF

      IF (.not. DEF_USE_SNICAR .or. present(urban_call)) THEN
         if (snl == 0) then
            ! calculate the nir fraction of absorbed solar.
            betaprime = (forc_soll+forc_solld)/max(1.e-5,forc_sols+forc_soll+forc_solsd+forc_solld)
            betavis = 0. ! The fraction of the visible (e.g. vis not nir from atm) sunlight
                         ! absorbed in ~1 m of water (the surface layer za_lake).
                         ! This is roughly the fraction over 700 nm but may depend on the details
                         ! of atmospheric radiative transfer.
                         ! As long as NIR = 700 nm and up, this can be zero.
            betaprime = betaprime + (1.0-betaprime)*betavis
         else
            ! or frozen but no snow layers or
            ! currently ignor the transmission of solar in snow and ice layers
            ! to be updated in the future version
            betaprime = 1.0
         end if

      ELSE
         ! calculate the nir fraction of absorbed solar.
         betaprime = (forc_soll+forc_solld)/max(1.e-5,forc_sols+forc_soll+forc_solsd+forc_solld)
         betavis = 0. ! The fraction of the visible (e.g. vis not nir from atm) sunlight
                      ! absorbed in ~1 m of water (the surface layer za_lake).
                      ! This is roughly the fraction over 700 nm but may depend on the details
                      ! of atmospheric radiative transfer.
                      ! As long as NIR = 700 nm and up, this can be zero.
         betaprime = betaprime + (1.0-betaprime)*betavis
      ENDIF

      call qsadv(t_grnd,forc_psrf,eg,degdT,qsatg,qsatgdT)
! potential temperatur at the reference height
      beta1=1.       ! -  (in computing W_*)
      zii = 1000.    ! m  (pbl height)
      thm = forc_t + 0.0098*forc_hgt_t  ! intermediate variable equivalent to
                                        ! forc_t*(pgcm/forc_psrf)**(rgas/cpair)
      th = forc_t*(100000./forc_psrf)**(rgas/cpair) ! potential T
      thv = th*(1.+0.61*forc_q)         ! virtual potential T
      ur = max(0.1,sqrt(forc_us*forc_us+forc_vs*forc_vs))   ! limit set to 0.1

! Initialization variables
      nmozsgn = 0
      obuold = 0.
      dth   = thm-t_grnd
      dqh   = forc_q-qsatg
      dthv  = dth*(1.+0.61*forc_q)+0.61*th*dqh
      zldis = forc_hgt_u-0.

! Roughness lengths, allow all roughness lengths to be prognostic
      ustar=0.06
      wc=0.5

    ! Kinematic viscosity of dry air (m2/s)- Andreas (1989) CRREL Rep. 89-11
      visa=1.326e-5*(1.+6.542e-3*(forc_t-tfrz) &
           + 8.301e-6*(forc_t-tfrz)**2 - 4.84e-9*(forc_t-tfrz)**3)

      cur = cur0 + curm * exp( max( -(fetch*grav/ur/ur)**(1./3.)/fcrit, & ! Fetch-limited
                                    -(z_lake(nl_lake)*grav)**0.5/ur ) )   ! depth-limited

      if(dthv.ge.0.) then
         um=max(ur,0.1)
      else
         um=sqrt(ur*ur+wc*wc)
      endif

      do i=1,5
         z0mg=0.013*ustar*ustar/grav+0.11*visa/ustar
         ustar=vonkar*um/log(zldis/z0mg)
      enddo

      call roughness_lake (snl,t_grnd,t_lake(1),lake_icefrac(1),forc_psrf,&
                           cur,ustar,z0mg,z0hg,z0qg)

      call moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu)

      if (snl == 0) then
         dzsur = dz_lake(1)/2.
      else
         dzsur = z_soisno(lb)-zi_soisno(lb-1)
      end if


      iter = 1
      del_T_grnd = 1.0    ! t_grnd diff
      convernum = 0       ! number of time when del_T_grnd <= 0.01


! ======================================================================
!*[3] Begin stability iteration and temperature and fluxes calculation
! ======================================================================

      IF( .not. DEF_USE_TWOSTREAM) then
      ! =====================================
         ITERATION : DO WHILE (iter <= itmax)
      ! =====================================

            t_grnd_bef = t_grnd

            if (t_grnd_bef > tfrz .and. t_lake(1) > tfrz .and. snl == 0) then
               tksur = savedtke1       !water molecular conductivity
               tsur = t_lake(1)
               htvp = hvap
            else if (snl == 0) then !frozen but no snow layers
               tksur = tkice        ! This is an approximation because the whole layer may not be frozen, and it is not
                                    ! accounting for the physical (but not nominal) expansion of the frozen layer.
               tsur = t_lake(1)
               htvp = hsub
            else
            ! need to calculate thermal conductivity of the top snow layer
               rhosnow = (wice_soisno(lb)+wliq_soisno(lb))/dz_soisno(lb)
               tksur = tkair + (7.75e-5*rhosnow + 1.105e-6*rhosnow*rhosnow)*(tkice-tkair)
               tsur = t_soisno(lb)
               htvp = hsub
            end if

   ! Evaluated stability-dependent variables using moz from prior iteration
            displax = 0.
            if (DEF_USE_CBL_HEIGHT) then
               call moninobuk_leddy(forc_hgt_u,forc_hgt_t,forc_hgt_q,displax,z0mg,z0hg,z0qg,obu,um, hpbl, &
                           ustar,fh2m,fq2m,fm10m,fm,fh,fq)
            else
               call moninobuk(forc_hgt_u,forc_hgt_t,forc_hgt_q,displax,z0mg,z0hg,z0qg,obu,um,&
                           ustar,fh2m,fq2m,fm10m,fm,fh,fq)
            endif

   ! Get derivative of fluxes with repect to ground temperature
            ram    = 1./(ustar*ustar/um)
            rah    = 1./(vonkar/fh*ustar)
            raw    = 1./(vonkar/fq*ustar)
            stftg3 = emg*stefnc*t_grnd_bef*t_grnd_bef*t_grnd_bef

            ax  = betaprime*sabg + emg*forc_frl + 3.*stftg3*t_grnd_bef &
               + forc_rhoair*cpair/rah*thm &
               - htvp*forc_rhoair/raw*(qsatg-qsatgdT*t_grnd_bef - forc_q) &
               + tksur*tsur/dzsur

            bx  = 4.*stftg3 + forc_rhoair*cpair/rah &
               + htvp*forc_rhoair/raw*qsatgdT + tksur/dzsur

            t_grnd = ax/bx

         !-----------------------------------------------------------------
         ! h_fin = betaprime*sabg + emg*forc_frl + 3.*stftg3*t_grnd_bef & !
         !     + forc_rhoair*cpair/rah*thm &                              !
         !     - htvp*forc_rhoair/raw*(qsatg-qsatgdT*t_grnd_bef - forc_q) !
         ! h_finDT = 4.*stftg3 + forc_rhoair*cpair/rah &                  !
         !     + htvp*forc_rhoair/raw*qsatgdT                             !
         ! del_T_grnd = t_grnd - t_grnd_bef                               !
         !----------------------------------------------------------------!

   ! surface fluxes of momentum, sensible and latent
   ! using ground temperatures from previous time step

            fseng = forc_rhoair*cpair*(t_grnd-thm)/rah
            fevpg = forc_rhoair*(qsatg+qsatgdT*(t_grnd-t_grnd_bef)-forc_q)/raw

            call qsadv(t_grnd,forc_psrf,eg,degdT,qsatg,qsatgdT)
            dth = thm-t_grnd
            dqh = forc_q-qsatg
            tstar = vonkar/fh*dth
            qstar = vonkar/fq*dqh
            thvstar = tstar*(1.+0.61*forc_q)+0.61*th*qstar
            zeta = zldis*vonkar*grav*thvstar/(ustar**2*thv)
            if(zeta >= 0.) then     !stable
            zeta = min(2.,max(zeta,1.e-6))
            else                    !unstable
            zeta = max(-100.,min(zeta,-1.e-6))
            endif
            obu = zldis/zeta
            if(zeta >= 0.)then
            um = max(ur,0.1)
            else
            if (DEF_USE_CBL_HEIGHT) then !//TODO: Shaofeng, 2023.05.18
               zii = max(5.*forc_hgt_u,hpbl)
            endif !//TODO: Shaofeng, 2023.05.18
            wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
            wc2 = beta1*beta1*(wc*wc)
            um = sqrt(ur*ur+wc2)
            endif

            call roughness_lake (snl,t_grnd,t_lake(1),lake_icefrac(1),forc_psrf,&
                                 cur,ustar,z0mg,z0hg,z0qg)

            iter = iter + 1
            del_T_grnd = abs(t_grnd - t_grnd_bef)

            if(iter .gt. itmin) then
               if(del_T_grnd <= dtmin) then
                  convernum = convernum + 1
               end if
               if(convernum >= 4) EXIT
            endif

      ! ===============================================
         END DO ITERATION   ! end of stability iteration
      ! ===============================================

   !*----------------------------------------------------------------------
   !*Zack Subin, 3/27/09
   !*Since they are now a function of whatever t_grnd was before cooling
   !*to freezing temperature, then this value should be used in the derivative correction term.
   !*Allow convection if ground temp is colder than lake but warmer than 4C, or warmer than
   !*lake which is warmer than freezing but less than 4C.
         tdmax = tfrz + 4.0
         if ( (snl < 0 .or. t_lake(1) <= tfrz) .and. t_grnd > tfrz) then
            t_grnd_bef = t_grnd
            t_grnd = tfrz
            fseng = forc_rhoair*cpair*(t_grnd-thm)/rah
            fevpg = forc_rhoair*(qsatg+qsatgdT*(t_grnd-t_grnd_bef)-forc_q)/raw
         else if ( (t_lake(1) > t_grnd .and. t_grnd > tdmax) .or. &
                  (t_lake(1) < t_grnd .and. t_lake(1) > tfrz .and. t_grnd < tdmax) ) then
                  ! Convective mixing will occur at surface
            t_grnd_bef = t_grnd
            t_grnd = t_lake(1)
            fseng = forc_rhoair*cpair*(t_grnd-thm)/rah
            fevpg = forc_rhoair*(qsatg+qsatgdT*(t_grnd-t_grnd_bef)-forc_q)/raw
         end if
      
   !*----------------------------------------------------------------------

   ! net longwave from ground to atmosphere
         stftg3 = emg*stefnc*t_grnd_bef*t_grnd_bef*t_grnd_bef
         olrg = (1.-emg)*forc_frl + emg*stefnc*t_grnd_bef**4 + 4.*stftg3*(t_grnd - t_grnd_bef)
         if (t_grnd > tfrz )then
            htvp = hvap
         else
            htvp = hsub
         end if

   !The actual heat flux from the ground interface into the lake, not including the light that penetrates the surface.
         fgrnd1 = betaprime*sabg + forc_frl - olrg - fseng - htvp*fevpg

         ! January 12, 2023 by Yongjiu Dai
         IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
            hs = sabg_lyr(lb) + forc_frl - olrg - fseng - htvp*fevpg
            dhsdT = 0.0
         ENDIF

   !------------------------------------------------------------
   ! Set up vector r and vectors a, b, c that define tridiagonal matrix
   ! snow and lake and soil layer temperature
   !------------------------------------------------------------

   !------------------------------------------------------------
   ! Lake density
   !------------------------------------------------------------
      ! IF(DEF_USE_TWOSTREAM) then

         do j = 1, nl_lake
            rhow(j) = (1.-lake_icefrac(j))*denh2o*(1.0-1.9549e-05*(abs(t_lake(j)-277.))**1.68) &
                        + lake_icefrac(j)*denice
            ! allow for ice fraction; assume constant ice density.
            ! this is not the correct average-weighting but that's OK because the density will only
            ! be used for convection for lakes with ice, and the ice fraction will dominate the
            ! density differences between layers.
            ! using this average will make sure that surface ice is treated properly during
            ! convective mixing.
         end do
      ! ENDIF
   !------------------------------------------------------------
   ! Diffusivity and implied thermal "conductivity" = diffusivity * cwat
   !------------------------------------------------------------

         do j = 1, nl_lake
            cv_lake(j) = dz_lake(j) * (cwat*(1.-lake_icefrac(j)) + cice_eff*lake_icefrac(j))
         end do

         call hConductivity_lake(nl_lake,snl,t_grnd,&
                                 z_lake,t_lake,lake_icefrac,rhow,&
                                 dlat,ustar,z0mg,lakedepth,depthcrit,tk_lake,savedtke1)

   !------------------------------------------------------------
   ! Set the thermal properties of the snow above frozen lake and underlying soil
   ! and check initial energy content.
   !------------------------------------------------------------

         lb = snl+1
         do i = 1, nl_soil
            vf_water(i) = wliq_soisno(i)/(dz_soisno(i)*denh2o)
            vf_ice(i) = wice_soisno(i)/(dz_soisno(i)*denice)
            CALL soil_hcap_cond(vf_gravels(i),vf_om(i),vf_sand(i),porsl(i),&
                              wf_gravels(i),wf_sand(i),k_solids(i),&
                              csol(i),dkdry(i),dksatu(i),dksatf(i),&
                              BA_alpha(i),BA_beta(i),&
                              t_soisno(i),vf_water(i),vf_ice(i),hcap(i),thk(i))
            cv_soisno(i) = hcap(i)*dz_soisno(i)
         enddo

   ! Snow heat capacity and conductivity
         if(lb <=0 )then
         do j = lb, 0
            cv_soisno(j) = cpliq*wliq_soisno(j) + cpice*wice_soisno(j)
            rhosnow = (wice_soisno(j)+wliq_soisno(j))/dz_soisno(j)
            thk(j) = tkair + (7.75e-5*rhosnow + 1.105e-6*rhosnow*rhosnow)*(tkice-tkair)
         enddo
         endif

   ! Thermal conductivity at the layer interface
         do i = lb, nl_soil-1

   ! the following consideration is try to avoid the snow conductivity
   ! to be dominant in the thermal conductivity of the interface.
   ! Because when the distance of bottom snow node to the interfacee
   ! is larger than that of interface to top soil node,
   ! the snow thermal conductivity will be dominant, and the result is that
   ! lees heat tranfer between snow and soil

   ! modified by Nan Wei, 08/25/2014
            if (i /= 0) then
               tk_soisno(i) = thk(i)*thk(i+1)*(z_soisno(i+1)-z_soisno(i)) &
                     /(thk(i)*(z_soisno(i+1)-zi_soisno(i))+thk(i+1)*(zi_soisno(i)-z_soisno(i)))
            else
               tk_soisno(i) = thk(i)
            end if
         end do
         tk_soisno(nl_soil) = 0.
         tktopsoil = thk(1)

   ! Sum cv_lake*t_lake for energy check
   ! Include latent heat term, and use tfrz as reference temperature
   ! to prevent abrupt change in heat content due to changing heat capacity with phase change.

         ! This will need to be over all soil / lake / snow layers. Lake is below.
         ocvts = 0.
         do j = 1, nl_lake
            ocvts = ocvts + cv_lake(j)*(t_lake(j)-tfrz) + cfus*dz_lake(j)*(1.-lake_icefrac(j))
         end do

         ! Now do for soil / snow layers
         do j = lb, nl_soil
            ocvts = ocvts + cv_soisno(j)*(t_soisno(j)-tfrz) + hfus*wliq_soisno(j)
            if (j == 1 .and. scv > 0. .and. j == lb) then
               ocvts = ocvts - scv*hfus
            end if
         end do
      ! ENDIF   ! End copy

      !------------------------------------------------------------
      ! Lake density
      !------------------------------------------------------------
      ELSE 
         do j = 1, nl_lake
            rhow(j) = (1.-lake_icefrac(j))*denh2o*(1.0-1.9549e-05*(abs(t_lake(j)-277.))**1.68) &
                        + lake_icefrac(j)*denice
            ! allow for ice fraction; assume constant ice density.
            ! this is not the correct average-weighting but that's OK because the density will only
            ! be used for convection for lakes with ice, and the ice fraction will dominate the
            ! density differences between layers.
            ! using this average will make sure that surface ice is treated properly during
            ! convective mixing.
         end do
      ENDIF
      ! Set up solar source terms (phix)
      ! Modified January 12, 2023 by Yongjiu Dai
      IF (.not. DEF_USE_SNICAR .or. present(urban_call)) THEN
         IF(DEF_USE_ORIGINAL)THEN
            write(*,*)'t_grnd',t_grnd ,' t_lake(1)', t_lake(1) , tfrz, lake_icefrac(1)
            if ((t_grnd > tfrz .and. t_lake(1) > tfrz .and. snl == 0)) then      !no snow cover, unfrozen layer lakes
               do j = 1, nl_lake
                  ! extinction coefficient from surface data (1/m), if no eta from surface data,
                  ! set eta, the extinction coefficient, according to L Hakanson, Aquatic Sciences, 1995
                  ! (regression of secchi depth with lake depth for small glacial basin lakes), and the
                  ! Poole & Atkins expression for extinction coeffient of 1.7 / secchi Depth (m).

                  eta = 1.1925*max(lakedepth,1.)**(-0.424)
                  zin  = z_lake(j) - 0.5*dz_lake(j)
                  zout = z_lake(j) + 0.5*dz_lake(j)
                  rsfin  = exp( -eta*max(  zin-za(idlak),0. ) )  ! the radiation within surface layer (z<za)
                  rsfout = exp( -eta*max( zout-za(idlak),0. ) )  ! is considered fixed at (1-beta)*sabg
                                                               ! i.e, max(z-za, 0)
                  ! Let rsfout for bottom layer go into soil.
                  ! This looks like it should be robust even for pathological cases,
                  ! like lakes thinner than za(idlak).

                  phi(j) = (rsfin-rsfout) * sabg * (1.-betaprime)
                  if (j == nl_lake) phi_soil = rsfout * sabg * (1.-betaprime)
               end do
            else if (snl == 0) then     !no snow-covered layers, but partially frozen
               phi(1) = sabg * (1.-betaprime)
               phi(2:nl_lake) = 0.
               phi_soil = 0.
            else   ! snow covered, this should be improved upon; Mironov 2002 suggests that SW can penetrate thin ice and may
                  ! cause spring convection.
               phi(:) = 0.
               phi_soil = 0.
            end if
            do i = 1, nl_lake
               write(*,*)i,phi(i)
            enddo 
            write(*,*)'phi_soil',phi_soil
         ! ENDIF
         ELSE IF(DEF_USE_ONEBANDICE)THEN
         ! IF(DEF_USE_ONEBANDICE)THEN
            if ((t_grnd > tfrz .and. t_lake(1) > tfrz .and. snl == 0)) then      !no snow cover, unfrozen layer lakes
               do j = 1, nl_lake
                  ! extinction coefficient from surface data (1/m), if no eta from surface data,
                  ! set eta, the extinction coefficient, according to L Hakanson, Aquatic Sciences, 1995
                  ! (regression of secchi depth with lake depth for small glacial basin lakes), and the
                  ! Poole & Atkins expression for extinction coeffient of 1.7 / secchi Depth (m).

                  etav = 0.239012444
                  zin  = z_lake(j) - 0.5*dz_lake(j)
                  zout = z_lake(j) + 0.5*dz_lake(j)
                  rsfin  = exp( -etav*max(  zin-za(idlak),0. ) )  ! the radiation within surface layer (z<za)
                  rsfout = exp( -etav*max( zout-za(idlak),0. ) )  ! is considered fixed at (1-beta)*sabg
                                                               ! i.e, max(z-za, 0)
                  ! Let rsfout for bottom layer go into soil.
                  ! This looks like it should be robust even for pathological cases,
                  ! like lakes thinner than za(idlak).

                  phi(j) = (rsfin-rsfout) * sabg * (1.-betaprime)
                  if (j == nl_lake) phi_soil = rsfout * sabg * (1.-betaprime)
               end do
            else if (snl == 0) then     !no snow-covered layers, but partially frozen
               do j = 1, nl_lake
                  etav = 0.239012444
                  etaiv = 1.5
                  zin  = z_lake(j) - 0.5*dz_lake(j)
                  zout = z_lake(j) + 0.5*dz_lake(j)
                  rsfin  = exp( -etav*max(  zin-za(idlak),0. ) )  ! the radiation within surface layer (z<za)
                  rsfout = exp( -etav*max( zout-za(idlak),0. ) )
                  rsfinv  = exp( -etaiv*max(  zin-za(idlak),0. ) ) 
                  rsfoutv = exp( -etaiv*max( zout-za(idlak),0. ) )
                  if(j == 1) then
                     phi(j) = sabg * (1.-betaprime) * (1. - (1-lake_icefrac(j)) * rsfout - &
                              lake_icefrac(j)  * rsfoutv  )
                  else
                     phi(j) = sabg * (1.-betaprime) * ((1-lake_icefrac(j-1)) * rsfin + lake_icefrac(j-1) * rsfinv  &
                                                      -  (1-lake_icefrac(j)) * rsfout  - lake_icefrac(j) * rsfoutv )
                  end if
                     if (j == nl_lake) phi_soil = (1-lake_icefrac(j)) * sabg * (1.-betaprime) * rsfout + &
                                                lake_icefrac(j)  * sabg * (1.-betaprime) * rsfoutv 
               end do
            else   ! snow covered, this should be improved upon; Mironov 2002 suggests that SW can penetrate thin ice and may
                  ! cause spring convection.
               phi(:) = 0.
               phi_soil = 0.
            end if
         ! ENDIF
         ELSE IF(DEF_USE_TWOBANDICE)THEN
            if ((t_grnd > tfrz .and. t_lake(1) > tfrz .and. snl == 0)) then      !no snow cover, unfrozen layer lakes
               do j = 1, nl_lake
                  ! extinction coefficient from surface data (1/m), if no eta from surface data,
                  ! set eta, the extinction coefficient, according to L Hakanson, Aquatic Sciences, 1995
                  ! (regression of secchi depth with lake depth for small glacial basin lakes), and the
                  ! Poole & Atkins expression for extinction coeffient of 1.7 / secchi Depth (m).
                  etav = 0.239012444
                  etan = 7
                  zin  = z_lake(j) - 0.5*dz_lake(j)
                  zout = z_lake(j) + 0.5*dz_lake(j)
                  rsfinvw  = exp( -etav*max(  zin-za(idlak),0. ) )  ! the radiation within surface layer (z<za)
                  rsfoutvw = exp( -etav*max( zout-za(idlak),0. ) )  ! is considered fixed at (1-beta)*sabg
                  rsfinnw  = exp( -etan*max(  zin-za(idlak),0. ) )
                  rsfoutnw  = exp( -etan*max( zout-za(idlak),0. ) )                                           ! i.e, max(z-za, 0)
                  ! Let rsfout for bottom layer go into soil.
                  ! This looks like it should be robust even for pathological cases,
                  ! like lakes thinner than za(idlak).

                  phi(j) = (1-betaprime)*((rsfinvw-rsfoutvw) * sabg * (1.-betaprime)+(rsfinnw-rsfoutnw) * sabg * betaprime)
                  if (j == nl_lake) phi_soil = (1-betaprime)*(rsfoutvw * sabg * (1.-betaprime)+ rsfoutnw * sabg * betaprime)
               end do
            else if (snl == 0) then     !no snow-covered layers, but partially frozen
               do j = 1, nl_lake
                  etav = 0.239012444
                  etan = 7
                  etaiv = 1.5
                  etain = 20
                  zin  = z_lake(j) - 0.5*dz_lake(j)
                  zout = z_lake(j) + 0.5*dz_lake(j)
                  rsfinvw  = exp( -etav*max(  zin-za(idlak),0. ) )  ! the radiation within surface layer (z<za)
                  rsfoutvw = exp( -etav*max( zout-za(idlak),0. ) )
                  rsfinnw  = exp( -etan*max(  zin-za(idlak),0. ) )  ! the radiation within surface layer (z<za)
                  rsfoutnw = exp( -etan*max( zout-za(idlak),0. ) )
                  rsfinv  = exp( -etaiv*max(  zin-za(idlak),0. ) ) 
                  rsfoutv = exp( -etaiv*max( zout-za(idlak),0. ) )
                  rsfinn  = exp( -etain*max(  zin-za(idlak),0. ) ) 
                  rsfoutn = exp( -etain*max( zout-za(idlak),0. ) )
                  if(j == 1) then
                     phi(j) = sabg *(1-betaprime)* (1. - (1-lake_icefrac(j)) * (rsfoutvw*(1-betaprime) + rsfoutnw*betaprime) &
                              - lake_icefrac(j)  * (rsfoutv * (1.-betaprime) + rsfoutn * betaprime))
                  else
                     phi(j) = sabg*(1-betaprime)* ((1-lake_icefrac(j-1)) * (rsfinvw *(1-betaprime) + rsfinnw *betaprime) &
                              + lake_icefrac(j-1) * (rsfinv * (1.-betaprime) + rsfinn * betaprime) &
                              -  (1-lake_icefrac(j)) * (rsfoutvw*(1-betaprime)+rsfoutnw*betaprime) &
                              - lake_icefrac(j) * (rsfoutv * (1.-betaprime) + rsfoutn * betaprime))
                  end if
                     if (j == nl_lake) phi_soil = (1-lake_icefrac(j)) * sabg *(1-betaprime)*( (1.-betaprime) * rsfoutvw + betaprime * rsfoutnw) &
                                                +lake_icefrac(j)  * sabg * (1-betaprime)*(rsfoutv * (1.-betaprime) + rsfoutn * betaprime)
               end do
            else   ! snow covered, this should be improved upon; Mironov 2002 suggests that SW can penetrate thin ice and may
                  ! cause spring convection.
               phi(:) = 0.
               phi_soil = 0.
            end if 


         ! ELSE IF(DEF_USE_TWOSTREAM)THEN
         !    calday = calendarday(idate)
         !    czen=orb_coszen(calday,dlon,dlat)
         !    coszen=max(czen,0.001)

         !    if (snl == 0) then
         !       do j = 1, nl_lake
         !          etatsa(j)=0.239012444*(1-lake_icefrac(j))+1.5*lake_icefrac(j)
         !       end do
         !       pi=3.14159265359  
         !       tao(1) = 0.
         !       do j = 2, nl_lake+1
         !          dtao(j-1) = etatsa(j-1)*dz_lake(j-1) 
         !          tao(j) = tao(j-1) + dtao(j-1)
         !          print*,'tao',j,tao(j)
         !       end do
         !       w0=0.3568144
         !       g=0.0442765
         !       f=sabgvd/coszen
         !       s1=f*w0/(4*pi)*(1+3*g*(1.0/(3**0.5))*coszen)
         !       s2=f*w0/(4*pi)*(1-3*g*(1.0/(3**0.5))*coszen)
         !       k=((1.0-w0)*(1.0-w0*g)/(1.0/3.0))**0.5
         !       ssa=((1.0-w0)/(1.0-w0*g))**0.5
         !       z1=-((1.0-w0*g)*(s2+s1))/(1.0/3.0)+(s2-s1)/(1.0/(3**0.5)*coszen) 
         !       z2=-((1.0-w0)*(s2-s1))/(1.0/3.0)+(s2+s1)/(1.0/(3**0.5)*coszen) 
         !       saa=z1*(coszen**2.0)/(1.0-(coszen**2.0)*(k**2.0)) !right
         !       sbe=z2*(coszen**2.0)/(1.0-(coszen**2.0)*(k**2.0)) !right
         !       se=(saa+sbe)/2.0
         !       sr=(saa-sbe)/2.0
         !       sv=(1.+ssa)/2.0
         !       su=(1.-ssa)/2.0
         !       tu=exp(-tao(nl_lake+1)/coszen)
         !       !print*,tu,coszen,-tao(nl_lake+1)/coszen,exp(-tao(nl_lake+1)/coszen)
         !       ! print*,'se',se, 'sr', sr,'sv',sv, 'su',su
         !       !F-(0)=sabgvda F+(nl_lake+1)=0
         !       ! ssh=(sv*sr*exp(k*tao(nl_lake+1))-su*se*exp(tu)-sv*sabgvd*exp(k*tao(nl_lake+1))/(2*pi*(1.0/(3.**0.5))))&
         !       !      /(su**2*exp(-k*tao(nl_lake+1))-sv**2*exp(k*tao(nl_lake+1)))
         !       ! ssk=sabgvd/(su*2*pi*(1.0/(3.**0.5)))-ssh*sv/su-sr/su
         !       !F-(0)-F+(0)=sabgvd F+(nl_lake+1)=0
         !       ssh=(sv*((sr-se)/(su-sv))*exp(k*tao(nl_lake+1))-sv*sabgvd*exp(k*tao(nl_lake+1))/(su-sv)/(2*pi*(1.0/(3.**0.5)))&
         !            -se*tu)/(sv*exp(k*tao(nl_lake+1))+su*exp(-k*tao(nl_lake+1)))
         !       ssk=ssh+sabgvd/(su-sv)/(2*pi*(1.0/(3.**0.5)))+(se-sr)/(su-sv)
         !       !F+(0)=sabgvdu F+(nl_lake+1)=0               
         !       ! ssh=(se*exp(k*tao(nl_lake+1))-sabgvdu*exp(k*tao(nl_lake+1))-se*exp(tu))&
         !       !      /(su*exp(k*tao(nl_lake+1))-su*exp(-k*tao(nl_lake+1)))
         !       ! ssk=(sabgvdu-ssh*su-se)/sv
         !       ! print*,'ssh',ssh,'ssk',ssk
         !       r1=(1.0/(1.0/(3.0**0.5)))*(1.0-(w0/2.0)*(1.0+g))
         !       r2=(w0/(2.*(1.0/(3.**0.5))))*(1.0-g)
         !       r3=(1.0/2.0)*(1.0-3.0*g*(1.0/(3**0.5))*coszen)
         !       k2=(r1**2-r2**2)**0.5 !right
         !       sv2=(1.0/2.0)*(1.0+(r1-r2)/k2) !right
         !       su2=(1.0/2.0)*(1.0-(r1-r2)/k2) !right
         !       se2=(r3*(1.0/coszen-r1)-r2*(1-r3))*(coszen**2)*w0*f/(1-(coszen**2)*(k2**2))
         !       sr2=-((1-r3)*(1.0/coszen+r1)+r2*r3)*(coszen**2)*w0*f/(1-(coszen**2)*(k2**2))
         !       !F-(0)=sabgvda F+(nl_lake+1)=0
         !       ! sh2=(sv2*sr2*exp(k2*tao(nl_lake+1))-sv2*sabgvda*exp(k2*tao(nl_lake+1))-su2*se2*tu)&
         !       !     /(su2**2*exp(-k2*tao(nl_lake+1))-sv2**2*exp(k2*tao(nl_lake+1)))
         !       ! sk2=(sabgvda-sv2*sh2-sr2)/su2
         !       !F-(0)-F+(0)=sabgvd F+(nl_lake+1)=0
         !       sh2=((sabgvd+se2-sr2)*sv2*exp(k2*tao(11))/(sv2-su2)-se2*tu)&
         !            /(sv2*exp(k2*tao(11))+su2*exp(-k2*tao(11)))
         !       sk2=sh2+sabgvd/(su2-sv2)+(se2-sr2)/(su2-sv2)

         !       do j = 1, nl_lake+1
         !          ttu=exp(-tao(j)/coszen)
         !          !print*,j,coszen*f*ttu
         !          !fangan1
         !          ! eup(j)=2.0*pi*(1.0/(3**0.5))*(ssk*su*exp(k*tao(j))+ssh*sv*exp(-k*tao(j))+sr*ttu)&
         !          !       +sabgv*ttu&
         !          !       -2.0*pi*(1.0/(3**0.5))*(ssk*sv*exp(k*tao(j))+ssh*su*exp(-k*tao(j))+se*ttu)
         !          ! eupd(j)=2.0*pi*(1.0/(3**0.5))*(ssk*su*exp(k*tao(j))+ssh*sv*exp(-k*tao(j))+sr*ttu)&
         !          !        -2.0*pi*(1.0/(3**0.5))*(ssk*sv*exp(k*tao(j))+ssh*su*exp(-k*tao(j))+se*ttu)
         !          ! !fangan2
         !          eup(j)=su2*sk2*exp(k2*tao(j))+sv2*sh2*exp(-k2*tao(j))+sr2*ttu&
         !                 +sabgv*ttu&
         !                 -(sv2*sk2*exp(k2*tao(j))+su2*sh2*exp(-k2*tao(j))+se2*ttu)
         !          eupd(j)=su2*sk2*exp(k2*tao(j))+sv2*sh2*exp(-k2*tao(j))+sr2*ttu&
         !                 -(sv2*sk2*exp(k2*tao(j))+su2*sh2*exp(-k2*tao(j))+se2*ttu)
         !          eupz(j)=sabgv*ttu                  
         !       end do 
         !       sum = 0
         !       do j = 1,nl_lake+1
         !       !    ttu=exp(-tao(j)/coszen)
         !          ! print*,'edown',j,su2*sk2*exp(k2*tao(j))+sv2*sh2*exp(-k2*tao(j))+sr2*ttu&
         !          !        +sabgv*ttu
         !       end do
         !       do j = 1, nl_lake+1
         !          ! print*,'eup',j,sv2*sk2*exp(k2*tao(j))+su2*sh2*exp(-k2*tao(j))+se2*ttu
         !       end do
         !       do j = 1, nl_lake
         !          phi(j)=eup(j)-eup(j+1)
         !          print*,j,phi(j)
         !          sum=sum+phi(j)
         !          if (j == nl_lake) phi_soil = eup(j+1)
         !       end do
         !       ! print*,'bj',sabgvd,eupd(1),sv2*sk2*exp(k2*tao(11))+su2*sh2*exp(-k2*tao(11))+se2*tu
         !       ! print*,'check energy',sum+sk2*su2*exp(k2*tao(11))+sh2*sv2*exp(-k2*tao(11))&    !sv2*sk2+su2*sh2+se2+
         !       !          +sr2*exp(-tao(11)/coszen)&
         !       !          +sabgv*exp(-tao(11)/coszen),&
         !       !          sabgv+sabgvd,sabgn+sabgnd,sabg*(1-betaprime),(sabgv+sabgvd+sabgn+sabgnd)*(1-betaprime),betaprime
         !    else   ! snow covered, this should be improved upon; Mironov 2002 suggests that SW can penetrate thin ice and may
         !          ! cause spring convection.
         !       phi(:) = 0.
         !       phi_soil = 0.
         !    end if 

            ! if (snl == 0) then      !no snow cover
            !    do j = 1, nl_lake
            !      etatsa1(j)=0.69*(1-lake_icefrac(j))+1.5*lake_icefrac(j)
            !    end do 
            !    pi=3.1415926
            !    uavg=0.5
            !    topc = sabg*(1.-betaprime)/(2*pi*uavg)
            !    print*, topc 
            !    tao(1) = 0.
            !    do j = 2, nl_lake+1
            !       dtao(j-1) = etatsa(j-1)*dz_lake(j-1) 
            !       tao(j) = tao(j-1) + dtao(j-1)
            !    end do
            !   w10=0.595356
            !    daT = (1-w0)**(0.5)/uavg 
            !    rou = (1-uavg*daT)/(1+uavg*daT)
            !    ! print*, T
            !    ! print*, rou
            !    ! print*, tao(nl_lake+1)
            !    D = exp(daT*tao(nl_lake+1))-(rou**2)*exp(-daT*tao(nl_lake+1))
            !    ! print*, D
            !    do j = 1, nl_lake+1
            !       eup(j)=-2*uavg*pi*topc*(1-rou)*(exp(daT*(tao(nl_lake+1)-tao(j)))+rou*exp(-daT*(tao(nl_lake+1)-tao(j))))/D
            !    end do
            !    sum=0
            !    do j = 1, nl_lake
            !       phi(j)=eup(j+1)-eup(j)
            !       sum=sum+phi(j)
            !       if (j == nl_lake) phi_soil = -eup(j+1)
            !    end do 
            !    sum = sum + phi_soil + 2*pi*uavg*topc*rou*(exp(daT*(tao(nl_lake+1)-tao(1)))-exp(-daT*(tao(nl_lake+1)-tao(1))))/D
            !    ! print*, sum,sabg*(1.-betaprime), 2*pi*uavg*topc*rou*(exp(T*(tao(nl_lake+1)-tao(1)))-exp(-T*(tao(nl_lake+1)-tao(1))))/D
            ! else   ! snow covered, this should be improved upon; Mironov 2002 suggests that SW can penetrate thin ice and may
            !       ! cause spring convection.
            !    phi(:) = 0.
            !    phi_soil = 0.
            ! end if 
         ! ENDIF   

         ELSE IF(DEF_USE_TWOSTREAM)THEN
            write(*,*)'twostream' 
            ! zero aerosol input arrays
            DO aer = 1, sno_nbr_aer
               DO i = maxsnl+1, 0
                  mss_cnc_aer_in(i,aer) = 0._r8
               ENDDO
            ENDDO
            if (snl /= 0) then 
               write(*,*)'snowsnl',snl 
            endif 
            write(*,*)'lake_icefrac',lake_icefrac(1)
            ss_alb_wtr_avg_cld(:)=(/0.0556923605439147, 0.000143872453054684,&
            3.88147493245108e-06, 3.51096969783593e-07, 1.64723778561854e-08/) 
            ext_cff_mss_wtr_avg_cld(:) = (/0.117858180005272, 7.08651632540221, 38.4097270605075,&
            296.897220484403, 1521.60471472379/)   ! ice refractive index
            write(*,*)'year',idate(1)
            write(*,*)'day',idate(2)
            calday = calendarday(idate)
            coszen=orb_coszen(calday,dlon,dlat)
            flg_snw_ice = 1
            snwcp_ice   = 0.0       !excess precipitation due to snow capping [kg m-2 s-1]
            do_capsnow  = .false.   !true => DO snow capping
            CALL AerosolMasses( deltim, snl ,do_capsnow ,&
            wice_soisno(:0),wliq_soisno(:0),snwcp_ice      ,snw_rds       ,&

            mss_bcpho     ,mss_bcphi       ,mss_ocpho      ,mss_ocphi     ,&
            mss_dst1      ,mss_dst2        ,mss_dst3       ,mss_dst4      ,&

            mss_cnc_bcphi ,mss_cnc_bcpho   ,mss_cnc_ocphi  ,mss_cnc_ocpho ,&
            mss_cnc_dst1  ,mss_cnc_dst2    ,mss_cnc_dst3   ,mss_cnc_dst4   )

            mss_cnc_aer_in(:,5) = mss_cnc_dst1(:)
            mss_cnc_aer_in(:,6) = mss_cnc_dst2(:)
            mss_cnc_aer_in(:,7) = mss_cnc_dst3(:)
            mss_cnc_aer_in(:,8) = mss_cnc_dst4(:)
            snw_rds_in(:) = nint(snw_rds(:))
            ! write(*,*)'flg_snw_ice1',flg_snw_ice
         !    coszen=max(czen,0.001)
         !--------------------------Define all layer---------------------------------------- 
            srf_lyr = 0
            nbr_lyr = nl_lake 
            do j = 1, nl_lake 
               if (lake_icefrac(j)>0.and.lake_icefrac(j)<1) then !Mixed ice-water layer
                  nbr_lyr = nl_lake + 1 
                  iwin = j ! The interface between ice and water
               endif 
            end do 
            ! lyr_typ(:) = 0. ! 1=snow, 2=ice, 3=water
            ! allocate(lyr_typ(snl:nbr_lyr))
            ! lyr_typ(snl:0) = 1 !snow layers
            ! write (*,*) 'layer1'
            dziw(:) = 0
            if (nbr_lyr==10.) then ! no mixed ice-water layer
               do j = 1,nbr_lyr
                  if (lake_icefrac(j)==1) then 
                     dziw(j) = dz_lake(j) ! depth of ice or water layer
                  elseif (lake_icefrac(j)==0.) then 
                     dziw(j) = dz_lake(j)! depth of ice or water layer
                  endif
               end do
            elseif (nbr_lyr==11.) then ! mixed ice-water layer
               do j = 1, nbr_lyr
                  if( j < iwin ) then ! above interface between ice and water
                     dziw(j) = dz_lake(j)
                  elseif ( j == iwin) then ! interface between ice and water
                     dziw(j) = dz_lake(j)*lake_icefrac(j)
                     dziw(j+1) = dz_lake(j)*(1-lake_icefrac(j))
                  elseif (j > iwin+1) then
                     dziw(j) = dz_lake(j-1)
                  endif 
               enddo 
            endif 
            ! do i = 1, nbr_lyr
            !    write (*,*) i,'dziw',dziw(i)
            ! enddo
            ! ziw(1) = 0 
            ziw(:) = 0
            ziwsum = 0
            do i = 1, nbr_lyr
               ziwsum = ziwsum + dziw(i)
               ziw(i) = ziwsum
            enddo
            do i = 1, nbr_lyr
               if(ziw(i)<=za(idlak)) then 
                  srf_lyr = i 
               endif 
            enddo 
            write(*,*)'srf_lyr',srf_lyr
            nl_ice = 0  ! number of ice layers
            nl_wat = nl_lake !number of water layers
            do j = 1, nl_lake
            ! write(*,*)'lake_icefrac' ,lake_icefrac(j)
               if (lake_icefrac(j) == 1.) then
                  nl_ice = j
                  nl_wat = nl_lake-j
               elseif (lake_icefrac(j)>0.and.lake_icefrac(j)<1.) then 
                  nl_ice = j 
                  nl_wat = nl_lake-j+1 
               endif
            end do 
            ! if(snl < 0 ) then 
            !    write(*,*)'snl',snl
            ! endif
            if(nl_ice > 0 ) then 
               write(*,*)'nl_ice',nl_ice
            endif
            ! write(*,*)'t_grnd',t_grnd ,' t_lake(1)', t_lake(1) , tfrz
         !-------------------------------find kfrsnl--------------------------------
            kfrsnl1 = nbr_lyr+5 !between the last snow layer and the first ice layer or air and ice
            kfrsnl2 = nbr_lyr+5 !between air and water 
            kfrsnl3 = nbr_lyr+5 !between the last ice layer and the first water layer
            kfrsnl4 = nbr_lyr+5 !between air and ice layer

            if (lake_icefrac(1)>0.and.snl==0) then 
               kfrsnl4 = 1
            endif
            if (lake_icefrac(1)>0.and.snl<0) then 
               kfrsnl1 = 1
            endif
            if (lake_icefrac(1)==0.and.snl == 0) then 
               kfrsnl2 = 1
            endif       
            do j = 1, nl_lake-1 
               if (lake_icefrac(j)>0.and.lake_icefrac(j+1)==0) then
                  kfrsnl3 = j+1
               endif
            end do 
            write(*,*)'kfrsnl1',kfrsnl1,'kfrsnl2',kfrsnl2,'kfrsnl3',kfrsnl3,'kfrsnl4',kfrsnl4
            ! Define constants
            pi = SHR_CONST_PI
            nint_snw_rds_min = nint(snw_rds_min) 
            ! nint_snw_rds_max = nint(snw_rds_max) !-----06/13--------

            ! always use Delta approximation for snow
            DELTA = 1

            !Gaussian integration angle and coefficients for diffuse radiation
            difgauspt(1:8)     & ! gaussian angles (radians)
                  = (/ 0.9894009_r8,  0.9445750_r8, &
                     0.8656312_r8,  0.7554044_r8, &
                     0.6178762_r8,  0.4580168_r8, &
                     0.2816036_r8,  0.0950125_r8/)
            difgauswt(1:8)     & ! gaussian weights
                  = (/ 0.0271525_r8,  0.0622535_r8, &
                     0.0951585_r8,  0.1246290_r8, &
                     0.1495960_r8,  0.1691565_r8, &
                     0.1826034_r8,  0.1894506_r8/)
            ! dif_a(1:5)     &
            !    = (/ 0.0633402270380342_r8,  0.06189314981880321_r8, &
            !          0.0611905910736522_r8,  0.0605842642044964_r8, &
            !          0.0617648953124473_r8/)
         
            ! dif_b(1:5)     &
            !    = (/  0.432870505704913_r8,  0.424519192762872_r8, &
            !          0.420382469844439_r8,  0.416853617451133_r8, &
            !          0.411700917011931_r8/)

            snw_shp_lcl(:) = snow_shape_sphere
            snw_fs_lcl(:)  = 0._r8
            snw_ar_lcl(:)  = 0._r8
            atm_type_index = atm_type_mid_latitude_winter

            ! Define snow grain shape
            if (trim(snow_shape) == 'sphere') then
                  snw_shp_lcl(:) = snow_shape_sphere
            elseif (trim(snow_shape) == 'spheroid') then
                  snw_shp_lcl(:) = snow_shape_spheroid
            elseif (trim(snow_shape) == 'hexagonal_plate') then
                  snw_shp_lcl(:) = snow_shape_hexagonal_plate
            elseif (trim(snow_shape) == 'koch_snowflake') then
                  snw_shp_lcl(:) = snow_shape_koch_snowflake
            ! else
            !       IF (p_is_master) THEN
            !          write(iulog,*) "snow_shape = ", snow_shape
            !          call abort
            !       ENDIF
            endif

            ! Define atmospheric type
            if (trim(snicar_atm_type) == 'default') then
                  atm_type_index = atm_type_default
            elseif (trim(snicar_atm_type) == 'mid-latitude_winter') then
                  atm_type_index = atm_type_mid_latitude_winter
            elseif (trim(snicar_atm_type) == 'mid-latitude_summer') then
                  atm_type_index = atm_type_mid_latitude_summer
            elseif (trim(snicar_atm_type) == 'sub-Arctic_winter') then
                  atm_type_index = atm_type_sub_Arctic_winter
            elseif (trim(snicar_atm_type) == 'sub-Arctic_summer') then
                  atm_type_index = atm_type_sub_Arctic_summer
            elseif (trim(snicar_atm_type) == 'summit_Greenland') then
                  atm_type_index = atm_type_summit_Greenland
            elseif (trim(snicar_atm_type) == 'high_mountain') then
                  atm_type_index = atm_type_high_mountain
            ! else
            !       IF (p_is_master) THEN
            !          write(iulog,*) "snicar_atm_type = ", snicar_atm_type
            !          call abort
            !       ENDIF
            endif

            ! Define ice parameter
            ice_ri = 3      ! ice refractive index dataset 
            do i = 1, nbr_lyr
               if(nbr_lyr==10) then 
                  rho_snw(1:nbr_lyr) = rhow(:)  ! snow density for each layer (units: kg/m3)
               elseif (i <= nl_ice) then 
                  rho_snw(i) = rhow(i)
               elseif (i > nl_ice) then
                  rho_snw(i) = rho_snw(i-1)
               endif              
            enddo
            ! write(*,*)'rho_snw'
            ! rds_snw(1:nl_ice) = 750  ! snow grain size for each snow layer or bubble radius for ice
            ! write(*,*)'rds_snw'
            
            snw_rds_lcl(1:nbr_lyr)    = 750
            ice_density_wgted = 916.999  !adjust
            write(*,*)'ice_properity start'
            if (idate(1)==2015.and.idate(2)>342) then 
               snw_rds_lcl(1:nbr_lyr) = nint(((30.0-idate(2)+343.0)/30.0)*974.66666667+ &
                                       ((idate(2)-343.0)/30.0)*1040.66666667)
               ice_density_wgted = ((30.0-idate(2)+342.0)/30.0)*850.+ &
                                    ((idate(2)-342.0)/30.0)*852.
            elseif (idate(1)==2016.and.idate(2)<9) then 
               snw_rds_lcl(1:nbr_lyr) = nint(((8.0-idate(2))/30.0)*974.66666667+ &
                                       ((22.0+idate(2))/30.0)*1040.66666667)
               ice_density_wgted = (((8.0-idate(2))/30.0))*850.+ &
                                    ((22.0+idate(2))/30.0)*852.
            elseif (idate(1)==2016.and.idate(2)>8.and.idate(2)<38) then 
               snw_rds_lcl(1:nbr_lyr) = nint(((30.0-idate(2)+8.0)/30.0)*1040.66666667+ &
                                       ((idate(2)-8.0)/30.0)*824.)
               ice_density_wgted = (((30.0-idate(2)+8.0)/30.0))*852.+ &
                                    ((idate(2)-8.0)/30.0)*845.5 
            elseif (idate(1)==2016.and.idate(2)>37.and.idate(2)<68) then 
               snw_rds_lcl(1:nbr_lyr) = nint(((30.0-idate(2)+37.0)/30.0)*824.+ &
                                       ((idate(2)-37.0)/30.0)*907.33333333)
               ice_density_wgted = (((30.0-idate(2)+37.0)/30.0))*845.5+ &
                                    ((idate(2)-37.0)/30.0)*794.5   
            elseif (idate(1)==2016.and.idate(2)>67.and.idate(2)<=98) then 
               snw_rds_lcl(1:nbr_lyr) = nint(((30.0-idate(2)+67.0)/30.0)*907.33333333+ &
                                       ((idate(2)-67.0)/30.0)*913.33333333)
               ice_density_wgted = (((30.0-idate(2)+67.0)/30.0))*794.5+ &
                                    ((idate(2)-67.0)/30.0)*803.75                             
            endif
            write(*,*)'ice_properity end',snw_rds_lcl(1),ice_density_wgted
            ! rho_air = 1.025 ! density of air [kg / m3] 
            !-------correct dimension later--------
            do i=maxsnl+1,nbr_lyr+1,1
               flx_abs_lcl(:,:)   = 0._r8
               flx_sum2(i)        = 0._r8
            enddo
            do i=maxsnl+1,nbr_lyr+1,1
               flx_abs(i,:) = 0._r8
            enddo 
            do i=maxsnl+1,nbr_lyr,1
               mss_cnc_aer_lcl(i,:) = 0._r8 
            enddo
            ! write(*,*)'flg_snw_ice2',flg_snw_ice
            ! set snow/ice mass to be used for RT:
            do flg_slr_in = 1, 2, 1
               ! write(*,*)'flg_slr_in',flg_slr_in
               if (flg_snw_ice == 1) then
                  h2osno_lcl = h2osno
               else
                  h2osno_lcl = wice_soisno(0)
               endif
               ! write(*,*) 'before if'
               ! Qualifier for computing snow RT:
               !  1) sunlight from atmosphere model
               !  2) minimum amount of snow on ground.
               !     Otherwise, set snow albedo to zero
               ! if ((coszen > 0._r8) .and. (h2osno_lcl > min_snw) ) then

               ! Set variables specific to ELM
               if (flg_snw_ice == 1) then
                  ! If there is snow, but zero snow layers, we must create a layer locally.
                  ! This layer is presumed to have the fresh snow effective radius.
                  if (snl > -1) then
                     flg_nosnl         =  1
                     snl_lcl           =  0
                     h2osno_ice_lcl(0) =  h2osno_lcl
                     h2osno_liq_lcl(0) =  0._r8
                     snw_rds_lcl(0)    =  nint_snw_rds_min
                  else
                     flg_nosnl         =  0
                     snl_lcl           =  snl
                     h2osno_liq_lcl(maxsnl+1:0) =  h2osno_liq(maxsnl+1:0)
                     h2osno_ice_lcl(maxsnl+1:0) =  h2osno_ice(maxsnl+1:0)
                     snw_rds_lcl(maxsnl+1:0)    =  snw_rds_in(maxsnl+1:0)
                     ! do i = 1, nbr_lyr
                     !    if (nbr_lyr == 10) then 
                     !       h2osno_liq_lcl(1:10) =  wliq_soisno(:)
                     !       h2osno_ice_lcl(1:10) =  wice_soisno(:)
                     !       ! snw_rds_lcl(1:10)    =  snw_rds(:)
                     !    elseif (nbr_lyr == 11.and.0 <i <= nl_ice) then 
                     !       h2osno_liq_lcl(i) = wliq_soisno(i)
                     !       h2osno_ice_lcl(i) =  wice_soisno(i)
                     !       ! snw_rds_lcl(i)    =  snw_rds(i)
                     !    elseif (nbr_lyr == 11.and.i > nl_ice) then
                     !       h2osno_liq_lcl(i) = wliq_soisno(i-1)
                     !       h2osno_ice_lcl(i) =  wice_soisno(i-1)
                     !       ! snw_rds_lcl(i)    =  snw_rds(i-1) 
                     !    endif
                     ! enddo                         
                  endif
                  ! write(*,*)'flg_snw_ice == 1'

                  snl_btm   = 0
                  snl_top   = snl_lcl+1
                  ! write(*,*)'snl_top_fir',snl_top,snl_lcl
                  ! Set variables specific to CSIM
               else
                  flg_nosnl         = 0
                  snl_lcl           = -1
                  ! do i = 1, nbr_lyr
                  !    if (nbr_lyr == 10) then 
                  !       h2osno_liq_lcl(1:10) =  wliq_soisno(:)
                  !       h2osno_ice_lcl(1:10) =  wice_soisno(:)
                  !       ! snw_rds_lcl(1:10)    =  snw_rds(:)
                  !    elseif (i <= nl_ice) then 
                  !       h2osno_liq_lcl(i) = wliq_soisno(i)
                  !       h2osno_ice_lcl(i) =  wice_soisno(i)
                  !       ! write(*,*)'h2osno_liq_water'
                  !       ! snw_rds_lcl(i)    =  snw_rds(i)
                  !    elseif (i > nl_ice) then
                  !       h2osno_liq_lcl(i) = wliq_soisno(i-1)
                  !       h2osno_ice_lcl(i) =  wice_soisno(i-1)
                  !       ! snw_rds_lcl(i)    =  snw_rds(i-1)  
                  !       ! write(*,*)'h2osno_liq_water'
                  !    endif
                  ! enddo     
                  snl_btm           = 0
                  snl_top           = 0
               endif ! end if flg_snw_ice == 1
               ! write(*,*) 'end if flg_snw_ice == 1'
               do i = snl_top,nbr_lyr,1
                  if (i <= 0) then 
                     if (snw_rds_lcl(i)<=30) then 
                        snw_rds_idx(i) = 1
                     elseif(snw_rds_lcl(i)<=1500) then 
                        snw_rds_idx(i) = snw_rds_lcl(i)-29
                     elseif (snw_rds_lcl(i)<=20000)then
                        snw_rds_idx(i) = int((snw_rds_lcl(i)-1500)/250)+1471
                     else
                        snw_rds_idx(i) = 1509
                     endif 
                  else 
                     if (snw_rds_lcl(i)<=10) then 
                        snw_rds_idx(i) = 1
                     elseif(snw_rds_lcl(i)<=1500) then 
                        snw_rds_idx(i) = snw_rds_lcl(i)-9
                     elseif (snw_rds_lcl(i)<=20000)then 
                        snw_rds_idx(i) = int((snw_rds_lcl(i)-1500)/250)+1491
                     else 
                        snw_rds_idx(i) = 1566
                     endif                      
                  endif
                  write(*,*)'snw_rds_idx',i,snw_rds_idx(i)
               enddo
               
#ifdef MODAL_AER
               !mgf++
               !
               ! Assume fixed BC effective radii of 100nm. This is close to
               ! the effective radius of 95nm (number median radius of
               ! 40nm) assumed for freshly-emitted BC in MAM.  Future
               ! implementations may prognose the BC effective radius in
               ! snow.
               rds_bcint_lcl(:)  =  100._r8
               rds_bcext_lcl(:)  =  100._r8
               !mgf--
#endif
               ! write(*,*)'aer'
               ! Set local aerosol array
               do j=1,sno_nbr_aer
                  mss_cnc_aer_lcl(maxsnl+1:0,j) = mss_cnc_aer_in(:,j)
                  ! write(*,*),'mss_cnc_aer_lcl',j,mss_cnc_aer_lcl(maxsnl+1:0,j)
               enddo
               ! write(*,*) 'Set local aerosol array'
               ! Set spectral underlying surface albedos to their corresponding VIS or NIR albedos
               albsfc_lcl(1)                       = albsfc(1)
               albsfc_lcl(nir_bnd_bgn:nir_bnd_end) = albsfc(2)

                  ! Error check for snow grain size:
                  ! IF (p_is_master) THEN
                  !    do i=snl_top,snl_btm,1
                  !    if ((snw_rds_lcl(i) < snw_rds_min_tbl) .or. (snw_rds_lcl(i) > snw_rds_max_tbl)) then
                  !       write (iulog,*) "SNICAR ERROR: snow grain radius of ", snw_rds_lcl(i), " out of bounds."
                  !       write (iulog,*) "flg_snw_ice= ", flg_snw_ice
                  !       write (iulog,*) " level: ", i, " snl(c)= ", snl_lcl
                  !       write (iulog,*) "h2osno(c)= ", h2osno_lcl
                  !       call abort
                  !    endif
                  !    enddo
                  ! ENDIF
               ! Incident flux weighting parameters
               !  - sum of all VIS bands must equal 1
               !  - sum of all NIR bands must equal 1
               !
               ! Spectral bands (5-band case)
               !  Band 1: 0.3-0.7um (VIS)
               !  Band 2: 0.7-1.0um (NIR)
               !  Band 3: 1.0-1.2um (NIR)
               !  Band 4: 1.2-1.5um (NIR)
               !  Band 5: 1.5-5.0um (NIR)
               !
               ! The following weights are appropriate for surface-incident flux in a mid-latitude winter atmosphere
               !
               ! 3-band weights
               write(*,*)atm_type_index
               if (numrad_snw==3) then
                  ! Direct:
                  if (flg_slr_in == 1) then
                     flx_wgt(1) = 1._r8
                     flx_wgt(2) = 0.66628670195247_r8
                     flx_wgt(3) = 0.33371329804753_r8
                  ! Diffuse:
                  elseif (flg_slr_in == 2) then
                     flx_wgt(1) = 1._r8
                     flx_wgt(2) = 0.77887652162877_r8
                     flx_wgt(3) = 0.22112347837123_r8
                  endif

                  ! 5-band weights
               elseif(numrad_snw==5) then
                  ! Direct:
                  if (flg_slr_in == 1) then
                     if (atm_type_index == atm_type_default) then
                        flx_wgt(1) = 1._r8
                        flx_wgt(2) = 0.49352158521175_r8
                        flx_wgt(3) = 0.18099494230665_r8
                        flx_wgt(4) = 0.12094898498813_r8
                        flx_wgt(5) = 0.20453448749347_r8
                     else
                        slr_zen = nint(acos(coszen) * 180._r8 / pi)
                        if (slr_zen>89) then
                              slr_zen = 89
                        endif
                        ! write(*,*)'in'
                        flx_wgt(1) = 1._r8
                        flx_wgt(2) = flx_wgt_dir(atm_type_index, slr_zen+1, 2)
                        flx_wgt(3) = flx_wgt_dir(atm_type_index, slr_zen+1, 3)
                        flx_wgt(4) = flx_wgt_dir(atm_type_index, slr_zen+1, 4)
                        flx_wgt(5) = flx_wgt_dir(atm_type_index, slr_zen+1, 5)
                        write(*,*)flx_wgt(2)
                     endif

                  ! Diffuse:
                  elseif (flg_slr_in == 2) then
                     if  (atm_type_index == atm_type_default) then
                        flx_wgt(1) = 1._r8
                        flx_wgt(2) = 0.58581507618433_r8
                        flx_wgt(3) = 0.20156903770812_r8
                        flx_wgt(4) = 0.10917889346386_r8
                        flx_wgt(5) = 0.10343699264369_r8
                     else
                        flx_wgt(1) = 1._r8
                        flx_wgt(2) = flx_wgt_dif(atm_type_index, 2)
                        flx_wgt(3) = flx_wgt_dif(atm_type_index, 3)
                        flx_wgt(4) = flx_wgt_dif(atm_type_index, 4)
                        flx_wgt(5) = flx_wgt_dif(atm_type_index, 5)
                     endif
                  endif
               endif ! end if numrad_snw 
               ! write(*,*) 'end if numrad_snw '
               ! Loop over snow spectral bands

               exp_min = exp(-argmax)
               do bnd_idx = 1,numrad_snw
               ! write(*,*)'flg_slr_in', flg_slr_in
                  ! mss_cnc_aer_lcl(:,:) = 0._r8
                  ! if ( (lnd_ice == 1) ) then
                  !    mss_cnc_aer_lcl(1:10,1) = ice_bc_conc_wgted
                  !    if (bnd_idx > 1) then
                  !       mss_cnc_aer_lcl(1:10,1) = 0._r8
                  !    endif
                  ! endif
               ! note that we can remove flg_dover since this algorithm is
               ! stable for mu_not > 0.01

               ! mu_not is cosine solar zenith angle above the fresnel level; make
               ! sure mu_not is large enough for stable and meaningful radiation
               ! solution: .01 is like sun just touching horizon with its lower edge
               ! equivalent to mu0 in sea-ice shortwave model ice_shortwave.F90
                  mu_not = max(coszen, cp01)
                  ! Pre-emptive error handling: aerosols can reap havoc on these absorptive bands.
                  ! Since extremely high soot concentrations have a negligible effect on these bands, zero them.
                  if ( (numrad_snw == 5).and.((bnd_idx == 5).or.(bnd_idx == 4)) ) then
                     mss_cnc_aer_lcl(:,:) = 0._r8
                  endif
                  if ( (numrad_snw == 3).and.(bnd_idx == 3) ) then
                     mss_cnc_aer_lcl(:,:) = 0._r8
                  endif
                  ! write(*,*)'aer'
                  ! Define local Mie parameters based on snow grain size and aerosol species,
                  !  retrieved from a lookup table.
                  if (flg_slr_in == 1) then
                     refindxwat_im(bnd_idx) = refindxwat_im_clr(bnd_idx,slr_zen+1)
                     refindxwat_re(bnd_idx) = refindxwat_re_clr(bnd_idx,slr_zen+1)
                     refindx_im(bnd_idx)    = refindx_im_clr(bnd_idx,slr_zen+1)
                     refindx_re(bnd_idx)    = refindx_re_clr(bnd_idx,slr_zen+1)
                     ! write(*,*)'read_rfidx'
                     do i=snl_top,nbr_lyr,1
                        ! write(*,*)'nl_ice',nl_ice
                        ! write(*,*) 'i', i
                        if (i <=0 ) then     !snow
                           ss_alb_snw_lcl(i)      = ss_alb_snw_avg_clr(bnd_idx,slr_zen+1,snw_rds_idx(i))
                           asm_prm_snw_lcl(i)     = asm_prm_snw_avg_clr(bnd_idx,slr_zen+1,snw_rds_idx(i))
                           ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_avg_clr(bnd_idx,slr_zen+1,snw_rds_idx(i))
                           ! write(*,*) 'exist snow dir'
                           ! write(*,*)i,'ss_alb_snw_lcl',ss_alb_snw_lcl(i)
                           ! write(*,*)i,'asm_prm_snw_lcl',asm_prm_snw_lcl(i)
                           ! write(*,*)i,'asm_prm_snw_lcl',asm_prm_snw_lcl(i)

                        else if (nl_ice>0.and.i<=nl_ice) then     !ice
                           sca_cff_vlm_airbbl_lcl(i) = sca_cff_vlm_avg_clr(bnd_idx,slr_zen+1,snw_rds_idx(i)) ! +CAW
                           ! write(*,*)'sca_cff_vlm_airbbl_lcl', sca_cff_vlm_airbbl_lcl(i)
                           asm_prm_snw_lcl(i)        = asm_prm_ice_avg_clr(bnd_idx,slr_zen+1,snw_rds_idx(i))     ! +CAW
                           ! write(*,*)'asm_prm_snw_lcl', asm_prm_snw_lcl(i)
                           abs_cff_mss_ice_lcl(i)    = abs_cff_mss_avg_clr(bnd_idx,slr_zen+1,snw_rds_idx(i))            ! +CAW
                           ! vlm_frac_air(i)           = ( ice_density_wgted -rho_snw(i)  ) / ice_density_wgted         ! +CAW
                           ! ext_cff_mss_snw_lcl(i)    = ((sca_cff_vlm_airbbl_lcl(i) * &
                           !          vlm_frac_air(i)) /rho_snw(i)) + abs_cff_mss_ice_lcl(i) ! +CAW
                           ! ss_alb_snw_lcl(i)         = ((sca_cff_vlm_airbbl_lcl(i) * &
                           !          vlm_frac_air(i)) /rho_snw(i)) / ext_cff_mss_snw_lcl(i) ! +CAW 
                           
                           vlm_frac_air(i)           = ( 917.0 - ice_density_wgted) / 917.0         ! +CAW
                           ext_cff_mss_snw_lcl(i)    = ((sca_cff_vlm_airbbl_lcl(i) * &
                                    vlm_frac_air(i)) /ice_density_wgted) + abs_cff_mss_ice_lcl(i) ! +CAW
                           ss_alb_snw_lcl(i)         = ((sca_cff_vlm_airbbl_lcl(i) * &
                                    vlm_frac_air(i)) /ice_density_wgted) / ext_cff_mss_snw_lcl(i) ! +CAW 
                           ! write(*,*)i,'ice_density_wgted -rho_snw' ,ice_density_wgted -rho_snw(i)
                           ! write(*,*)i, 'vlm_frac_air' , vlm_frac_air(i),'ext_cff_mss_snw_lcl'  ,ext_cff_mss_snw_lcl(i),&
                           !               'ss_alb_snw_lcl' , ss_alb_snw_lcl(i)                   
                           ! write(*,*) 'exist ice dir'
                        else    !water
                           ! ss_alb_snw_lcl(i)      = 0.002
                           ! asm_prm_snw_lcl(i)     = 0.0442765
                           ! ext_cff_mss_snw_lcl(i) = 0.12/rho_snw(i)
                           ss_alb_snw_lcl(i)      = ss_alb_wtr_avg_clr(bnd_idx,slr_zen+1)
                           if(bnd_idx == 1) then 
                              ext_cff_mss_wtr_avg_clr(bnd_idx,:) = 0.3
                              ss_alb_snw_lcl(bnd_idx) = 0.001/0.3
                           endif 
                           asm_prm_snw_lcl(i)     = 0.0442765
                           ext_cff_mss_snw_lcl(i) = ext_cff_mss_wtr_avg_clr(bnd_idx,slr_zen+1)/rho_snw(i)
                           ! write(*,*) 'water dir'
                        endif                            
                     enddo
                     ! write(*,*)'read dir iops'
                  elseif (flg_slr_in == 2) then
                     refindx_im(bnd_idx)    = refindx_im_cld(bnd_idx)
                     refindx_re(bnd_idx)    = refindx_re_cld(bnd_idx)
                     refindxwat_im(bnd_idx) = refindxwat_im_cld(bnd_idx)
                     refindxwat_re(bnd_idx) = refindxwat_re_cld(bnd_idx)
                     ! write(*,*)'diffuse'
                     do i=snl_top,nbr_lyr,1
                        if (i <= 0) then
                           ss_alb_snw_lcl(i)      = ss_alb_snw_avg(bnd_idx,snw_rds_idx(i))
                           asm_prm_snw_lcl(i)     = asm_prm_snw_avg(bnd_idx,snw_rds_idx(i))
                           ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_avg(bnd_idx,snw_rds_idx(i))
                           ! write(*,*) 'exist snow dif'
                           ! write(*,*)i,'ss_alb_snw_lcl',ss_alb_snw_lcl(i)
                           ! write(*,*)i,'asm_prm_snw_lcl',asm_prm_snw_lcl(i)
                           ! write(*,*)i,'asm_prm_snw_lcl',asm_prm_snw_lcl(i)
                        else if (i> 0.and.i<=nl_ice) then                               
                           sca_cff_vlm_airbbl_lcl(i) = sca_cff_vlm_avg(bnd_idx,snw_rds_idx(i)) ! +CAW
                           asm_prm_snw_lcl(i)        = asm_prm_ice_avg(bnd_idx,snw_rds_idx(i))     ! +CAW
                           abs_cff_mss_ice_lcl(i)    = abs_cff_mss_avg(bnd_idx,snw_rds_idx(i))            ! +CAW
                           ! vlm_frac_air(i)           = (ice_density_wgted-rho_snw(i) ) / ice_density_wgted         ! +CAW
                           ! ext_cff_mss_snw_lcl(i)    = ((sca_cff_vlm_airbbl_lcl(i) * &
                           !          vlm_frac_air(i)) /rho_snw(i)) + abs_cff_mss_ice_lcl(i) ! +CAW
                           ! ss_alb_snw_lcl(i)         = ((sca_cff_vlm_airbbl_lcl(i) * &
                           !          vlm_frac_air(i)) /rho_snw(i)) / ext_cff_mss_snw_lcl(i) ! +CAW 
                           vlm_frac_air(i)           = (917.0-ice_density_wgted) / 917.0         ! +CAW
                           ext_cff_mss_snw_lcl(i)    = ((sca_cff_vlm_airbbl_lcl(i) * &
                                    vlm_frac_air(i)) /ice_density_wgted) + abs_cff_mss_ice_lcl(i) ! +CAW
                           ss_alb_snw_lcl(i)         = ((sca_cff_vlm_airbbl_lcl(i) * &
                                    vlm_frac_air(i)) /ice_density_wgted) / ext_cff_mss_snw_lcl(i) ! +CAW 
                           ! write(*,*)i,'ice_density_wgted -rho_snw' ,ice_density_wgted -rho_snw(i)
                           ! write(*,*)i, 'vlm_frac_air' , vlm_frac_air(i),'ext_cff_mss_snw_lcl'  ,ext_cff_mss_snw_lcl(i),&
                           !                'ss_alb_snw_lcl' , ss_alb_snw_lcl(i)   
                           ! write(*,*) 'exist ice dif'
                        else    !water
                           ! ss_alb_snw_lcl(i)      = 0.002
                           ! asm_prm_snw_lcl(i)     = 0.0442765
                           ! ext_cff_mss_snw_lcl(i) = 0.12/rho_snw(i)
                           ss_alb_snw_lcl(i)      = ss_alb_wtr_avg_cld(bnd_idx)
                           if(bnd_idx == 1) then 
                              ext_cff_mss_wtr_avg_cld(bnd_idx) = 0.3
                              ss_alb_snw_lcl(bnd_idx) = 0.001/0.3
                           endif 
                           asm_prm_snw_lcl(i)     = 0.0442765
                           ext_cff_mss_snw_lcl(i) = ext_cff_mss_wtr_avg_cld(bnd_idx)/rho_snw(i)   
                        endif                        
                        
                     enddo
                     ! write(*,*)'read dif iops'
                  endif  
                  ! write(*,*)'read iops'
                  temp1      = (refindx_re(bnd_idx)*refindx_re(bnd_idx))-(refindx_im(bnd_idx)*refindx_im(bnd_idx))+(sin(acos(coszen))*sin(acos(coszen)))
                  temp2      = (refindx_re(bnd_idx)*refindx_re(bnd_idx))-(refindx_im(bnd_idx)*refindx_im(bnd_idx))-(sin(acos(coszen))*sin(acos(coszen)))
                  nreal      = (sqrt(2._r8)/2._r8)* sqrt(temp1 + sqrt( (temp2*temp2) + (4*refindx_re(bnd_idx)*refindx_re(bnd_idx)*refindx_im(bnd_idx)*refindx_im(bnd_idx)) ) )
                  ! write(*,*)bnd_idx,'nreal_ice',nreal
                  temp1      = (refindxwat_re(bnd_idx)*refindxwat_re(bnd_idx))-(refindxwat_im(bnd_idx)*refindxwat_im(bnd_idx))+(sin(acos(coszen))*sin(acos(coszen)))
                  temp2      = (refindxwat_re(bnd_idx)*refindxwat_re(bnd_idx))-(refindxwat_im(bnd_idx)*refindxwat_im(bnd_idx))-(sin(acos(coszen))*sin(acos(coszen)))
                  nreal_wtr  = (sqrt(2._r8)/2._r8)* sqrt(temp1 + sqrt( (temp2*temp2) + (4*refindxwat_re(bnd_idx)*refindxwat_re(bnd_idx)*refindxwat_im(bnd_idx)*refindxwat_im(bnd_idx)) ) )
                  ! write(*,*)bnd_idx,'nreal_wtr',nreal_wtr  


            ! Calculate the asymetry factors under different snow grain shapes
                  if (snl_top < 0 ) then 
                     do i=snl_top,0,1
                        if(snw_shp_lcl(i) == snow_shape_spheroid) then ! spheroid
                           diam_ice = 2._r8*snw_rds_lcl(i)
                           if(snw_fs_lcl(i) == 0._r8) then
                              fs_sphd = 0.929_r8
                           else
                              fs_sphd = snw_fs_lcl(i)
                           endif
                           fs_hex = 0.788_r8
                           if(snw_ar_lcl(i) == 0._r8) then
                              AR_tmp = 0.5_r8
                           else
                              AR_tmp = snw_ar_lcl(i)
                           endif
                           g_ice_Cg_tmp = g_b0 * ((fs_sphd/fs_hex)**g_b1) * (diam_ice**g_b2)
                           gg_ice_F07_tmp = g_F07_c0 + g_F07_c1 * AR_tmp + g_F07_c2 * (AR_tmp**2)
                        elseif(snw_shp_lcl(i) == snow_shape_hexagonal_plate) then ! hexagonal plate
                           diam_ice = 2._r8*snw_rds_lcl(i)
                           if(snw_fs_lcl(i) == 0._r8) then
                              fs_hex0 = 0.788_r8
                           else
                              fs_hex0 = snw_fs_lcl(i)
                           endif
                           fs_hex = 0.788_r8
                           if(snw_ar_lcl(i) == 0._r8) then
                              AR_tmp = 2.5_r8
                           else
                              AR_tmp = snw_ar_lcl(i)
                           endif
                           g_ice_Cg_tmp = g_b0 * ((fs_hex0/fs_hex)**g_b1) * (diam_ice**g_b2)
                           gg_ice_F07_tmp = g_F07_p0 + g_F07_p1 * log(AR_tmp) + g_F07_p2 * ((log(AR_tmp))**2)
                        elseif(snw_shp_lcl(i) == snow_shape_koch_snowflake) then ! Koch snowflake
                           diam_ice = 2._r8 * snw_rds_lcl(i) /0.544_r8
                           if(snw_fs_lcl(i) == 0._r8) then
                              fs_koch = 0.712_r8
                           else
                              fs_koch = snw_fs_lcl(i)
                           endif
                           fs_hex = 0.788_r8
                           if(snw_ar_lcl(i) == 0._r8) then
                              AR_tmp = 2.5_r8
                           else
                              AR_tmp = snw_ar_lcl(i)
                           endif
                           g_ice_Cg_tmp = g_b0 * ((fs_koch/fs_hex)**g_b1) * (diam_ice**g_b2)
                           gg_ice_F07_tmp = g_F07_p0 + g_F07_p1 * log(AR_tmp) + g_F07_p2 * ((log(AR_tmp))**2)
                        endif

                        ! Linear interpolation for calculating the asymetry factor at band_idx.
                        if(snw_shp_lcl(i) > 1) then
                           if(bnd_idx == 1) then
                              g_Cg_intp   = (g_ice_Cg_tmp(2)-g_ice_Cg_tmp(1))/(1.055_r8-0.475_r8)*(0.5_r8-0.475_r8) +g_ice_Cg_tmp(1)
                              gg_F07_intp = (gg_ice_F07_tmp(2)-gg_ice_F07_tmp(1))/(1.055_r8-0.475_r8)*(0.5_r8-0.475_r8)+gg_ice_F07_tmp(1)
                           elseif(bnd_idx == 2) then
                              g_Cg_intp   = (g_ice_Cg_tmp(2)-g_ice_Cg_tmp(1))/(1.055_r8-0.475_r8)*(0.85_r8-0.475_r8)+g_ice_Cg_tmp(1)
                              gg_F07_intp = (gg_ice_F07_tmp(2)-gg_ice_F07_tmp(1))/(1.055_r8-0.475_r8)*(0.85_r8-0.475_r8)+gg_ice_F07_tmp(1)
                           elseif(bnd_idx == 3) then
                              g_Cg_intp   = (g_ice_Cg_tmp(3)-g_ice_Cg_tmp(2))/(1.655_r8-1.055_r8)*(1.1_r8-1.055_r8)&
                                             +g_ice_Cg_tmp(2)
                              gg_F07_intp = (gg_ice_F07_tmp(3)-gg_ice_F07_tmp(2))/(1.655_r8-1.055_r8)*(1.1_r8-1.055_r8)&
                                          +gg_ice_F07_tmp(2)
                           elseif(bnd_idx == 4) then
                              g_Cg_intp   = (g_ice_Cg_tmp(3)-g_ice_Cg_tmp(2))/(1.655_r8-1.055_r8)*(1.35_r8-1.055_r8)&
                                             +g_ice_Cg_tmp(2)
                              gg_F07_intp = (gg_ice_F07_tmp(3)-gg_ice_F07_tmp(2))/(1.655_r8-1.055_r8)*(1.35_r8-1.055_r8)&
                                          +gg_ice_F07_tmp(2)
                           elseif(bnd_idx == 5) then
                              g_Cg_intp   = (g_ice_Cg_tmp(6)-g_ice_Cg_tmp(5))/(3.75_r8-3.0_r8)*(3.25_r8-3.0_r8)&
                                             +g_ice_Cg_tmp(5)
                              gg_F07_intp = (gg_ice_F07_tmp(6)-gg_ice_F07_tmp(5))/(3.75_r8-3.0_r8)*(3.25_r8-3.0_r8)&
                                          +gg_ice_F07_tmp(5)
                           endif
                              g_ice_F07 = gg_F07_intp + (1._r8 - gg_F07_intp) / ss_alb_snw_lcl(i) / 2._r8
                              g_ice = g_ice_F07 * g_Cg_intp
                              asm_prm_snw_lcl(i) = g_ice
                        endif

                        if(asm_prm_snw_lcl(i) > 0.99_r8) then
                           asm_prm_snw_lcl(i) = 0.99_r8
                        endif
                     enddo 
                  endif   !snl_top < 0
                  ! write(*,*)'snl_top < 0'

         !H. Wang
                  ! aerosol species 1 optical properties
                  ! ss_alb_aer_lcl(1)        = ss_alb_bc1(bnd_idx)
                  ! asm_prm_aer_lcl(1)       = asm_prm_bc1(bnd_idx)
                  ! ext_cff_mss_aer_lcl(1)   = ext_cff_mss_bc1(bnd_idx)

                  ! aerosol species 2 optical properties
                  ! ss_alb_aer_lcl(2)        = ss_alb_bc2(bnd_idx)
                  ! asm_prm_aer_lcl(2)       = asm_prm_bc2(bnd_idx)
                  ! ext_cff_mss_aer_lcl(2)   = ext_cff_mss_bc2(bnd_idx)
         !H. Wang
                  ! aerosol species 3 optical properties
                  ss_alb_aer_lcl(3)        = ss_alb_oc1(bnd_idx)
                  asm_prm_aer_lcl(3)       = asm_prm_oc1(bnd_idx)
                  ext_cff_mss_aer_lcl(3)   = ext_cff_mss_oc1(bnd_idx)

                  ! aerosol species 4 optical properties
                  ss_alb_aer_lcl(4)        = ss_alb_oc2(bnd_idx)
                  asm_prm_aer_lcl(4)       = asm_prm_oc2(bnd_idx)
                  ext_cff_mss_aer_lcl(4)   = ext_cff_mss_oc2(bnd_idx)

                  ! aerosol species 5 optical properties
                  ss_alb_aer_lcl(5)        = ss_alb_dst1(bnd_idx)
                  asm_prm_aer_lcl(5)       = asm_prm_dst1(bnd_idx)
                  ext_cff_mss_aer_lcl(5)   = ext_cff_mss_dst1(bnd_idx)

                  ! aerosol species 6 optical properties
                  ss_alb_aer_lcl(6)        = ss_alb_dst2(bnd_idx)
                  asm_prm_aer_lcl(6)       = asm_prm_dst2(bnd_idx)
                  ext_cff_mss_aer_lcl(6)   = ext_cff_mss_dst2(bnd_idx)

                  ! aerosol species 7 optical properties
                  ss_alb_aer_lcl(7)        = ss_alb_dst3(bnd_idx)
                  asm_prm_aer_lcl(7)       = asm_prm_dst3(bnd_idx)
                  ext_cff_mss_aer_lcl(7)   = ext_cff_mss_dst3(bnd_idx)

                  ! aerosol species 8 optical properties
                  ss_alb_aer_lcl(8)        = ss_alb_dst4(bnd_idx)
                  asm_prm_aer_lcl(8)       = asm_prm_dst4(bnd_idx)
                  ext_cff_mss_aer_lcl(8)   = ext_cff_mss_dst4(bnd_idx)  
                  ! write(*,*)'aer'
                  ! 1. snow and aerosol layer column mass (L_snw, L_aer [kg/m^2])
                  ! 2. optical Depths (tau_snw, tau_aer)
                  ! 3. weighted Mie properties (tau, omega, g)

                  ! Weighted Mie parameters of each layer
                  ! write(*,*)'snl_top',snl_top
                  do i=snl_top,nbr_lyr,1
#ifdef MODAL_AER      
                     !mgf++ within-ice and external BC optical properties
                     !
                     ! Lookup table indices for BC optical properties,
                     ! dependent on snow grain size and BC particle
                     ! size.

                     ! valid for 25 < snw_rds < 1625 um:
                     if (snw_rds_lcl(i) < 125) then
                        tmp1 = snw_rds_lcl(i)/50
                        idx_bcint_icerds = nint(tmp1)
                     elseif (snw_rds_lcl(i) < 175) then
                        idx_bcint_icerds = 2
                     else
                        tmp1 = (snw_rds_lcl(i)/250)+2
                        idx_bcint_icerds = nint(tmp1)
                     endif

                     ! valid for 25 < bc_rds < 525 nm
                     idx_bcint_nclrds = nint(rds_bcint_lcl(i)/50)
                     idx_bcext_nclrds = nint(rds_bcext_lcl(i)/50)

                     ! check bounds:
                     if (idx_bcint_icerds < idx_bcint_icerds_min) idx_bcint_icerds = idx_bcint_icerds_min
                     if (idx_bcint_icerds > idx_bcint_icerds_max) idx_bcint_icerds = idx_bcint_icerds_max
                     if (idx_bcint_nclrds < idx_bc_nclrds_min) idx_bcint_nclrds = idx_bc_nclrds_min
                     if (idx_bcint_nclrds > idx_bc_nclrds_max) idx_bcint_nclrds = idx_bc_nclrds_max
                     if (idx_bcext_nclrds < idx_bc_nclrds_min) idx_bcext_nclrds = idx_bc_nclrds_min
                     if (idx_bcext_nclrds > idx_bc_nclrds_max) idx_bcext_nclrds = idx_bc_nclrds_max

                     ! retrieve absorption enhancement factor for within-ice BC
                     enh_fct = bcenh(bnd_idx,idx_bcint_nclrds,idx_bcint_icerds)

                     ! get BC optical properties (moved from above)
                     ! aerosol species 1 optical properties (within-ice BC)
                     ss_alb_aer_lcl(1)        = ss_alb_bc1(bnd_idx,idx_bcint_nclrds)
                     asm_prm_aer_lcl(1)       = asm_prm_bc1(bnd_idx,idx_bcint_nclrds)
                     ext_cff_mss_aer_lcl(1)   = ext_cff_mss_bc1(bnd_idx,idx_bcint_nclrds)*enh_fct

                     ! aerosol species 2 optical properties (external BC)
                     ss_alb_aer_lcl(2)        = ss_alb_bc2(bnd_idx,idx_bcext_nclrds)
                     asm_prm_aer_lcl(2)       = asm_prm_bc2(bnd_idx,idx_bcext_nclrds)
                     ext_cff_mss_aer_lcl(2)   = ext_cff_mss_bc2(bnd_idx,idx_bcext_nclrds)

#else
                     ! bulk aerosol treatment (BC optical properties independent
                     ! of BC and ice grain size)
                     ! aerosol species 1 optical properties (within-ice BC)
                     ss_alb_aer_lcl(1)        = ss_alb_bc1(bnd_idx)
                     asm_prm_aer_lcl(1)       = asm_prm_bc1(bnd_idx)
                     ext_cff_mss_aer_lcl(1)   = ext_cff_mss_bc1(bnd_idx)

                     ! aerosol species 2 optical properties
                     ss_alb_aer_lcl(2)        = ss_alb_bc2(bnd_idx)
                     asm_prm_aer_lcl(2)       = asm_prm_bc2(bnd_idx)
                     ext_cff_mss_aer_lcl(2)   = ext_cff_mss_bc2(bnd_idx)
#endif

                     ! Calculate single-scattering albedo for internal mixing of dust-snow
                     if (use_dust_snow_internal_mixing) then
                           if (bnd_idx < 4) then
                              C_dust_total = mss_cnc_aer_lcl(i,5) + mss_cnc_aer_lcl(i,6) &
                                          + mss_cnc_aer_lcl(i,7) + mss_cnc_aer_lcl(i,8)
                              C_dust_total = C_dust_total * 1.0E+06_r8
                              if(C_dust_total > 0._r8) then
                                 if (flg_slr_in == 1) then
                                    R_1_omega_tmp = dust_clear_d0(bnd_idx) &
                                                + dust_clear_d2(bnd_idx)*(C_dust_total**dust_clear_d1(bnd_idx))
                                 else
                                    R_1_omega_tmp = dust_cloudy_d0(bnd_idx) &
                                                + dust_cloudy_d2(bnd_idx)*(C_dust_total**dust_cloudy_d1(bnd_idx))
                                 endif
                                 ss_alb_snw_lcl(i) = 1.0_r8 - (1.0_r8 - ss_alb_snw_lcl(i)) *R_1_omega_tmp
                              endif
                           endif
                           do j = 5,8,1
                              ss_alb_aer_lcl(j)        = 0._r8
                              asm_prm_aer_lcl(j)       = 0._r8
                              ext_cff_mss_aer_lcl(j)   = 0._r8
                           enddo
                     endif

                     !mgf--
                     if (i>0) then 
                        L_snw(i)   = rho_snw(i)*dziw(i)
                     else 
                        L_snw(i)   = wice_soisno(i) + wliq_soisno(i)
                        ! write(*,*)i,'wice_soisno',wice_soisno(i),'wice_soisno',wliq_soisno(i)
                     endif
                     ! write(*,*)i,'h2osno_ice_lcl',h2osno_ice_lcl(i)
                     ! write(*,*)i,'h2osno_liq_lcl',h2osno_liq_lcl(i)
                     tau_snw(i) = L_snw(i)*ext_cff_mss_snw_lcl(i)
                     ! write(*,*)i,'tau_snw',tau_snw(i)
                     !write(*,*)i,'ss_alb_snw_lcl',ss_alb_snw_lcl(i)
                     !write(*,*)i,'asm_prm_snw_lcl',asm_prm_snw_lcl(i)

                     do j=1,sno_nbr_aer
                        if (use_dust_snow_internal_mixing .and. (j >= 5)) then
                           L_aer(i,j)  = 0._r8
                        else
                           L_aer(i,j)   = L_snw(i)*mss_cnc_aer_lcl(i,j)
                        endif
                        if (L_aer(i,j)/= 0 ) then 
                           write(*,*)'mss_cnc_aer_lcl',j,mss_cnc_aer_lcl(i,j)
                        endif
                        ! write(*,*)'mss_cnc_aer_lcl',mss_cnc_aer_lcl(i,j)
                        tau_aer(i,j) = L_aer(i,j)*ext_cff_mss_aer_lcl(j)
                        ! write(*,*)j,'ext_cff_mss_aer_lcl',ext_cff_mss_aer_lcl(j)
                     enddo
                     ! if (i<=0) then 
                     !     L_snw(i)   = h2osno_ice_lcl(i)+h2osno_liq_lcl(i)
                     ! else if(i<=nl_ice) then 
                     !     L_snw(i)   = rho_snw(i)*dziw(i) - sum(L_aer(i,:))
                     ! else 
                     !     L_snw(i)   = ziw(i)
                     ! endif 

                     ! tau_snw(i) = L_snw(i)*ext_cff_mss_snw_lcl(i)

                     tau_sum   = 0._r8
                     omega_sum = 0._r8
                     g_sum     = 0._r8

                     do j=1,sno_nbr_aer
                        tau_sum    = tau_sum + tau_aer(i,j)
                        omega_sum  = omega_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j))
                        g_sum      = g_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j)*asm_prm_aer_lcl(j))
                     enddo
                     ! write(*,*)'tau_sum',tau_sum
                     ! write(*,*)'omega_sum',omega_sum
                     ! write(*,*)'g_sum',g_sum
                     ! if (i<=nl_ice) then 
                     tau(i)    = tau_sum + tau_snw(i)
                     omega(i)  = (1/tau(i))*(omega_sum+(ss_alb_snw_lcl(i)*tau_snw(i)))
                     g(i)      = (1/(tau(i)*omega(i)))*(g_sum+ (asm_prm_snw_lcl(i)*ss_alb_snw_lcl(i)*tau_snw(i)))
                     ! write(*,*)i,'tau',tau(i),'omega',omega(i),'g',g(i)
                     ! else 
                     !     tau(i)    = tau_sum
                     !     omega(i)  = (1/tau(i))*(omega_sum+(ss_alb_snw_lcl(i)*tau_snw(i)))
                     !     g(i)      = (1/(tau(i)*omega(i)))*(g_sum+ (asm_prm_snw_lcl(i)*ss_alb_snw_lcl(i)*tau_snw(i)))
                     ! endif 
                  enddo ! endWeighted Mie parameters of each layer  
                  ! write(*,*)'endWeighted Mie parameters of each layer'
                  ! DELTA transformations, if requested
                  if (DELTA == 1) then
                     do i=snl_top,nbr_lyr,1
                        g_star(i)     = g(i)/(1+g(i))
                        omega_star(i) = ((1-(g(i)**2))*omega(i)) / (1-(omega(i)*(g(i)**2)))
                        tau_star(i)   = (1-(omega(i)*(g(i)**2)))*tau(i)
                     enddo
                  else
                     do i=snl_top,nbr_lyr,1
                        g_star(i)     = g(i)
                        omega_star(i) = omega(i)
                        tau_star(i)   = tau(i)
                     enddo
                  endif
                  ! write(*,*)'DELTA transformations, if requested'

                  ! Begin radiative transfer solver
                  ! Given input vertical profiles of optical properties, evaluate the
                  ! monochromatic Delta-Eddington adding-doubling solution

                  ! note that trndir, trntdr, trndif, rupdir, rupdif, rdndif
                  ! are variables at the layer interface,
                  ! for snow with layers rangeing from snl_top to snl_btm
                  ! there are snl_top to snl_btm+1 layer interface
                  snl_btm_itf = nbr_lyr + 1

                  do i = snl_top,snl_btm_itf,1
                     trndir(i) = c0
                     trntdr(i) = c0
                     trndif(i) = c0
                     rupdir(i) = c0
                     rupdif(i) = c0
                     rdndif(i) = c0
                  enddo
                  ! write(*,*)'snl_btm_itf'
                  ! initialize top interface of top layer
                  trndir(snl_top) = c1
                  trntdr(snl_top) = c1
                  trndif(snl_top) = c1
                  rdndif(snl_top) = c0
                  ! write(*,*)'snl_top',snl_top
                  ! begin main level loop
                  ! for layer interfaces except for the very bottom
                  ! write(*,*)'nbr_lyr',nbr_lyr
                  do i = snl_top,nbr_lyr,1
                     ! write(*,*)'i',i

                  ! initialize all layer apparent optical properties to 0
                     rdir  (i) = c0
                     rdif_a(i) = c0
                     rdif_b(i) = c0
                     tdir  (i) = c0
                     tdif_a(i) = c0
                     tdif_b(i) = c0
                     trnlay(i) = c0
                     ! write(*,*)i,'check',trntdr(i),trmin

                  ! compute next layer Delta-eddington solution only if total transmission
                  ! of radiation to the interface just above the layer exceeds trmin.

                     if (trntdr(i) > trmin ) then
                        ! calculation over layers with penetrating radiation

                        ! delta-transformed single-scattering properties
                        ! of this layer
                        ts = tau_star(i)
                        ws = omega_star(i)
                        gs = g_star(i)
                        ! write(*,*)i,'ts',ts,'ws',ws,'gs',gs

                        ! Delta-Eddington solution expressions
                        ! write(*,*)'c3*(c1-ws)',c3*(c1-ws)
                        ! write(*,*)'c1 - ws*gs',c1 - ws*gs
                        ! write(*,*)'c3*(c1-ws)*(c1 - ws*gs)',c3*(c1-ws)*(c1 - ws*gs)
                        ! write(*,*)'lm',sqrt(c3*(c1-ws)*(c1 - ws*gs))
                        lm = sqrt(c3*(c1-ws)*(c1 - ws*gs))  !lm = el(ws,gs)
                        ! write(*,*)'lm',lm                        
                        ue = c1p5*(c1 - ws*gs)/lm           !ue = u(ws,gs,lm)
                        ! write(*,*)'ue',ue
                        extins = max(exp_min, exp(-lm*ts))
                        ! write(*,*)'extins',extins
                        ne = ((ue+c1)*(ue+c1)/extins) - ((ue-c1)*(ue-c1)*extins) !ne = n(ue,extins)
                        ! write(*,*)'lm ue extins ne',lm,ue,extins,ne
                        ! first calculation of rdif, tdif using Delta-Eddington formulas
                        ! rdif_a(k) = (ue+c1)*(ue-c1)*(c1/extins - extins)/ne
                        rdif_a(i) = (ue**2-c1)*(c1/extins - extins)/ne
                        tdif_a(i) = c4*ue/ne
                        ! evaluate rdir,tdir for direct beam
                        ! trnlay(i) = max(exp_min, exp(-ts/mu0n))
                        trnlay(i) = max(exp_min, exp(-ts/mu_not))

                        ! Delta-Eddington solution expressions
                        ! alpha(w,uu,gg,e) = p75*w*uu*((c1 + gg*(c1-w))/(c1 - e*e*uu*uu))
                        ! agamm(w,uu,gg,e) = p5*w*((c1 + c3*gg*(c1-w)*uu*uu)/(c1-e*e*uu*uu))
                        ! alp = alpha(ws,mu_not,gs,lm)
                        ! gam = agamm(ws,mu_not,gs,lm)
                        alp = cp75*ws*mu_not*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu_not*mu_not))
                        ! write(*,*)'cp75*ws*mu_not',cp75*ws*mu_not
                        ! write(*,*)'c1 + gs*(c1-ws)',c1 + gs*(c1-ws)
                        ! write(*,*)'c1 - lm*lm*mu_not*mu_not',c1 - lm*lm*mu_not*mu_not
                        gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu_not*mu_not)/(c1-lm*lm*mu_not*mu_not))
                        ! alp = cp75*ws*mu0n*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu0n*mu0n))
                        ! gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu0n*mu0n)/(c1-lm*lm*mu0n*mu0n))
                        apg = alp + gam
                        amg = alp - gam
                        ! write(*,*),'apg',apg,'amg',amg
                        rdir(i) = apg*rdif_a(i) +  amg*(tdif_a(i)*trnlay(i) - c1)
                        tdir(i) = apg*tdif_a(i) + (amg* rdif_a(i)-apg+c1)*trnlay(i)
                        ! write(*,*)i,'rdir',rdir(i),'tdir',tdir(i)

                        ! recalculate rdif,tdif using direct angular integration over rdir,tdir,
                        ! since Delta-Eddington rdif formula is not well-behaved (it is usually
                        ! biased low and can even be negative); use ngmax angles and gaussian
                        ! integration for most accuracy:
                        R1 = rdif_a(i) ! use R1 as temporary
                        T1 = tdif_a(i) ! use T1 as temporary
                        swt = c0
                        smr = c0
                        smt = c0
                        do ng=1,ngmax
                           mu  = difgauspt(ng)
                           gwt = difgauswt(ng)
                           swt = swt + mu*gwt
                           trn = max(exp_min, exp(-ts/mu))
                           ! alp = alpha(ws,mu,gs,lm)
                           ! gam = agamm(ws,mu,gs,lm)
                           alp = cp75*ws*mu*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu*mu))
                           gam = cp5*ws*((c1 + c3*gs*(c1-ws)*mu*mu)/(c1-lm*lm*mu*mu))
                           apg = alp + gam
                           amg = alp - gam
                           rdr = apg*R1 + amg*T1*trn - amg
                           tdr = apg*T1 + amg*R1*trn - apg*trn + trn
                           smr = smr + mu*rdr*gwt
                           smt = smt + mu*tdr*gwt
                        enddo      ! ng
                        rdif_a(i) = smr/swt
                        tdif_a(i) = smt/swt
                        ! write(*,*)i,'rdif_a',rdif_a(i),'tdif_a',tdif_a(i)
                        ! homogeneous layer
                        rdif_b(i) = rdif_a(i)
                        tdif_b(i) = tdif_a(i)  

                        if( i == kfrsnl1 .or. i == kfrsnl2 .or. i == kfrsnl3.or. i == kfrsnl4) then
                           ! compute fresnel reflection and transmission
                           ! amplitudes for two polarizations: 1=perpendicular
                           ! and 2=parallel to he plane containing incident,
                           ! reflected and refracted rays.
                           if ( i == kfrsnl1.or. i == kfrsnl4) then 
                              nreal = nreal 
                           elseif (i == kfrsnl2) then 
                              nreal = nreal_wtr
                           elseif (i == kfrsnl3) then 
                              nreal = nreal_wtr/nreal
                           endif
                           ! write(*,*)'nreal',nreal
                           !mu0n      = sqrt(c1-((c1-(coszen(c_idx))**2)/(refindx*refindx))) !SZA under the refractive boundary 
                           !mu0n      = sqrt(c1-((c1-(coszen(c_idx))**2)/nreal))
                           mu0n      = cos(asin(sin(acos(coszen))/nreal))
                           mu0= mu_not
                           if (isnan(nreal) .or. isnan(mu0n)) then
                              write(iulog,*) "CAW c",bnd_idx,"NAN"
                              write(iulog,*) "CAW c",bnd_idx,"mu0n=", mu0n
                              write(iulog,*) "CAW c",bnd_idx,"nreal=", nreal
                              write(iulog,*) "CAW c",bnd_idx,"temp2=", temp2
                              write(iulog,*) "CAW c",bnd_idx,"temp1=", temp1
                              !   call endrun(decomp_index=c_idx, elmlevel=namec, msg=errmsg(__FILE__, __LINE__))
                           endif


                           !write(iulog,*) "CAW c",c_idx,"bnd",bnd_idx,"mu0n", mu0n                                
                           !write(iulog,*) "CAW c",c_idx,"bnd",bnd_idx,"nreal", nreal
                           !write(iulog,*) "CAW c",c_idx,"bnd",bnd_idx,"dif_a(bnd_idx)", dif_a(bnd_idx)
                           !R1 = (mu0 - refindx*mu0n) / &
                           !     (mu0 + refindx*mu0n)
                           !R2 = (refindx*mu0 - mu0n) / &
                           !     (refindx*mu0 + mu0n)
                           !T1 = c2*mu0 / &
                           !     (mu0 + refindx*mu0n)
                           !T2 = c2*mu0 / &
                           !     (refindx*mu0 + mu0n)
                           R1 = (mu0 - nreal*mu0n) / &
                                    (mu0 + nreal*mu0n)
                           R2 = (nreal*mu0 - mu0n) / &
                                    (nreal*mu0 + mu0n)
                           T1 = c2*mu0 / &
                                    (mu0 + nreal*mu0n)
                           T2 = c2*mu0 / &
                                    (nreal*mu0 + mu0n)
                        ! write(iulog,*) "CAW c",c_idx,"bnd",bnd_idx,"mu0", mu0
                        ! write(iulog,*) "CAW c",c_idx,"bnd",bnd_idx,"mu_not", mu_not
                        ! write(iulog,*) "CAW c",c_idx,"bnd",bnd_idx,"mu0n", mu0n
                        ! write(iulog,*) "CAW c",c_idx,"bnd",bnd_idx,"R1",R1
                        ! write(iulog,*) "CAW c",c_idx,"bnd",bnd_idx,"R2",R2
                        ! write(iulog,*) "CAW c",c_idx,"bnd",bnd_idx,"T1",T1
                        ! write(iulog,*) "CAW c",c_idx,"bnd",bnd_idx,"T2",T2


                           ! unpolarized light for direct beam
                           Rf_dir_a = cp5 * (R1*R1 + R2*R2)
                           !Tf_dir_a = cp5 * (T1*T1 + T2*T2)*refindx*mu0n/mu0
                           Tf_dir_a = cp5 * (T1*T1 + T2*T2)*nreal*mu0n/mu0
                           
                           !write(iulog,*) "CAW c",c_idx,"bnd",bnd_idx,"Rf_dir_a",Rf_dir_a
                           !write(iulog,*) "CAW c",c_idx,"bnd",bnd_idx,"Tf_dir_a",Tf_dir_a
                           ! precalculated diffuse reflectivities and
                           ! transmissivities for incident radiation above and
                           ! below fresnel layer, using the direct albedos and
                           ! accounting for complete internal reflection from
                           ! below; precalculated because high order number of
                           ! gaussian points (~256) is required for convergence:
                           
                           ! above
                           !Rf_dif_a = cp063
                           Rf_dif_a = dif_a(bnd_idx)
                           Tf_dif_a = c1 - Rf_dif_a
                           ! below
                           !Rf_dif_b = cp455
                           Rf_dif_b = dif_b(bnd_idx)
                           Tf_dif_b = c1 - Rf_dif_b


                           ! the i = kfrsnl layer properties are updated to
                           ! combined the fresnel (refractive) layer, always
                           ! taken to be above he present layer i (i.e. be the
                           ! top interface):
                           rintfc   = c1 / (c1-Rf_dif_b*rdif_a(i))
                           ! write(*,*)i,'kfrsnl_rintfc',rintfc,'bnd_idx',bnd_idx,'flg_slr_in',flg_slr_in
                           tdir(i)   = Tf_dir_a*tdir(i) + &
                              Tf_dir_a*rdir(i) * &
                              Rf_dif_b*rintfc*tdif_a(i)
                           rdir(i)   = Rf_dir_a + &
                              Tf_dir_a*rdir(i) * &
                              rintfc*Tf_dif_b
                           rdif_a(i) = Rf_dif_a + &
                              Tf_dif_a*rdif_a(i) * &
                              rintfc*Tf_dif_b
                           rdif_b(i) = rdif_b(i) + &
                              tdif_b(i)*Rf_dif_b * &
                              rintfc*tdif_a(i)
                           tdif_a(i) = tdif_a(i)*rintfc*Tf_dif_a
                           tdif_b(i) = tdif_b(i)*rintfc*Tf_dif_b

                           ! update trnlay to include fresnel transmission
                           trnlay(i) = Tf_dir_a*trnlay(i)
                           !write(iulog,*) "CAW c",c_idx,"bnd",bnd_idx,"layer",i,"tdif_b",tdif_b(i)  
                           !write(iulog,*) "CAW c",c_idx,"bnd",bnd_idx,"layer",i,"tdif_a",tdif_a(i)                      


                        endif      ! i = kfrsnl

                     endif ! trntdr(k) > trmin
                     ! write(*,*)'trntdr(k) > trmin'

                     ! Calculate the solar beam transmission, total transmission, and
                     ! reflectivity for diffuse radiation from below at interface i,
                     ! the top of the current layer k:
                     !
                     !              layers       interface
                     !
                     !       ---------------------  i-1
                     !                i-1
                     !       ---------------------  i
                     !                 i
                     !       ---------------------

                     trndir(i+1) = trndir(i)*trnlay(i)
                     refkm1      = c1/(c1 - rdndif(i)*rdif_a(i))
                     tdrrdir     = trndir(i)*rdir(i)
                     tdndif      = trntdr(i) - trndir(i)
                     trntdr(i+1) = trndir(i)*tdir(i) + &
                              (tdndif + tdrrdir*rdndif(i))*refkm1*tdif_a(i)
                     rdndif(i+1) = rdif_b(i) + &
                              (tdif_b(i)*rdndif(i)*refkm1*tdif_a(i))
                     trndif(i+1) = trndif(i)*refkm1*tdif_a(i)
                     ! write(*,*)i,'refkm1',refkm1,'rdif_b',rdif_b(i)

                  enddo       ! i    end main level loop
                  ! write(*,*)'end main level loop'
                  ! compute reflectivity to direct and diffuse radiation for layers
                  ! below by adding succesive layers starting from the underlying
                  ! ground and working upwards:
                  !
                  !              layers       interface
                  !
                  !       ---------------------  i
                  !                 i
                  !       ---------------------  i+1
                  !                i+1
                  !       ---------------------

                  ! set the underlying ground albedo == albedo of near-IR
                  ! unless bnd_idx == 1, for visible
                  rupdir(snl_btm_itf) = albsfc(2)
                  rupdif(snl_btm_itf) = albsfc(2)
                  if (bnd_idx == 1) then
                     rupdir(snl_btm_itf) = albsfc(1)
                     rupdif(snl_btm_itf) = albsfc(1)
                  endif

                  do i=nbr_lyr,snl_top,-1
                     ! interface scattering
                     refkp1        = c1/( c1 - rdif_b(i)*rupdif(i+1))
                     ! dir from top layer plus exp tran ref from lower layer, interface
                     ! scattered and tran thru top layer from below, plus diff tran ref
                     ! from lower layer with interface scattering tran thru top from below
                     rupdir(i) = rdir(i) &
                              + (        trnlay(i)  *rupdir(i+1) &
                              +  (tdir(i)-trnlay(i))*rupdif(i+1))*refkp1*tdif_b(i)
                     ! dif from top layer from above, plus dif tran upwards reflected and
                     ! interface scattered which tran top from below
                     rupdif(i) = rdif_a(i) + tdif_a(i)*rupdif(i+1)*refkp1*tdif_b(i)
                     ! write(*,*)i,'rupdir',rupdir(i),'rupdif',rupdif(i), &
                     !          'rdif_a',rdif_a(i),'tdif_a',tdif_a(i),'tdif_b',tdif_b(i),'refkp1',refkp1
                     
                  enddo       ! i
                  ! write(*,*)'rupdif'
                  ! net flux (down-up) at each layer interface from the
                  ! snow top (i = snl_top) to bottom interface above land (i = snl_btm_itf)
                  ! the interface reflectivities and transmissivities required
                  ! to evaluate interface fluxes are returned from solution_dEdd;
                  ! now compute up and down fluxes for each interface, using the
                  ! combined layer properties at each interface:
                  !
                  !              layers       interface
                  !
                  !       ---------------------  i
                  !                 i
                  !       ---------------------

                  do i = snl_top, snl_btm_itf
                     ! interface scattering
                     refk          = c1/(c1 - rdndif(i)*rupdif(i))
                     ! dir tran ref from below times interface scattering, plus diff
                     ! tran and ref from below times interface scattering
                     fdirup(i) = (trndir(i)*rupdir(i) + &
                                       (trntdr(i)-trndir(i))  &
                                       *rupdif(i))*refk
                     ! dir tran plus total diff trans times interface scattering plus
                     ! dir tran with up dir ref and down dif ref times interface scattering
                     fdirdn(i) = trndir(i) + (trntdr(i) &
                                    - trndir(i) + trndir(i)  &
                                    *rupdir(i)*rdndif(i))*refk
                     ! diffuse tran ref from below times interface scattering
                     fdifup(i) = trndif(i)*rupdif(i)*refk
                     ! diffuse tran times interface scattering
                     fdifdn(i) = trndif(i)*refk
                     ! write(*,*)i,'fdirup',fdirup(i),'fdirdn',fdirdn(i),'fdifup',fdifup(i),'fdifdn',fdifdn(i)
                     ! write(*,*)i
                     ! write(*,*)i
                     ! write(*,*)i
                     
                     ! netflux, down - up
                     ! dfdir = fdirdn - fdirup
                     dfdir(i) = trndir(i) &
                                 + (trntdr(i)-trndir(i)) * (c1 - rupdif(i)) * refk &
                                 -  trndir(i)*rupdir(i)  * (c1 - rdndif(i)) * refk
                     if (dfdir(i) < puny) dfdir(i) = c0
                     ! write(*,*)i,'dfdir',dfdir(i)
                     ! dfdif = fdifdn - fdifup
                     dfdif(i) = trndif(i) * (c1 - rupdif(i)) * refk
                     if (dfdif(i) < puny) dfdif(i) = c0
                     ! write(*,*)i,'dfdif',dfdif(i)
                  enddo       ! i
                  ! write(*,*)'fdifup'
                  if (flg_slr_in == 1) then 
                     albedo = rupdir(snl_top)
                     dftmp = dfdir
                  else 
                     albedo = rupdif(snl_top)
                     dftmp = dfdif 
                  endif 
                  albout_lcl(bnd_idx) = albedo
               
                  do i = snl_top, nl_lake 
                     if (nbr_lyr == 10 .or. i < kfrsnl3-1) then 
                        F_abs(i) = dftmp(i)-dftmp(i+1)
                     elseif (i == kfrsnl3-1) then
                        F_abs(i) = dftmp(i)-dftmp(i+2)
                     elseif (i > kfrsnl3-1) then
                        F_abs(i) = dftmp(i+1)-dftmp(i+2)
                     endif 
                     flx_abs_lcl(i,bnd_idx) = F_abs(i)
                     ! write(*,*)i,'flx_abs_lcl',flx_abs_lcl(i,bnd_idx)
                     ! flx_abs_btm(bnd_idx) = dftmp(snl_btm_itf)
                  enddo

                  if (flg_slr_in == 1) then 
                     flx_abs_up_lcl(bnd_idx) =(trndir(snl_top)*rupdir(snl_top) + (trntdr(snl_top)-trndir(snl_top))  &
                                          *rupdif(snl_top))*refk
                  else
                     flx_abs_up_lcl(bnd_idx) = trndif(snl_top) *  rupdif(snl_top)* refk 
                  endif
                  if (flg_slr_in == 1) then 
                     flx_abs_dn_lcl(bnd_idx) =trndir(snl_btm_itf) + (trntdr(snl_btm_itf) - trndir(snl_btm_itf) + trndir(snl_btm_itf)  &
                                       *rupdir(snl_btm_itf)*rdndif(snl_btm_itf))*refk
                     ! flx_abs_dn_lcl(bnd_idx) = dftmp(snl_btm_itf)
                  else
                     flx_abs_dn_lcl(bnd_idx) = trndif(snl_btm_itf)*refk 
                     ! flx_abs_dn_lcl(bnd_idx) = dftmp(snl_sbtm_itf)
                  endif 

                  ! albout_lcl(bnd_idx) = (fdirup(snl_top)+fdifup(snl_top))/(fdirdn(snl_top)+fdifdn(snl_top))
                  
               enddo   !bnd_idx
               ! Weight output NIR albedo appropriately
               albout(1) = albout_lcl(1)
               flx_sum         = 0._r8
               do bnd_idx= nir_bnd_bgn,nir_bnd_end
                  flx_sum = flx_sum + flx_wgt(bnd_idx)*albout_lcl(bnd_idx)
               enddo
               albout(2) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
               ! write(*,*)'albout',albout(:)
               flx_abs(:,1) = flx_abs_lcl(:,1)
               flx_abs_up(1) = flx_abs_up_lcl(1)
               flx_abs_dn(1) = flx_abs_dn_lcl(1)
               flx_sum = 0._r8
               do bnd_idx= nir_bnd_bgn,nir_bnd_end
                  flx_sum = flx_sum + flx_wgt(bnd_idx)*flx_abs_up_lcl(bnd_idx)
               enddo
               flx_abs_up(2) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
               ! write(*,*)'flx_abs_up',flx_abs_up(:),flg_slr_in
               flx_sum = 0._r8
               do bnd_idx= nir_bnd_bgn,nir_bnd_end
                  flx_sum = flx_sum + flx_wgt(bnd_idx)*flx_abs_dn_lcl(bnd_idx)
               enddo
               flx_abs_dn(2) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
               do i=snl_top,nl_lake,1
                  flx_sum = 0._r8
                  do bnd_idx= nir_bnd_bgn,nir_bnd_end
                     flx_sum = flx_sum + flx_wgt(bnd_idx)*flx_abs_lcl(i,bnd_idx)
                  enddo
                  flx_abs(i,2) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
                  !  write(*,*)'sumflx_wgt',sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))
               enddo
               if (flg_slr_in == 1) then 
                  flx_absd(:,:) = flx_abs(:,:)
                  flx_abs_upd(:) = flx_abs_up(:)
                  flx_abs_dnd(:) = flx_abs_dn(:)
                  ! write(*,*)'direct_vis',flx_absd(:,1),flx_abs_upd(1),flx_abs_dnd(1)
                  ! write(*,*)'direct_nir',flx_absd(:,2),flx_abs_upd(2),flx_abs_dnd(2)
               elseif (flg_slr_in == 2) then 
                  flx_absi(:,:) = flx_abs(:,:)
                  flx_abs_upi(:) = flx_abs_up(:)
                  flx_abs_dni(:) = flx_abs_dn(:)
                  ! write(*,*)'diffuse_vis',flx_absi(:,1),flx_abs_upi(1),flx_abs_dni(1)
                  ! write(*,*)'diffuse_nir',flx_absi(:,2),flx_abs_upi(2),flx_abs_dni(2)
               endif

            enddo ! flg_slr_in
            ! write(*,*)'direct_vis',sum(flx_absd(:,1))+flx_abs_dnd(1),1.*(1.-albout(1))
            ! write(*,*)'direct_nir',sum(flx_absd(:,2))+flx_abs_dnd(2),1.*(1.-albout(1))
            ! write(*,*)'diffuse_vis',sum(flx_absi(:,1))+flx_abs_dni(1),1.*(1.-albout(2))
            ! write(*,*)'diffuse_nir',sum(flx_absi(:,2))+flx_abs_dni(2),1.*(1.-albout(2))

            phi_up = flx_abs_upd(1)*forc_sols+flx_abs_upd(2)*forc_soll &
                     +flx_abs_upi(1)*forc_solsd+flx_abs_upi(2)*forc_solld
            phi_srf = 0 
            if (snl<0.or.lake_icefrac(1)== 1 ) then 
               phi_srf = flx_absd(snl_top,1)*forc_sols+flx_absd(snl_top,2)*forc_soll&
                        +flx_absi(snl_top,1)*forc_solsd+flx_absi(snl_top,2)*forc_solld
            elseif (srf_lyr == 1) then 
               phi_srf = flx_absd(srf_lyr,1)*forc_sols+flx_absd(srf_lyr,2)*forc_soll&
                           +flx_absi(srf_lyr,1)*forc_solsd+flx_absi(srf_lyr,2)*forc_solld
            elseif(srf_lyr > 1) then 
               do i = 1, srf_lyr
                  phi_srf = phi_srf+flx_absd(i,1)*forc_sols+flx_absd(i,2)*forc_soll&
                              +flx_absi(i,1)*forc_solsd+flx_absi(i,2)*forc_solld
               enddo
            endif
            ! write(*,*)'srf',srf_lyr,'phi_srf',phi_srf
            if(forc_sols+forc_soll+forc_solsd+forc_solld==0) then 
               srf_prop = 1
            else
               srf_prop = phi_srf/(forc_sols+forc_soll+forc_solsd+forc_solld-phi_up)
            endif
            ! write(*,*)'srf_prop',srf_prop
            if(snl<0.or.lake_icefrac(1)==1)then 
               phi(snl_top) = 0
               do i = snl_top+1, nl_lake
                  phi(i) = flx_absd(i,1)*forc_sols+flx_absd(i,2)*forc_soll&
                           +flx_absi(i,1)*forc_solsd+flx_absi(i,2)*forc_solld
               enddo
            else 
               do i = snl_top, nl_lake
                  if(i>0.and.i<=srf_lyr) then 
                     ! write(*,*)'i>0.and.i<=srf_lyr',i
                     phi(i)= 0 
                  else 
                     phi(i) = flx_absd(i,1)*forc_sols+flx_absd(i,2)*forc_soll&
                           +flx_absi(i,1)*forc_solsd+flx_absi(i,2)*forc_solld
                  endif 
               enddo 
            endif
            ! do i = snl_top,nl_lake
            !    write(*,*) i, phi(i)
            ! enddo
            phi_soil = (flx_abs_dnd(1)*forc_sols+flx_abs_dnd(2)*forc_soll &
                        +flx_abs_dni(1)*forc_solsd+flx_abs_dni(2)*forc_solld)
            ! write(*,*)'in',forc_sols+forc_soll+forc_solsd+forc_solld,'sum',phi_up+phi_soil+sum(phi(:))+phi_srf

   ! ======================================================================
   !*[3] Begin stability iteration and temperature and fluxes calculation
   ! ======================================================================


      ! =====================================
            ITERATION1 : DO WHILE (iter <= itmax)
            ! =====================================
      
               t_grnd_bef = t_grnd
      
               if (t_grnd_bef > tfrz .and. t_lake(1) > tfrz .and. snl == 0) then
                  tksur = savedtke1       !water molecular conductivity
                  tsur = t_lake(1)
                  htvp = hvap
               else if (snl == 0) then !frozen but no snow layers
                  tksur = tkice        ! This is an approximation because the whole layer may not be frozen, and it is not
                                       ! accounting for the physical (but not nominal) expansion of the frozen layer.
                  tsur = t_lake(1)
                  htvp = hsub
               else
                  ! need to calculate thermal conductivity of the top snow layer
                  rhosnow = (wice_soisno(lb)+wliq_soisno(lb))/dz_soisno(lb)
                  tksur = tkair + (7.75e-5*rhosnow + 1.105e-6*rhosnow*rhosnow)*(tkice-tkair)
                  tsur = t_soisno(lb)
                  htvp = hsub
               end if
      
      ! Evaluated stability-dependent variables using moz from prior iteration
               displax = 0.
               if (DEF_USE_CBL_HEIGHT) then
                  call moninobuk_leddy(forc_hgt_u,forc_hgt_t,forc_hgt_q,displax,z0mg,z0hg,z0qg,obu,um, hpbl, &
                              ustar,fh2m,fq2m,fm10m,fm,fh,fq)
               else
                  call moninobuk(forc_hgt_u,forc_hgt_t,forc_hgt_q,displax,z0mg,z0hg,z0qg,obu,um,&
                              ustar,fh2m,fq2m,fm10m,fm,fh,fq)
               endif
      
      ! Get derivative of fluxes with repect to ground temperature
               ram    = 1./(ustar*ustar/um)
               rah    = 1./(vonkar/fh*ustar)
               raw    = 1./(vonkar/fq*ustar)
               stftg3 = emg*stefnc*t_grnd_bef*t_grnd_bef*t_grnd_bef
      
               ax  = srf_prop*(forc_sols+forc_soll+forc_solsd+forc_solld-phi_up) + emg*forc_frl + 3.*stftg3*t_grnd_bef &
                     + forc_rhoair*cpair/rah*thm &
                     - htvp*forc_rhoair/raw*(qsatg-qsatgdT*t_grnd_bef - forc_q) &
                     + tksur*tsur/dzsur
      
               bx  = 4.*stftg3 + forc_rhoair*cpair/rah &
                     + htvp*forc_rhoair/raw*qsatgdT + tksur/dzsur
      
               t_grnd = ax/bx
      
               !-----------------------------------------------------------------
               ! h_fin = betaprime*sabg + emg*forc_frl + 3.*stftg3*t_grnd_bef & !
               !     + forc_rhoair*cpair/rah*thm &                              !
               !     - htvp*forc_rhoair/raw*(qsatg-qsatgdT*t_grnd_bef - forc_q) !
               ! h_finDT = 4.*stftg3 + forc_rhoair*cpair/rah &                  !
               !     + htvp*forc_rhoair/raw*qsatgdT                             !
               ! del_T_grnd = t_grnd - t_grnd_bef                               !
               !----------------------------------------------------------------!
      
      ! surface fluxes of momentum, sensible and latent
      ! using ground temperatures from previous time step
      
               fseng = forc_rhoair*cpair*(t_grnd-thm)/rah
               fevpg = forc_rhoair*(qsatg+qsatgdT*(t_grnd-t_grnd_bef)-forc_q)/raw
      
               call qsadv(t_grnd,forc_psrf,eg,degdT,qsatg,qsatgdT)
               dth = thm-t_grnd
               dqh = forc_q-qsatg
               tstar = vonkar/fh*dth
               qstar = vonkar/fq*dqh
               thvstar = tstar*(1.+0.61*forc_q)+0.61*th*qstar
               zeta = zldis*vonkar*grav*thvstar/(ustar**2*thv)
               if(zeta >= 0.) then     !stable
                  zeta = min(2.,max(zeta,1.e-6))
               else                    !unstable
                  zeta = max(-100.,min(zeta,-1.e-6))
               endif
               obu = zldis/zeta
               if(zeta >= 0.)then
                  um = max(ur,0.1)
               else
                  if (DEF_USE_CBL_HEIGHT) then !//TODO: Shaofeng, 2023.05.18
                     zii = max(5.*forc_hgt_u,hpbl)
                  endif !//TODO: Shaofeng, 2023.05.18
                  wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
                  wc2 = beta1*beta1*(wc*wc)
                  um = sqrt(ur*ur+wc2)
               endif
      
               call roughness_lake (snl,t_grnd,t_lake(1),lake_icefrac(1),forc_psrf,&
                                    cur,ustar,z0mg,z0hg,z0qg)
      
               iter = iter + 1
               del_T_grnd = abs(t_grnd - t_grnd_bef)
      
               if(iter .gt. itmin) then
                  if(del_T_grnd <= dtmin) then
                     convernum = convernum + 1
                  end if
                  if(convernum >= 4) EXIT
               endif
         
            ! ===============================================
            END DO ITERATION1   ! end of stability iteration
            ! ===============================================
      
      !*----------------------------------------------------------------------
      !*Zack Subin, 3/27/09
      !*Since they are now a function of whatever t_grnd was before cooling
      !*to freezing temperature, then this value should be used in the derivative correction term.
      !*Allow convection if ground temp is colder than lake but warmer than 4C, or warmer than
      !*lake which is warmer than freezing but less than 4C.
            tdmax = tfrz + 4.0
            if ( (snl < 0 .or. t_lake(1) <= tfrz) .and. t_grnd > tfrz) then
               t_grnd_bef = t_grnd
               t_grnd = tfrz
               fseng = forc_rhoair*cpair*(t_grnd-thm)/rah
               fevpg = forc_rhoair*(qsatg+qsatgdT*(t_grnd-t_grnd_bef)-forc_q)/raw
            else if ( (t_lake(1) > t_grnd .and. t_grnd > tdmax) .or. &
                     (t_lake(1) < t_grnd .and. t_lake(1) > tfrz .and. t_grnd < tdmax) ) then
                     ! Convective mixing will occur at surface
               t_grnd_bef = t_grnd
               t_grnd = t_lake(1)
               fseng = forc_rhoair*cpair*(t_grnd-thm)/rah
               fevpg = forc_rhoair*(qsatg+qsatgdT*(t_grnd-t_grnd_bef)-forc_q)/raw
            end if
            
      !*----------------------------------------------------------------------
      
      ! net longwave from ground to atmosphere
            stftg3 = emg*stefnc*t_grnd_bef*t_grnd_bef*t_grnd_bef
            olrg = (1.-emg)*forc_frl + emg*stefnc*t_grnd_bef**4 + 4.*stftg3*(t_grnd - t_grnd_bef)
            if (t_grnd > tfrz )then
            htvp = hvap
            else
            htvp = hsub
            end if
      
      !The actual heat flux from the ground interface into the lake, not including the light that penetrates the surface.
            fgrnd1 = srf_prop*(forc_sols+forc_soll+forc_solsd+forc_solld-phi_up) + forc_frl - olrg - fseng - htvp*fevpg
      
            ! January 12, 2023 by Yongjiu Dai
            IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
               hs = sabg_lyr(lb) + forc_frl - olrg - fseng - htvp*fevpg
               dhsdT = 0.0
            ENDIF
      
      !------------------------------------------------------------
      ! Set up vector r and vectors a, b, c that define tridiagonal matrix
      ! snow and lake and soil layer temperature
      !------------------------------------------------------------
      
      
      !------------------------------------------------------------
      ! Diffusivity and implied thermal "conductivity" = diffusivity * cwat
      !------------------------------------------------------------
      
            do j = 1, nl_lake
               cv_lake(j) = dz_lake(j) * (cwat*(1.-lake_icefrac(j)) + cice_eff*lake_icefrac(j))
            end do
      
            call hConductivity_lake(nl_lake,snl,t_grnd,&
                                    z_lake,t_lake,lake_icefrac,rhow,&
                                    dlat,ustar,z0mg,lakedepth,depthcrit,tk_lake,savedtke1)
      
      !------------------------------------------------------------
      ! Set the thermal properties of the snow above frozen lake and underlying soil
      ! and check initial energy content.
      !------------------------------------------------------------
      
            lb = snl+1
            do i = 1, nl_soil
               vf_water(i) = wliq_soisno(i)/(dz_soisno(i)*denh2o)
               vf_ice(i) = wice_soisno(i)/(dz_soisno(i)*denice)
               CALL soil_hcap_cond(vf_gravels(i),vf_om(i),vf_sand(i),porsl(i),&
                                    wf_gravels(i),wf_sand(i),k_solids(i),&
                                    csol(i),dkdry(i),dksatu(i),dksatf(i),&
                                    BA_alpha(i),BA_beta(i),&
                                    t_soisno(i),vf_water(i),vf_ice(i),hcap(i),thk(i))
               cv_soisno(i) = hcap(i)*dz_soisno(i)
            enddo
      
      ! Snow heat capacity and conductivity
            if(lb <=0 )then
               do j = lb, 0
                  cv_soisno(j) = cpliq*wliq_soisno(j) + cpice*wice_soisno(j)
                  rhosnow = (wice_soisno(j)+wliq_soisno(j))/dz_soisno(j)
                  thk(j) = tkair + (7.75e-5*rhosnow + 1.105e-6*rhosnow*rhosnow)*(tkice-tkair)
               enddo
            endif
      
      ! Thermal conductivity at the layer interface
            do i = lb, nl_soil-1
      
      ! the following consideration is try to avoid the snow conductivity
      ! to be dominant in the thermal conductivity of the interface.
      ! Because when the distance of bottom snow node to the interfacee
      ! is larger than that of interface to top soil node,
      ! the snow thermal conductivity will be dominant, and the result is that
      ! lees heat tranfer between snow and soil
      
      ! modified by Nan Wei, 08/25/2014
               if (i /= 0) then
                  tk_soisno(i) = thk(i)*thk(i+1)*(z_soisno(i+1)-z_soisno(i)) &
                        /(thk(i)*(z_soisno(i+1)-zi_soisno(i))+thk(i+1)*(zi_soisno(i)-z_soisno(i)))
               else
                  tk_soisno(i) = thk(i)
               end if
            end do
            tk_soisno(nl_soil) = 0.
            tktopsoil = thk(1)
      
      ! Sum cv_lake*t_lake for energy check
      ! Include latent heat term, and use tfrz as reference temperature
      ! to prevent abrupt change in heat content due to changing heat capacity with phase change.
      
            ! This will need to be over all soil / lake / snow layers. Lake is below.
            ocvts = 0.
            do j = 1, nl_lake
               ocvts = ocvts + cv_lake(j)*(t_lake(j)-tfrz) + cfus*dz_lake(j)*(1.-lake_icefrac(j))
            end do
      
            ! Now do for soil / snow layers
            do j = lb, nl_soil
               ocvts = ocvts + cv_soisno(j)*(t_soisno(j)-tfrz) + hfus*wliq_soisno(j)
               if (j == 1 .and. scv > 0. .and. j == lb) then
                  ocvts = ocvts - scv*hfus
               end if
            end do
      
         ENDIF                                              
         
      ELSE

         do j = 1, nl_lake
            ! extinction coefficient from surface data (1/m), if no eta from surface data,
            ! set eta, the extinction coefficient, according to L Hakanson, Aquatic Sciences, 1995
            ! (regression of secchi depth with lake depth for small glacial basin lakes), and the
            ! Poole & Atkins expression for extinction coeffient of 1.7 / secchi Depth (m).

            eta = 1.1925*max(lakedepth,1.)**(-0.424)
            zin  = z_lake(j) - 0.5*dz_lake(j)
            zout = z_lake(j) + 0.5*dz_lake(j)
            rsfin  = exp( -eta*max(  zin-za(idlak),0. ) )  ! the radiation within surface layer (z<za)
            rsfout = exp( -eta*max( zout-za(idlak),0. ) )  ! is considered fixed at (1-beta)*sabg
                                                            ! i.e, max(z-za, 0)
            ! Let rsfout for bottom layer go into soil.
            ! This looks like it should be robust even for pathological cases,
            ! like lakes thinner than za(idlak).

            phi(j) = (rsfin-rsfout) * sabg_lyr(1) * (1.-betaprime)
            if (j == nl_lake) phi_soil = rsfout * sabg_lyr(1) * (1.-betaprime)
         end do
      ENDIF
      
      phix(:) = 0.
      phix(snl+1:nl_lake)= phi(snl+1:nl_lake)         !lake layer
      phix(nl_lake+1) = phi_soil               !top soil layer
   !   write(*,*)'check phix',phix(:) 
      ! Set up interface depths(zx), and temperatures (tx).

      do j = lb, nl_lake+nl_soil
         jprime = j - nl_lake
         if (j <= 0) then                      !snow layer
            zx(j) = z_soisno(j)
            tx(j) = t_soisno(j)
         else if (j <= nl_lake) then           !lake layer
            zx(j) = z_lake(j)
            tx(j) = t_lake(j)
         else                                  !soil layer
            zx(j) = z_lake(nl_lake) + dz_lake(nl_lake)/2. + z_soisno(jprime)
            tx(j) = t_soisno(jprime)
         end if
      end do



      tx_bef = tx
      ! do j = lb, nl_lake+nl_soil
      !    write(*,*)j,'tx',tx(j),'tx_bef',tx_bef(j),'-',tx(j)-tx_bef(j)
      ! enddo
      

! Heat capacity and resistance of snow without snow layers (<1cm) is ignored during diffusion,
! but its capacity to absorb latent heat may be used during phase change.

      do j = lb, nl_lake+nl_soil
         jprime = j - nl_lake

         ! heat capacity [J/(m2 K)]
         if (j <= 0) then                      !snow layer
            cvx(j) = cv_soisno(j)
         else if (j <= nl_lake) then           !lake layer
            cvx(j) = cv_lake(j)
         else                                  !soil layer
            cvx(j) = cv_soisno(jprime)
         end if

! Determine interface thermal conductivities at layer interfaces [W/(m K)]

         if (j < 0) then                       !non-bottom snow layer
            tkix(j) = tk_soisno(j)
         else if (j == 0) then                 !bottom snow layer
            dzp = zx(j+1) - zx(j)
            tkix(j) = tk_lake(1)*tk_soisno(j)*dzp &
                    /(tk_soisno(j)*z_lake(1) + tk_lake(1)*(-zx(j)))
                               ! tk_soisno(0) is the conductivity at the middle of that layer
         else if (j < nl_lake) then            !non-bottom lake layer
            tkix(j) = (tk_lake(j)*tk_lake(j+1) * (dz_lake(j+1)+dz_lake(j))) &
                    / (tk_lake(j)*dz_lake(j+1) + tk_lake(j+1)*dz_lake(j))
         else if (j == nl_lake) then           !bottom lake layer
            dzp = zx(j+1) - zx(j)
            tkix(j) = (tktopsoil*tk_lake(j)*dzp &
                    / (tktopsoil*dz_lake(j)/2. + tk_lake(j)*z_soisno(1)))
         else !soil layer
            tkix(j) = tk_soisno(jprime)
         end if

      end do
      ! do j = lb, nl_lake+nl_soil
      !    write(*,*)'cvx',cvx(j)
      ! enddo
      
! Determine heat diffusion through the layer interface and factor used in computing
! tridiagonal matrix and set up vector r and vectors a, b, c that define tridiagonal
! matrix and solve system

      do j = lb, nl_lake+nl_soil
         factx(j) = deltim/cvx(j)
         if (j < nl_lake+nl_soil) then         !top or interior layer
            fnx(j) = tkix(j)*(tx(j+1)-tx(j))/(zx(j+1)-zx(j))
         else                                  !bottom soil layer
            fnx(j) = 0. !not used
         end if
      end do

      IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
         if (lb <= 0) then                        ! snow covered
            do j = lb, 1
               if (j == lb) then                  ! top snow layer
                  dzp  = zx(j+1)-zx(j)
                  a(j) = 0.0
                  b(j) = 1. + (1.-cnfac)*factx(j)*tkix(j)/dzp
                  c(j) = -(1.-cnfac)*factx(j)* tkix(j)/dzp
                  r(j) = tx_bef(j)+ factx(j)*(hs - dhsdT*tx_bef(j) + cnfac*fnx(j))
               else if (j <= 0) then              ! non-top snow layers
                  dzm  = (zx(j)-zx(j-1))
                  dzp  = (zx(j+1)-zx(j))
                  a(j) =   - (1.-cnfac)*factx(j)* tkix(j-1)/dzm
                  b(j) = 1.+ (1.-cnfac)*factx(j)*(tkix(j)/dzp + tkix(j-1)/dzm)
                  c(j) =   - (1.-cnfac)*factx(j)* tkix(j)/dzp
                  r(j) = tx_bef(j) + cnfac*factx(j)*(fnx(j) - fnx(j-1)) + factx(j)*sabg_lyr(j)
               else                               ! snow covered top lake layer
                  dzm  = (zx(j)-zx(j-1))
                  dzp  = (zx(j+1)-zx(j))
                  a(j) =   - (1.-cnfac)*factx(j)* tkix(j-1)/dzm
                  b(j) = 1.+ (1.-cnfac)*factx(j)*(tkix(j)/dzp + tkix(j-1)/dzm)
                  c(j) =   - (1.-cnfac)*factx(j)* tkix(j)/dzp
                  r(j) = tx_bef(j) + cnfac*factx(j)*(fnx(j) - fnx(j-1)) + factx(j)*(phix(j) + betaprime*sabg_lyr(j))
               endif
            enddo
         else
            j = 1                                 ! no snow covered top lake layer
            dzp  = zx(j+1)-zx(j)
            a(j) = 0.0
            b(j) = 1. + (1.-cnfac)*factx(j)*tkix(j)/dzp
            c(j) = -(1.-cnfac)*factx(j)* tkix(j)/dzp
            r(j) = tx_bef(j)+ factx(j)*(cnfac*fnx(j)+ phix(j)+fgrnd1)
         endif

         do j = 2, nl_lake+nl_soil
            if (j < nl_lake+nl_soil) then         ! middle lake and soil layers
               dzm  = (zx(j)-zx(j-1))
               dzp  = (zx(j+1)-zx(j))
               a(j) =   - (1.-cnfac)*factx(j)* tkix(j-1)/dzm
               b(j) = 1.+ (1.-cnfac)*factx(j)*(tkix(j)/dzp + tkix(j-1)/dzm)
               c(j) =   - (1.-cnfac)*factx(j)* tkix(j)/dzp
               r(j) = tx_bef(j) + cnfac*factx(j)*(fnx(j) - fnx(j-1)) + factx(j)*phix(j)
            else                                  ! bottom soil layer
               dzm  = (zx(j)-zx(j-1))
               a(j) =   - (1.-cnfac)*factx(j)*tkix(j-1)/dzm
               b(j) = 1.+ (1.-cnfac)*factx(j)*tkix(j-1)/dzm
               c(j) = 0.
               r(j) = tx_bef(j) - cnfac*factx(j)*fnx(j-1)
            end if
         end do
         
      ! January 12, 2023

      ELSE
         do j = lb, nl_lake+nl_soil
            if (j == lb) then                     ! top layer
                dzp  = zx(j+1)-zx(j)
                a(j) = 0.0
                b(j) = 1. + (1.-cnfac)*factx(j)*tkix(j)/dzp
                c(j) = -(1.-cnfac)*factx(j)* tkix(j)/dzp
                r(j) = tx_bef(j)+ factx(j)*(cnfac*fnx(j)+ phix(j)+fgrnd1)
            else if (j < nl_lake+nl_soil) then    ! middle layer
               dzm  = (zx(j)-zx(j-1))
               dzp  = (zx(j+1)-zx(j))
               a(j) =   - (1.-cnfac)*factx(j)* tkix(j-1)/dzm
               b(j) = 1.+ (1.-cnfac)*factx(j)*(tkix(j)/dzp + tkix(j-1)/dzm)
               c(j) =   - (1.-cnfac)*factx(j)* tkix(j)/dzp
               r(j) = tx_bef(j) + cnfac*factx(j)*(fnx(j) - fnx(j-1)) + factx(j)*phix(j)
            else                                  ! bottom soil layer
               dzm  = (zx(j)-zx(j-1))
               a(j) =   - (1.-cnfac)*factx(j)*tkix(j-1)/dzm
               b(j) = 1.+ (1.-cnfac)*factx(j)*tkix(j-1)/dzm
               c(j) = 0.
               r(j) = tx_bef(j) - cnfac*factx(j)*fnx(j-1)
            end if
         end do
         ! write(*,*)'fgrnd12',fgrnd1
      ENDIF
      
!------------------------------------------------------------
! Solve for tdsolution
!------------------------------------------------------------

         nl_sls = abs(snl) + nl_lake + nl_soil

         call tridia (nl_sls, a(lb:), b(lb:), c(lb:), r(lb:), tx(lb:))

      ! do j = lb, nl_lake+nl_soil
         ! write(*,*)j,'tx',tx(j),'tx_bef',tx_bef(j)
      ! enddo
      do j = lb, nl_lake + nl_soil
         jprime = j - nl_lake
         if (j < 1) then               ! snow layer
            t_soisno(j) = tx(j)
         else if (j <= nl_lake) then   ! lake layer
            t_lake(j) = tx(j)
         else                          ! soil layer
            t_soisno(jprime) = tx(j)
         end if
      end do

! additional variables for CWRF or other atmospheric models
      emis = emg
      z0m  = z0mg
      zol  = zeta
      rib  = min(5.,zol*ustar**2/(vonkar*vonkar/fh*um**2))

! radiative temperature
      trad = (olrg/stefnc)**0.25

! solar absorption below the surface.
      IF(DEF_USE_TWOSTREAM) THEN
         fgrnd = forc_sols+forc_soll+forc_solsd+forc_solld - phi_up + forc_frl - olrg - fseng - htvp*fevpg
      ELSE 
         fgrnd = sabg + forc_frl - olrg - fseng - htvp*fevpg
      ENDIF
      taux = -forc_rhoair*forc_us/ram
      tauy = -forc_rhoair*forc_vs/ram
      fsena = fseng
      fevpa = fevpg
      lfevpa = htvp*fevpg

! 2 m height air temperature and specific humidity
      tref = thm + vonkar/fh*dth * (fh2m/vonkar - fh/vonkar)
      qref = forc_q + vonkar/fq*dqh * (fq2m/vonkar - fq/vonkar)

! calculate sublimation, frosting, dewing
      qseva = 0.
      qsubl = 0.
      qsdew = 0.
      qfros = 0.

      if (fevpg >= 0.0) then
         if(lb < 0)then
            qseva = min(wliq_soisno(lb)/deltim, fevpg)
            qsubl = fevpg - qseva
         else
            qseva = min((1.-lake_icefrac(1))*1000.*dz_lake(1)/deltim, fevpg)
            qsubl = fevpg - qseva
         endif
      else
         if (t_grnd < tfrz) then
            qfros = abs(fevpg)
         else
            qsdew = abs(fevpg)
         end if
      end if

      ! write(*,*)'t_grnd',t_grnd
      ! do j = lb, nl_lake + nl_soil
         ! write(*,*)j,'tx',tx(j),'tx_bef',tx_bef(j)
      ! enddo


#if(defined CoLMDEBUG)
      ! sum energy content and total energy into lake for energy check. any errors will be from the
      !     tridiagonal solution.
      esum1 = 0.0
      esum2 = 0.0
      do j = lb, nl_lake + nl_soil
         esum1 = esum1 + (tx(j)-tx_bef(j))*cvx(j)
         esum2 = esum2 + (tx(j)-tfrz)*cvx(j)
         ! write(*,*)j,'esum',esum1,(tx(j)-tx_bef(j)),(tx(j)-tx_bef(j))*cvx(j)
      end do
      ! write(*,*)'esum1',esum1
                           ! fgrnd includes all the solar radiation absorbed in the lake,
      ! IF(DEF_USE_TWOSTREAM) THEN
      !    IF (.not. DEF_USE_SNICAR .or. present(urban_call)) THEN
      !       if (snl == 0) then
      !          errsoi = esum1/deltim - fgrnd + 2*pi*uavg*topc*rou*(exp(daT*(tao(nl_lake+1)-tao(1)))-exp(-daT*(tao(nl_lake+1)-tao(1))))/D
      !       else 
      !          errsoi = esum1/deltim - fgrnd
      !       end if 
      !    ELSE 
      !       errsoi = esum1/deltim - fgrnd + 2*pi*uavg*topc*rou*(exp(daT*(tao(nl_lake+1)-tao(1)))-exp(-daT*(tao(nl_lake+1)-tao(1))))/D
      !    END IF
      ! ELSE
         errsoi = esum1/deltim - fgrnd 
         ! print*,'errsoi',errsoi
      ! END IF 
      if(abs(errsoi) > 0.1) then
         write(6,*)'energy conservation error in LAND WATER COLUMN during tridiagonal solution,', &
                   'error (W/m^2):', errsoi, fgrnd
      end if
#endif
      

!------------------------------------------------------------
!*[4] Phase change
!------------------------------------------------------------

      sm = 0.0
      xmf = 0.0
      imelt_soisno(:) = 0
      imelt_lake(:) = 0

      IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
         wice_soisno_bef(lb:0) = wice_soisno(lb:0)
      ENDIF

      ! Check for case of snow without snow layers and top lake layer temp above freezing.

      if (snl == 0 .and. scv > 0. .and. t_lake(1) > tfrz) then
         heatavail = (t_lake(1) - tfrz) * cv_lake(1)
         melt = min(scv, heatavail/hfus)
         heatrem = max(heatavail - melt*hfus, 0.) !catch small negative value to keep t at tfrz
         t_lake(1) = tfrz + heatrem/(cv_lake(1))

         snowdp = max(0., snowdp*(1. - melt/scv))
         scv = scv - melt

         if (scv < 1.e-12) scv = 0.        ! prevent tiny residuals
         if (snowdp < 1.e-12) snowdp = 0.  ! prevent tiny residuals
         sm = sm + melt/deltim
         xmf = xmf + melt*hfus
      end if

      ! Lake phase change
      do j = 1,nl_lake
         if (t_lake(j) > tfrz .and. lake_icefrac(j) > 0.) then ! melting
            imelt_lake(j) = 1
            heatavail = (t_lake(j) - tfrz) * cv_lake(j)
            melt = min(lake_icefrac(j)*denh2o*dz_lake(j), heatavail/hfus)
                       !denh2o is used because layer thickness is not adjusted for freezing
            heatrem = max(heatavail - melt*hfus, 0.)  !catch small negative value to keep t at tfrz
         else if (t_lake(j) < tfrz .and. lake_icefrac(j) < 1.) then !freezing
            imelt_lake(j) = 2
            heatavail = (t_lake(j) - tfrz) * cv_lake(j)
            melt = max(-(1.-lake_icefrac(j))*denh2o*dz_lake(j), heatavail/hfus)
                       !denh2o is used because layer thickness is not adjusted for freezing
            heatrem = min(heatavail - melt*hfus, 0.)  !catch small positive value to keep t at tfrz
         end if
         ! Update temperature and ice fraction.
         if (imelt_lake(j) > 0) then
            lake_icefrac(j) = lake_icefrac(j) - melt/(denh2o*dz_lake(j))
            if (lake_icefrac(j) > 1.-1.e-12) lake_icefrac(j) = 1.  ! prevent tiny residuals
            if (lake_icefrac(j) < 1.e-12)    lake_icefrac(j) = 0.  ! prevent tiny residuals
            cv_lake(j) = cv_lake(j) + melt*(cpliq-cpice)           ! update heat capacity
            t_lake(j) = tfrz + heatrem/cv_lake(j)
            xmf = xmf + melt*hfus
         end if
      end do

      ! snow & soil phase change. currently, does not do freezing point depression.
      do j = snl+1,nl_soil
         if (t_soisno(j) > tfrz .and. wice_soisno(j) > 0.) then ! melting
            imelt_soisno(j) = 1
            heatavail = (t_soisno(j) - tfrz) * cv_soisno(j)
            melt = min(wice_soisno(j), heatavail/hfus)
            heatrem = max(heatavail - melt*hfus, 0.) !catch small negative value to keep t at tfrz
            if (j <= 0) sm = sm + melt/deltim
         else if (t_soisno(j) < tfrz .and. wliq_soisno(j) > 0.) then !freezing
            imelt_soisno(j) = 2
            heatavail = (t_soisno(j) - tfrz) * cv_soisno(j)
            melt = max(-wliq_soisno(j), heatavail/hfus)
            heatrem = min(heatavail - melt*hfus, 0.) !catch small positive value to keep t at tfrz
         end if

         ! Update temperature and soil components.
         if (imelt_soisno(j) > 0) then
            wice_soisno(j) = wice_soisno(j) - melt
            wliq_soisno(j) = wliq_soisno(j) + melt
            if (wice_soisno(j) < 1.e-12) wice_soisno(j) = 0. ! prevent tiny residuals
            if (wliq_soisno(j) < 1.e-12) wliq_soisno(j) = 0. ! prevent tiny residuals
            cv_soisno(j) = cv_soisno(j) + melt*(cpliq-cpice) ! update heat capacity
            t_soisno(j) = tfrz + heatrem/cv_soisno(j)
            xmf = xmf + melt*hfus
         end if
      end do
      !------------------------------------------------------------

      IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
         !for SNICAR: layer freezing mass flux (positive):
         DO j = lb, 0
            IF (imelt_soisno(j)==2 .and. j<1) THEN
               snofrz(j) = max(0._r8,(wice_soisno(j)-wice_soisno_bef(j)))/deltim
            ENDIF
         ENDDO
      ENDIF

#if (defined CoLMDEBUG)
      ! second energy check and water check. now check energy balance before and after phase
      ! change, considering the possibility of changed heat capacity during phase change, by
      ! using initial heat capacity in the first step, final heat capacity in the second step,
      ! and differences from tfrz only to avoid enthalpy correction for (cpliq-cpice)*melt*tfrz.
      ! also check soil water sum.
      do j = 1, nl_lake
         esum2 = esum2 - (t_lake(j)-tfrz)*cv_lake(j)
      end do

      do j = lb, nl_soil
         esum2 = esum2 - (t_soisno(j)-tfrz)*cv_soisno(j)
      end do

      esum2 = esum2 - xmf
      errsoi = esum2/deltim

      if(abs(errsoi) > 0.1) then
      write(6,*) 'energy conservation error in LAND WATER COLUMN during phase change, error (W/m^2):', errsoi
      end if

#endif
      
!------------------------------------------------------------
!*[5] Convective mixing: make sure fracice*dz is conserved, heat content c*dz*T is conserved, and
! all ice ends up at the top. Done over all lakes even if frozen.
! Either an unstable density profile or ice in a layer below an incompletely frozen layer will trigger.
!------------------------------------------------------------

      ! recalculate density
      do j = 1, nl_lake
         rhow(j) = (1.-lake_icefrac(j))*1000.*(1.0-1.9549e-05*(abs(t_lake(j)-277.))**1.68) &
                     + lake_icefrac(j)*denice
      end do

      do j = 1, nl_lake-1
         qav = 0.
         nav = 0.
         iceav = 0.

         if (rhow(j)>rhow(j+1) .or. (lake_icefrac(j)<1.0 .and. lake_icefrac(j+1)>0.)) then
            do i = 1, j+1
               qav = qav + dz_lake(i)*(t_lake(i)-tfrz) * &
                       ((1. - lake_icefrac(i))*cwat + lake_icefrac(i)*cice_eff)
               iceav = iceav + lake_icefrac(i)*dz_lake(i)
               nav = nav + dz_lake(i)
            end do

            qav = qav/nav
            iceav = iceav/nav
            !if the average temperature is above freezing, put the extra energy into the water.
            !if it is below freezing, take it away from the ice.
            if (qav > 0.) then
               tav_froz = 0. !Celsius
               tav_unfr = qav / ((1. - iceav)*cwat)
            else if (qav < 0.) then
               tav_froz = qav / (iceav*cice_eff)
               tav_unfr = 0. !Celsius
            else
               tav_froz = 0.
               tav_unfr = 0.
            end if
         end if

         if (nav > 0.) then
            do i = 1, j+1

            !put all the ice at the top.
            !if the average temperature is above freezing, put the extra energy into the water.
            !if it is below freezing, take it away from the ice.
            !for the layer with both ice & water, be careful to use the average temperature
            !that preserves the correct total heat content given what the heat capacity of that
            !layer will actually be.

               if (i == 1) zsum = 0.
               if ((zsum+dz_lake(i))/nav <= iceav) then
                  lake_icefrac(i) = 1.
                  t_lake(i) = tav_froz + tfrz
               else if (zsum/nav < iceav) then
                  lake_icefrac(i) = (iceav*nav - zsum) / dz_lake(i)
                  ! Find average value that preserves correct heat content.
                  t_lake(i) = ( lake_icefrac(i)*tav_froz*cice_eff &
                              + (1. - lake_icefrac(i))*tav_unfr*cwat ) &
                              / ( lake_icefrac(i)*cice_eff + (1-lake_icefrac(i))*cwat ) + tfrz
               else
                  lake_icefrac(i) = 0.
                  t_lake(i) = tav_unfr + tfrz
               end if
               zsum = zsum + dz_lake(i)

               rhow(i) = (1.-lake_icefrac(i))*1000.*(1.-1.9549e-05*(abs(t_lake(i)-277.))**1.68) &
                           + lake_icefrac(i)*denice
            end do
         end if
      end do
      
!------------------------------------------------------------
!*[6] Re-evaluate thermal properties and sum energy content.
!------------------------------------------------------------
      ! for lake
      do j = 1, nl_lake
         cv_lake(j) = dz_lake(j) * (cwat*(1.-lake_icefrac(j)) + cice_eff*lake_icefrac(j))
      end do

      ! do as above to sum energy content
      ncvts = 0.
      do j = 1, nl_lake
         ncvts = ncvts + cv_lake(j)*(t_lake(j)-tfrz) + cfus*dz_lake(j)*(1.-lake_icefrac(j))
      end do

      do j = lb, nl_soil
         ncvts = ncvts + cv_soisno(j)*(t_soisno(j)-tfrz) + hfus*wliq_soisno(j)
         if (j == 1 .and. scv > 0. .and. j == lb) then
            ncvts = ncvts - scv*hfus
         end if
      end do
      ! IF(DEF_USE_TWOSTREAM) THEN
      !    IF (.not. DEF_USE_SNICAR .or. present(urban_call)) THEN
      !       if (snl == 0) then
      !          errsoi = (ncvts-ocvts)/deltim - fgrnd + 2*pi*uavg*topc*rou*(exp(daT*(tao(nl_lake+1)-tao(1)))-exp(-daT*(tao(nl_lake+1)-tao(1))))/D
      !       else 
      !          errsoi = (ncvts-ocvts)/deltim - fgrnd
      !       end if 
      !    ELSE 
      !       errsoi = (ncvts-ocvts)/deltim - fgrnd + 2*pi*uavg*topc*rou*(exp(daT*(tao(nl_lake+1)-tao(1)))-exp(-daT*(tao(nl_lake+1)-tao(1))))/D
      !    END IF
      ! ELSE 
         ! check energy conservation.
         errsoi = (ncvts-ocvts)/deltim - fgrnd 
      ! END IF
      if (abs(errsoi) < 0.10) then
         fseng = fseng - errsoi
         fsena = fseng
         fgrnd = fgrnd + errsoi
         errsoi = 0.
      else
         print*, "energy conservation error in LAND WATER COLUMN during convective mixing", errsoi,fgrnd,ncvts,ocvts
      end if
      write(*,*)'check final'
  end subroutine laketem



  subroutine snowwater_lake ( &
             ! "in" arguments
             ! ---------------------------
             maxsnl      , nl_soil     , nl_lake   , deltim       ,&
             ssi         , wimp        , porsl     , pg_rain      ,&
             pg_snow     , dz_lake     , imelt     , fiold        ,&
             qseva       , qsubl       , qsdew     , qfros        ,&

             ! "inout" arguments
             ! ---------------------------
             z_soisno    , dz_soisno   , zi_soisno , t_soisno     ,&
             wice_soisno , wliq_soisno , t_lake    , lake_icefrac ,&
             qout_snowb  ,                                         &
             fseng       , fgrnd       , snl       , scv          ,&
             snowdp      , sm          , forc_us   , forc_vs      ,&
! SNICAR model variables
             forc_aer    ,&
             mss_bcpho   , mss_bcphi   , mss_ocpho , mss_ocphi    ,&
             mss_dst1    , mss_dst2    , mss_dst3  , mss_dst4     ,&
! END SNICAR model variables
             urban_call  )

!-----------------------------------------------------------------------------------------------
! Calculation of Lake Hydrology. Lake water mass is kept constant. The soil is simply maintained at
! volumetric saturation if ice melting frees up pore space.
!
! Called:
!    -> snowwater:                  change of snow mass and snow water onto soil
!    -> snowcompaction:             compaction of snow layers
!    -> combinesnowlayers:          combine snow layers that are thinner than minimum
!    -> dividesnowlayers:           subdivide snow layers that are thicker than maximum
!
! Initial: Yongjiu Dai, December, 2012
!                          April, 2014
! REVISIONS:
! Nan Wei, 06/2018: update snow hydrology above lake
! Yongjiu Dai, 01/2023: added for SNICAR model effects for snowwater,
! combinesnowlayers, dividesnowlayers processes by calling snowwater_snicar(),
! SnowLayersCombine_snicar, SnowLayersDivide_snicar()
!-----------------------------------------------------------------------------------------------

  use MOD_Precision
  use MOD_Const_Physical, only : denh2o, denice, hfus, tfrz, cpliq, cpice
  use MOD_SoilSnowHydrology
  use MOD_SnowLayersCombineDivide

  implicit none

! ------------- in/inout/out variables -----------------------------------------

  integer, INTENT(in) :: maxsnl  ! maximum number of snow layers
  integer, INTENT(in) :: nl_soil ! number of soil layers
  integer, INTENT(in) :: nl_lake ! number of soil layers

  real(r8), INTENT(in) :: deltim             ! seconds in a time step (sec)
  real(r8), INTENT(in) :: ssi                ! irreducible water saturation of snow
  real(r8), INTENT(in) :: wimp               ! water impremeable if porosity less than wimp
  real(r8), INTENT(in) :: porsl(1:nl_soil)   ! volumetric soil water at saturation (porosity)

  real(r8), INTENT(in) :: pg_rain            ! rainfall incident on ground [mm/s]
  real(r8), INTENT(in) :: pg_snow            ! snowfall incident on ground [mm/s]
  real(r8), INTENT(in) :: dz_lake(1:nl_lake) ! layer thickness for lake (m)

  integer,  INTENT(in) :: imelt(maxsnl+1:0)  ! signifies if node in melting (imelt = 1)
  real(r8), INTENT(in) :: fiold(maxsnl+1:0)  ! fraction of ice relative to the total water content at the previous time step
  real(r8), INTENT(in) :: qseva              ! ground surface evaporation rate (mm h2o/s)
  real(r8), INTENT(in) :: qsubl              ! sublimation rate from snow pack (mm H2O /s) [+]
  real(r8), INTENT(in) :: qsdew              ! surface dew added to snow pack (mm H2O /s) [+]
  real(r8), INTENT(in) :: qfros              ! ground surface frosting formation (mm H2O /s) [+]

  real(r8), INTENT(inout) :: z_soisno   (maxsnl+1:nl_soil) ! layer depth  (m)
  real(r8), INTENT(inout) :: dz_soisno  (maxsnl+1:nl_soil) ! layer thickness depth (m)
  real(r8), INTENT(inout) :: zi_soisno    (maxsnl:nl_soil) ! interface depth (m)
  real(r8), INTENT(inout) :: t_soisno   (maxsnl+1:nl_soil) ! snow temperature (Kelvin)
  real(r8), INTENT(inout) :: wice_soisno(maxsnl+1:nl_soil) ! ice lens (kg/m2)
  real(r8), INTENT(inout) :: wliq_soisno(maxsnl+1:nl_soil) ! liquid water (kg/m2)
  real(r8), INTENT(inout) :: t_lake      (1:nl_lake) ! lake temperature (Kelvin)
  real(r8), INTENT(inout) :: lake_icefrac(1:nl_lake) ! mass fraction of lake layer that is frozen
  real(r8), INTENT(inout) :: qout_snowb ! rate of water out of snow bottom (mm/s)

  real(r8), INTENT(inout) :: fseng  ! total sensible heat flux (W/m**2) [+ to atm]
  real(r8), INTENT(inout) :: fgrnd  ! heat flux into snow / lake (W/m**2) [+ = into soil]

  integer , INTENT(inout) :: snl    ! number of snow layers
  real(r8), INTENT(inout) :: scv    ! snow water (mm H2O)
  real(r8), INTENT(inout) :: snowdp ! snow height (m)
  real(r8), INTENT(inout) :: sm     ! rate of snow melt (mm H2O /s)

  real(r8), intent(in) :: forc_us
  real(r8), intent(in) :: forc_vs

! SNICAR model variables
! Aerosol Fluxes (Jan. 07, 2023 by Yongjiu Dai)
  real(r8), intent(in) :: forc_aer ( 14 )  ! aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]

  logical, optional, intent(in) :: urban_call   ! whether it is a urban CALL

  real(r8), INTENT(inout) :: &
        mss_bcpho (maxsnl+1:0), &! mass of hydrophobic BC in snow  (col,lyr) [kg]
        mss_bcphi (maxsnl+1:0), &! mass of hydrophillic BC in snow (col,lyr) [kg]
        mss_ocpho (maxsnl+1:0), &! mass of hydrophobic OC in snow  (col,lyr) [kg]
        mss_ocphi (maxsnl+1:0), &! mass of hydrophillic OC in snow (col,lyr) [kg]
        mss_dst1  (maxsnl+1:0), &! mass of dust species 1 in snow  (col,lyr) [kg]
        mss_dst2  (maxsnl+1:0), &! mass of dust species 2 in snow  (col,lyr) [kg]
        mss_dst3  (maxsnl+1:0), &! mass of dust species 3 in snow  (col,lyr) [kg]
        mss_dst4  (maxsnl+1:0)   ! mass of dust species 4 in snow  (col,lyr) [kg]
! Aerosol Fluxes (Jan. 07, 2023)
! END SNICAR model variables

! ------------- other local variables -----------------------------------------
  integer  j          ! indices
  integer lb          ! lower bound of array

  real(r8) xmf        ! snow melt heat flux (W/m**2)

  real(r8) sumsnowice ! sum of snow ice if snow layers found above unfrozen lake [kg/m&2]
  real(r8) sumsnowliq ! sum of snow liquid if snow layers found above unfrozen lake [kg/m&2]
  logical  unfrozen   ! true if top lake layer is unfrozen with snow layers above
  real(r8) heatsum    ! used in case above [J/m^2]
  real(r8) heatrem    ! used in case above [J/m^2]

  real(r8) a, b, c, d
  real(r8) wice_lake(1:nl_lake)  ! ice lens (kg/m2)
  real(r8) wliq_lake(1:nl_lake)  ! liquid water (kg/m2)
  real(r8) t_ave, frac_
!-----------------------------------------------------------------------

      ! for runoff calculation (assumed no mass change in the land water bodies)
      lb = snl + 1
      qout_snowb = 0.0

      ! ----------------------------------------------------------
      !*[1] snow layer on frozen lake
      ! ----------------------------------------------------------
      if (snl < 0) then
         lb = snl + 1
         IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
            call snowwater_SNICAR (lb,deltim,ssi,wimp,&
                         pg_rain,qseva,qsdew,qsubl,qfros,&
                         dz_soisno(lb:0),wice_soisno(lb:0),wliq_soisno(lb:0),qout_snowb,    &
                         forc_aer,&
                         mss_bcpho(lb:0), mss_bcphi(lb:0), mss_ocpho(lb:0), mss_ocphi(lb:0),&
                         mss_dst1(lb:0),  mss_dst2(lb:0),  mss_dst3(lb:0),  mss_dst4(lb:0) )
         ELSE
            call snowwater (lb,deltim,ssi,wimp,&
                         pg_rain,qseva,qsdew,qsubl,qfros,&
                         dz_soisno(lb:0),wice_soisno(lb:0),wliq_soisno(lb:0),qout_snowb)
         ENDIF

         ! Natural compaction and metamorphosis.
         lb = snl + 1
         call snowcompaction (lb,deltim, &
                              imelt(lb:0),fiold(lb:0),t_soisno(lb:0),&
                              wliq_soisno(lb:0),wice_soisno(lb:0),forc_us,forc_vs,dz_soisno(lb:0))

         ! Combine thin snow elements
         lb = maxsnl + 1
         IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
            call snowlayerscombine_SNICAR (lb, snl,&
                                 z_soisno(lb:1),dz_soisno(lb:1),zi_soisno(lb-1:0),&
                                 wliq_soisno(lb:1),wice_soisno(lb:1), t_soisno(lb:1),scv,snowdp, &
                                 mss_bcpho(lb:0), mss_bcphi(lb:0), mss_ocpho(lb:0), mss_ocphi(lb:0),&
                                 mss_dst1(lb:0),  mss_dst2(lb:0),  mss_dst3(lb:0),  mss_dst4(lb:0))
         ELSE
            call snowlayerscombine (lb, snl,&
                                 z_soisno(lb:1),dz_soisno(lb:1),zi_soisno(lb-1:0),&
                                 wliq_soisno(lb:1),wice_soisno(lb:1),&
                                 t_soisno(lb:1),scv,snowdp)
         ENDIF

         ! Divide thick snow elements
         if (snl < 0) then
            IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
               call snowlayersdivide_SNICAR (lb,snl,z_soisno(lb:0),dz_soisno(lb:0),zi_soisno(lb-1:0),&
                                 wliq_soisno(lb:0),wice_soisno(lb:0),t_soisno(lb:0)     ,&
                                 mss_bcpho(lb:0), mss_bcphi(lb:0), mss_ocpho(lb:0), mss_ocphi(lb:0),&
                                 mss_dst1(lb:0),  mss_dst2(lb:0),  mss_dst3(lb:0),  mss_dst4(lb:0) )
            ELSE
               call snowlayersdivide (lb,snl,z_soisno(lb:0),dz_soisno(lb:0),zi_soisno(lb-1:0),&
                                 wliq_soisno(lb:0),wice_soisno(lb:0),t_soisno(lb:0))
            ENDIF
         endif

      ! ----------------------------------------------------------
      !*[2] check for single completely unfrozen snow layer over lake.
      !     Modeling this ponding is unnecessary and can cause instability after the timestep
      !     when melt is completed, as the temperature after melt can be excessive
      !     because the fluxes were calculated with a fixed ground temperature of freezing, but the
      !     phase change was unable to restore the temperature to freezing.  (Zack Subnin 05/2010)
      ! ----------------------------------------------------------

         if (snl == -1 .and. wice_soisno(0) == 0.) then
            ! Remove layer
            ! Take extra heat of layer and release to sensible heat in order to maintain energy conservation.
            heatrem = cpliq*wliq_soisno(0)*(t_soisno(0) - tfrz)
            fseng = fseng + heatrem/deltim
            fgrnd = fgrnd - heatrem/deltim

            snl = 0
            scv = 0.
            snowdp = 0.
         end if

      endif

      ! ----------------------------------------------------------
      !*[3] check for snow layers above lake with unfrozen top layer. Mechanically,
      !     the snow will fall into the lake and melt or turn to ice. If the top layer has
      !     sufficient heat to melt the snow without freezing, then that will be done.
      !     Otherwise, the top layer will undergo freezing, but only if the top layer will
      !     not freeze completely. Otherwise, let the snow layers persist and melt by diffusion.
      ! ----------------------------------------------------------

      if (t_lake(1) > tfrz .and. snl < 0 .and. lake_icefrac(1) < 0.001) then ! for unfrozen lake
         unfrozen = .true.
      else
         unfrozen = .false.
      end if

      sumsnowice = 0.
      sumsnowliq = 0.
      heatsum = 0.0
      do j = snl+1,0
         if (unfrozen) then
            sumsnowice = sumsnowice + wice_soisno(j)
            sumsnowliq = sumsnowliq + wliq_soisno(j)
            heatsum = heatsum + wice_soisno(j)*cpice*(tfrz-t_soisno(j)) &
                             + wliq_soisno(j)*cpliq*(tfrz-t_soisno(j))
         end if
      end do

      if (unfrozen) then
         ! changed by weinan as the subroutine newsnow_lake
         ! Remove snow and subtract the latent heat from the top layer.

         t_ave = tfrz - heatsum/(sumsnowice*cpice + sumsnowliq*cpliq)

         a = heatsum
         b = sumsnowice*hfus
         c = (t_lake(1) - tfrz)*cpliq*denh2o*dz_lake(1)
         d = denh2o*dz_lake(1)*hfus

         ! all snow melt
           if (c>=a+b)then
              t_lake(1) = (cpliq*(denh2o*dz_lake(1)*t_lake(1) + (sumsnowice+sumsnowliq)*tfrz) - a - b) / &
                        (cpliq*(denh2o*dz_lake(1) + sumsnowice+ sumsnowice))
              sm = sm + scv/deltim
              scv = 0.
              snowdp = 0.
              snl = 0
        ! lake partially freezing to melt all snow
           else if(c+d >= a+b)then
              t_lake(1) = tfrz
              sm = sm + scv/deltim
              scv = 0.
              snowdp = 0.
              snl = 0
              lake_icefrac(1) = (a+b-c)/d

        !  snow do not melt while all lake freezing
        !    else if(c+d < a) then
        !     t_lake(1) = (c+d + cpice*(sumsnowice*t_ave+denh2o*dz_lake(1)*tfrz) + cpliq*sumsnowliq*t_ave)/&
        !                (cpice*(sumsnowice+denh2o*dz_lake(1))+cpliq*sumsnowliq)
        !     lake_icefrac(1) = 1.0
            end if
      end if


      ! ----------------------------------------------------------
      !*[4] Soil water and ending water balance
      ! ----------------------------------------------------------
      ! Here this consists only of making sure that soil is saturated even as it melts and
      ! pore space opens up. Conversely, if excess ice is melting and the liquid water exceeds the
      ! saturation value, then remove water.

      do j = 1, nl_soil
         a = wliq_soisno(j)/(dz_soisno(j)*denh2o) + wice_soisno(j)/(dz_soisno(j)*denice)

         if (a < porsl(j)) then
            wliq_soisno(j) = max( 0., (porsl(j)*dz_soisno(j) - wice_soisno(j)/denice)*denh2o )
            wice_soisno(j) = max( 0., (porsl(j)*dz_soisno(j) - wliq_soisno(j)/denh2o)*denice )
         else
            wliq_soisno(j) = max(0., wliq_soisno(j) - (a - porsl(j))*denh2o*dz_soisno(j) )
            wice_soisno(j) = max( 0., (porsl(j)*dz_soisno(j) - wliq_soisno(j)/denh2o)*denice )
         end if

         if (wliq_soisno(j) > porsl(j)*denh2o*dz_soisno(j)) then
             wliq_soisno(j) = porsl(j)*denh2o*dz_soisno(j)
             wice_soisno(j) = 0.0
         endif
      end do


  end subroutine snowwater_lake



  subroutine roughness_lake (snl,t_grnd,t_lake,lake_icefrac,forc_psrf,&
                             cur,ustar,z0mg,z0hg,z0qg)

!-----------------------------------------------------------------------
! DESCRIPTION:
! Calculate lake surface roughness
!
! Original:
! The Community Land Model version 4.5 (CLM4.5)
!
! Revisions:
! Yongjiu Dai, Nan Wei, 01/2018
!-----------------------------------------------------------------------

  use MOD_Precision
  use MOD_Const_Physical, only : tfrz,vonkar,grav

  IMPLICIT NONE

  integer,  INTENT(in) :: snl       ! number of snow layers
  real(r8), INTENT(in) :: t_grnd    ! ground temperature
  real(r8), INTENT(in) :: t_lake(1) ! surface lake layer temperature [K]
  real(r8), INTENT(in) :: lake_icefrac(1) ! surface lake layer ice mass fraction [0-1]
  real(r8), INTENT(in) :: forc_psrf ! atmosphere pressure at the surface [pa]

  real(r8), INTENT(in) :: cur       ! Charnock parameter (-)
  real(r8), INTENT(in) :: ustar     ! u* in similarity theory [m/s]

  real(r8), INTENT(out) :: z0mg     ! roughness length over ground, momentum [m]
  real(r8), INTENT(out) :: z0hg     ! roughness length over ground, sensible heat [m]
  real(r8), INTENT(out) :: z0qg     ! roughness length over ground, latent heat [m]

  real(r8), parameter :: cus = 0.1       ! empirical constant for roughness under smooth flow
  real(r8), parameter :: kva0 = 1.51e-5  ! kinematic viscosity of air (m^2/s) at 20C and 1.013e5 Pa
  real(r8), parameter :: prn = 0.713     ! Prandtl # for air at neutral stability
  real(r8), parameter :: sch = 0.66      ! Schmidt # for water in air at neutral stability

  real(r8) kva    ! kinematic viscosity of air at ground temperature and forcing pressure
  real(r8) sqre0  ! root of roughness Reynolds number

      if (t_grnd > tfrz .and. t_lake(1) > tfrz .and. snl == 0) then
          kva = kva0 * (t_grnd/293.15)**1.5 * 1.013e5/forc_psrf ! kinematic viscosity of air
          z0mg = max(cus*kva/max(ustar,1.e-4),cur*ustar*ustar/grav) ! momentum roughness length
          z0mg = max(z0mg, 1.0e-5) ! This limit is redundant with current values.
          sqre0 = (max(z0mg*ustar/kva,0.1))**0.5   ! square root of roughness Reynolds number
          z0hg = z0mg * exp( -vonkar/prn*( 4.*sqre0 - 3.2) ) ! SH roughness length
          z0qg = z0mg * exp( -vonkar/sch*( 4.*sqre0 - 4.2) ) ! LH roughness length
          z0qg = max(z0qg, 1.0e-5)  ! Minimum allowed roughness length for unfrozen lakes
          z0hg = max(z0hg, 1.0e-5)  ! set low so it is only to avoid floating point exceptions
      else if (snl == 0) then ! frozen lake with ice, and no snow cover
          z0mg = 0.001              ! z0mg won't have changed
          z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
          z0qg = z0hg
      else                          ! use roughness over snow
          z0mg = 0.0024             ! z0mg won't have changed
          z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
          z0qg = z0hg
      end if

  end subroutine roughness_lake



  subroutine hConductivity_lake(nl_lake,snl,t_grnd,&
                                z_lake,t_lake,lake_icefrac,rhow,&
                                dlat,ustar,z0mg,lakedepth,depthcrit,tk_lake, savedtke1)

! -------------------------------------------------------------------------
! Diffusivity and implied thermal "conductivity" = diffusivity * cwat
! -------------------------------------------------------------------------

  use MOD_Precision
  use MOD_Const_Physical, only : tfrz,tkwat,tkice,tkair,&
                                vonkar,grav,cpliq,cpice,cpair,denh2o,denice

  IMPLICIT NONE

  integer, INTENT(in) :: nl_lake  ! number of soil layers
  integer, INTENT(in) :: snl      ! number of snow layers
  real(r8), INTENT(in) :: t_grnd  ! ground surface temperature [k]
  real(r8), INTENT(in) :: z_lake(nl_lake)       ! lake node depth (middle point of layer) (m)
  real(r8), INTENT(in) :: t_lake(nl_lake)       ! lake temperature (kelvin)
  real(r8), INTENT(in) :: lake_icefrac(nl_lake) ! lake mass fraction of lake layer that is frozen
  real(r8), INTENT(in) :: rhow(nl_lake)         ! density of water (kg/m**3)

  real(r8), INTENT(in) :: dlat      ! latitude (radians)
  real(r8), INTENT(in) :: ustar     ! u* in similarity theory [m/s]
  real(r8), INTENT(in) :: z0mg      ! roughness length over ground, momentum [m]
  real(r8), INTENT(in) :: lakedepth ! column lake depth (m)
  real(r8), INTENT(in) :: depthcrit ! (m) Depth beneath which to enhance mixing


  real(r8), INTENT(out) :: tk_lake(nl_lake) ! thermal conductivity at layer node [W/(m K)]
  real(r8), INTENT(out) :: savedtke1      ! top level eddy conductivity (W/mK)

! local
  real(r8) kme(nl_lake) ! molecular + eddy diffusion coefficient (m**2/s)
  real(r8) cwat   ! specific heat capacity of water (j/m**3/kelvin)
  real(r8) den    ! used in calculating ri
  real(r8) drhodz ! d [rhow] /dz (kg/m**4)
  real(r8) fangkm ! (m^2/s) extra diffusivity based on Fang & Stefan 1996
  real(r8) ke     ! eddy diffusion coefficient (m**2/s)
  real(r8) km     ! molecular diffusion coefficient (m**2/s)
  real(r8) ks     ! coefficient for calculation of decay of eddy diffusivity with depth
  real(r8) n2     ! brunt-vaisala frequency (/s**2)
  real(r8) num    ! used in calculating ri
  real(r8) ri     ! richardson number
  real(r8) tkice_eff  ! effective conductivity since layer depth is constant
  real(r8) tmp    !
  real(r8) u2m    ! 2 m wind speed (m/s
  real(r8) ws     ! surface friction velocity (m/s)

  real(r8), parameter :: mixfact = 5. ! Mixing enhancement factor.
  real(r8), parameter :: p0 = 1.      ! neutral value of turbulent prandtl number


  integer j

! -------------------------------------------------------------------

      cwat = cpliq*denh2o
      tkice_eff = tkice * denice/denh2o ! effective conductivity since layer depth is constant
      km = tkwat/cwat                   ! a constant (molecular diffusivity)
      u2m = max(0.1,ustar/vonkar*log(2./z0mg))
      ws = 1.2e-03 * u2m
      ks = 6.6 * sqrt( abs(sin(dlat)) ) * (u2m**(-1.84))

      do j = 1, nl_lake-1
         drhodz = (rhow(j+1)-rhow(j)) / (z_lake(j+1)-z_lake(j))
         n2 = max(7.5e-5, grav / rhow(j) * drhodz)
         num = 40. * n2 * (vonkar*z_lake(j))**2
         tmp = -2.*ks*z_lake(j)        ! to avoid underflow computing
         if(tmp < -40.) tmp = -40.     !
         den = max( (ws**2) * exp(tmp), 1.e-10 )
         ri = ( -1. + sqrt( max(1.+num/den, 0.) ) ) / 20.

         if ((t_grnd > tfrz .and. t_lake(1) > tfrz .and. snl == 0) ) then
            tmp = -ks*z_lake(j)        ! to avoid underflow computing
            if(tmp < -40.) tmp = -40.  !
            ke = vonkar*ws*z_lake(j)/p0 * exp(tmp) / (1.+37.*ri*ri)
            kme(j) = km + ke

            fangkm = 1.039e-8_r8 * max(n2,7.5e-5)**(-0.43)  ! Fang & Stefan 1996, citing Ellis et al 1991
            kme(j) = kme(j) + fangkm

            if (lakedepth >= depthcrit) then
               kme(j) = kme(j) * mixfact    ! Mixing enhancement factor for lake deep than 25m.
            end if
            tk_lake(j) = kme(j)*cwat
         else
            kme(j) = km
            fangkm = 1.039e-8 * max(n2,7.5e-5)**(-0.43)
            kme(j) = kme(j) + fangkm
            if (lakedepth >= depthcrit) then
                kme(j) = kme(j) * mixfact
            end if
            tk_lake(j) = kme(j)*cwat*tkice_eff / ((1.-lake_icefrac(j))*tkice_eff &
                       + kme(j)*cwat*lake_icefrac(j))
         end if
      end do

      kme(nl_lake) = kme(nl_lake-1)
       savedtke1 = kme(1)*cwat

      if ((t_grnd > tfrz .and. t_lake(1) > tfrz .and. snl == 0) ) then
         tk_lake(nl_lake) = tk_lake(nl_lake-1)
      else
         tk_lake(nl_lake) = kme(nl_lake)*cwat*tkice_eff / ( (1.-lake_icefrac(nl_lake))*tkice_eff &
                    + kme(nl_lake)*cwat*lake_icefrac(nl_lake) )
      end if

  end subroutine hConductivity_lake

  !--------------------------------------------------------------------------------------
  subroutine LakeOptics_init( fsnowoptics, foptd, foptr )

   USE MOD_NetCDFSerial

   IMPLICIT NONE

   character(len=256), intent(in) :: fsnowoptics     ! snow optical properties file name
   character(len=256), intent(in) :: foptd     ! direct optical properties file name
   character(len=256), intent(in) :: foptr     ! diffuse optical properties file name
   character(len= 32) :: subname = 'LakeOptics_init' ! subroutine name
   integer :: atm_type_index                         ! index for atmospheric type

   logical :: readvar  ! determine if variable was read from NetCDF file
!-----------------------------------------------------------------------

   readvar = .true.

   atm_type_index = atm_type_mid_latitude_winter
   ! Define atmospheric type
   if (trim(snicar_atm_type) == 'default') then
     atm_type_index = atm_type_default
   elseif (trim(snicar_atm_type) == 'mid-latitude_winter') then
     atm_type_index = atm_type_mid_latitude_winter
   elseif (trim(snicar_atm_type) == 'mid-latitude_summer') then
     atm_type_index = atm_type_mid_latitude_summer
   elseif (trim(snicar_atm_type) == 'sub-Arctic_winter') then
     atm_type_index = atm_type_sub_Arctic_winter
   elseif (trim(snicar_atm_type) == 'sub-Arctic_summer') then
     atm_type_index = atm_type_sub_Arctic_summer
   elseif (trim(snicar_atm_type) == 'summit_Greenland') then
     atm_type_index = atm_type_summit_Greenland
   elseif (trim(snicar_atm_type) == 'high_mountain') then
     atm_type_index = atm_type_high_mountain
   ! else
   !    IF (p_is_master) THEN
   !       write(iulog,*) "snicar_atm_type = ", snicar_atm_type
   !       call abort
   !    ENDIF
   endif

   !
   ! Open optics file:
   ! IF (p_is_master) THEN
   !    write(iulog,*) 'Attempting to read snow optical properties .....'
   !    write(iulog,*) subname,trim(fsnowoptics)
   ! ENDIF

   ! direct-beam snow and ice Mie parameters:
   CALL ncio_read_bcast_serial (foptr, 'ss_alb_snw_avg', ss_alb_snw_avg_clr)
   CALL ncio_read_bcast_serial (foptr, 'asm_prm_snw_avg', asm_prm_snw_avg_clr)
   CALL ncio_read_bcast_serial (foptr, 'ext_cff_mss_snw_avg', ext_cff_mss_snw_avg_clr)
   CALL ncio_read_bcast_serial (foptr, 'sca_cff_vlm_avg', sca_cff_vlm_avg_clr)
   CALL ncio_read_bcast_serial (foptr, 'asm_prm_ice_avg', asm_prm_ice_avg_clr)
   CALL ncio_read_bcast_serial (foptr, 'abs_cff_mss_avg', abs_cff_mss_avg_clr)
   CALL ncio_read_bcast_serial (foptr, 'ext_cff_mss_wtr_avg', ext_cff_mss_wtr_avg_clr)
   CALL ncio_read_bcast_serial (foptr, 'ss_alb_wtr_avg', ss_alb_wtr_avg_clr)

   ! direct-beam ice and water complex refractive index:
   CALL ncio_read_bcast_serial (foptr, 'rfidx_ice_im_avg_Pic16', refindx_im_clr)
   CALL ncio_read_bcast_serial (foptr, 'rfidx_ice_re_avg_Pic16', refindx_re_clr)
   CALL ncio_read_bcast_serial (foptr, 'rfidx_wtr_im_avg', refindxwat_im_clr)
   CALL ncio_read_bcast_serial (foptr, 'rfidx_wtr_re_avg', refindxwat_re_clr)

   ! diffuse snow and ice Mie parameters
   CALL ncio_read_bcast_serial (foptd, 'ss_alb_snw_avg', ss_alb_snw_avg)
   CALL ncio_read_bcast_serial (foptd, 'asm_prm_snw_avg', asm_prm_snw_avg)
   CALL ncio_read_bcast_serial (foptd, 'ext_cff_mss_snw_avg', ext_cff_mss_snw_avg)
   CALL ncio_read_bcast_serial (foptd, 'sca_cff_vlm_avg', sca_cff_vlm_avg)
   CALL ncio_read_bcast_serial (foptd, 'asm_prm_ice_avg', asm_prm_ice_avg)
   CALL ncio_read_bcast_serial (foptd, 'abs_cff_mss_avg', abs_cff_mss_avg)

   ! direct-beam ice and water complex refractive index:
   CALL ncio_read_bcast_serial (foptd, 'rfidx_ice_im_avg_Pic16', refindx_im_cld)
   CALL ncio_read_bcast_serial (foptd, 'rfidx_ice_re_avg_Pic16', refindx_re_cld)
   CALL ncio_read_bcast_serial (foptd, 'rfidx_wtr_im_avg', refindxwat_im_cld)
   CALL ncio_read_bcast_serial (foptd, 'rfidx_wtr_re_avg', refindxwat_re_cld) 

   !!! Direct and diffuse flux under different atmospheric conditions
   ! Direct-beam incident spectral flux:
   CALL ncio_read_bcast_serial (fsnowoptics, 'flx_wgt_dir', flx_wgt_dir)

   ! Diffuse incident spectral flux:
   CALL ncio_read_bcast_serial (fsnowoptics, 'flx_wgt_dif', flx_wgt_dif)

#ifdef MODAL_AER
  ! size-dependent BC parameters and BC enhancement factors
!   IF (p_is_master) THEN
!      write(iulog,*) 'Attempting to read optical properties for within-ice BC (modal aerosol treatment) ...'
!   ENDIF
  !
  ! BC species 1 Mie parameters
  !
  CALL ncio_read_bcast_serial (fsnowoptics, 'ss_alb_bc_mam', ss_alb_bc1)
  CALL ncio_read_bcast_serial (fsnowoptics, 'asm_prm_bc_mam', asm_prm_bc1)
  CALL ncio_read_bcast_serial (fsnowoptics, 'ext_cff_mss_bc_mam', ext_cff_mss_bc1)
  !
  ! BC species 2 Mie parameters (identical, before enhancement factors applied)
  !
  CALL ncio_read_bcast_serial (fsnowoptics, 'ss_alb_bc_mam', ss_alb_bc2)
  CALL ncio_read_bcast_serial (fsnowoptics, 'asm_prm_bc_mam', asm_prm_bc2)
  CALL ncio_read_bcast_serial (fsnowoptics, 'ext_cff_mss_bc_mam', ext_cff_mss_bc2)
  !
  ! size-dependent BC absorption enhancement factors for within-ice BC
  CALL ncio_read_bcast_serial (fsnowoptics, 'bcint_enh_mam', bcenh)
  !
#else
  ! bulk aerosol treatment
  ! BC species 1 Mie parameters
  CALL ncio_read_bcast_serial (fsnowoptics, 'ss_alb_bcphil', ss_alb_bc1)
  CALL ncio_read_bcast_serial (fsnowoptics, 'asm_prm_bcphil', asm_prm_bc1)
  CALL ncio_read_bcast_serial (fsnowoptics, 'ext_cff_mss_bcphil', ext_cff_mss_bc1)

  !
  ! BC species 2 Mie parameters
  !
  CALL ncio_read_bcast_serial (fsnowoptics, 'ss_alb_bcphob', ss_alb_bc2)
  CALL ncio_read_bcast_serial (fsnowoptics, 'asm_prm_bcphob', asm_prm_bc2)
  CALL ncio_read_bcast_serial (fsnowoptics, 'ext_cff_mss_bcphob', ext_cff_mss_bc2)
  !
#endif
  !
  ! OC species 1 Mie parameters
  CALL ncio_read_bcast_serial (fsnowoptics, 'ss_alb_ocphil', ss_alb_oc1)
  CALL ncio_read_bcast_serial (fsnowoptics, 'asm_prm_ocphil', asm_prm_oc1)
  CALL ncio_read_bcast_serial (fsnowoptics, 'ext_cff_mss_ocphil', ext_cff_mss_oc1)
  !
  ! OC species 2 Mie parameters
  CALL ncio_read_bcast_serial (fsnowoptics, 'ss_alb_ocphob', ss_alb_oc2)
  CALL ncio_read_bcast_serial (fsnowoptics, 'asm_prm_ocphob', asm_prm_oc2)
  CALL ncio_read_bcast_serial (fsnowoptics, 'ext_cff_mss_ocphob', ext_cff_mss_oc2)
  !
  ! dust species 1 Mie parameters
  CALL ncio_read_bcast_serial (fsnowoptics, 'ss_alb_dust01', ss_alb_dst1)
  CALL ncio_read_bcast_serial (fsnowoptics, 'asm_prm_dust01', asm_prm_dst1)
  CALL ncio_read_bcast_serial (fsnowoptics, 'ext_cff_mss_dust01', ext_cff_mss_dst1)
  !
  ! dust species 2 Mie parameters
  CALL ncio_read_bcast_serial (fsnowoptics, 'ss_alb_dust02', ss_alb_dst2)
  CALL ncio_read_bcast_serial (fsnowoptics, 'asm_prm_dust02', asm_prm_dst2)
  CALL ncio_read_bcast_serial (fsnowoptics, 'ext_cff_mss_dust02', ext_cff_mss_dst2)
  !
  ! dust species 3 Mie parameters
  CALL ncio_read_bcast_serial (fsnowoptics, 'ss_alb_dust03', ss_alb_dst3)
  CALL ncio_read_bcast_serial (fsnowoptics, 'asm_prm_dust03', asm_prm_dst3)
  CALL ncio_read_bcast_serial (fsnowoptics, 'ext_cff_mss_dust03', ext_cff_mss_dst3)
  !
  ! dust species 4 Mie parameters
  CALL ncio_read_bcast_serial (fsnowoptics, 'ss_alb_dust04', ss_alb_dst4)
  CALL ncio_read_bcast_serial (fsnowoptics, 'asm_prm_dust04', asm_prm_dst4)
  CALL ncio_read_bcast_serial (fsnowoptics, 'ext_cff_mss_dust04', ext_cff_mss_dst4)
  !
  !

!   IF (p_is_master) THEN
!      write(iulog,*) 'Successfully read snow optical properties'
!   ENDIF


  ! print some diagnostics:
!   IF (p_is_master) THEN
!      write (iulog,*) 'SNICAR: Mie single scatter albedos for direct-beam ice, rds=100um: ', &
!         ss_alb_snw_drc(71,1), ss_alb_snw_drc(71,2), ss_alb_snw_drc(71,3),     &
!         ss_alb_snw_drc(71,4), ss_alb_snw_drc(71,5)
!      write (iulog,*) 'SNICAR: Mie single scatter albedos for diffuse ice, rds=100um: ',     &
!         ss_alb_snw_dfs(71,1), ss_alb_snw_dfs(71,2), ss_alb_snw_dfs(71,3),     &
!         ss_alb_snw_dfs(71,4), ss_alb_snw_dfs(71,5)
!      if (DO_SNO_OC) then
!         write (iulog,*) 'SNICAR: Including OC aerosols from snow radiative transfer calculations'
!      else
!         write (iulog,*) 'SNICAR: Excluding OC aerosols from snow radiative transfer calculations'
!      endif
!   ENDIF
  !
! #ifdef MODAL_AER
! !   IF (p_is_master) THEN
! !      ! unique dimensionality for modal aerosol optical properties
! !      write (iulog,*) 'SNICAR: Subset of Mie single scatter albedos for BC: ', &
! !         ss_alb_bc1(1,1), ss_alb_bc1(1,2), ss_alb_bc1(2,1), ss_alb_bc1(5,1), ss_alb_bc1(1,10), ss_alb_bc2(1,10)
! !      write (iulog,*) 'SNICAR: Subset of Mie mass extinction coefficients for BC: ', &
! !         ext_cff_mss_bc2(1,1), ext_cff_mss_bc2(1,2), ext_cff_mss_bc2(2,1), ext_cff_mss_bc2(5,1), ext_cff_mss_bc2(1,10),&
! !         ext_cff_mss_bc1(1,10)
! !      write (iulog,*) 'SNICAR: Subset of Mie asymmetry parameters for BC: ', &
! !         asm_prm_bc1(1,1), asm_prm_bc1(1,2), asm_prm_bc1(2,1), asm_prm_bc1(5,1), asm_prm_bc1(1,10), asm_prm_bc2(1,10)
! !      write (iulog,*) 'SNICAR: Subset of BC absorption enhancement factors: ', &
! !         bcenh(1,1,1), bcenh(1,2,1), bcenh(1,1,2), bcenh(2,1,1), bcenh(5,10,1), bcenh(5,1,8), bcenh(5,10,8)
! !   ENDIF
! #else
! !   IF (p_is_master) THEN
! !      write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophillic BC: ', &
! !         ss_alb_bc1(1), ss_alb_bc1(2), ss_alb_bc1(3), ss_alb_bc1(4), ss_alb_bc1(5)
! !      write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophobic BC: ', &
! !         ss_alb_bc2(1), ss_alb_bc2(2), ss_alb_bc2(3), ss_alb_bc2(4), ss_alb_bc2(5)
! !   ENDIF
! #endif

!   IF (p_is_master) THEN
!      if (DO_SNO_OC) then
!         write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophillic OC: ', &
!            ss_alb_oc1(1), ss_alb_oc1(2), ss_alb_oc1(3), ss_alb_oc1(4), ss_alb_oc1(5)
!         write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophobic OC: ', &
!            ss_alb_oc2(1), ss_alb_oc2(2), ss_alb_oc2(3), ss_alb_oc2(4), ss_alb_oc2(5)
!      endif

!      write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 1: ', &
!         ss_alb_dst1(1), ss_alb_dst1(2), ss_alb_dst1(3), ss_alb_dst1(4), ss_alb_dst1(5)
!      write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 2: ', &
!         ss_alb_dst2(1), ss_alb_dst2(2), ss_alb_dst2(3), ss_alb_dst2(4), ss_alb_dst2(5)
!      write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 3: ', &
!         ss_alb_dst3(1), ss_alb_dst3(2), ss_alb_dst3(3), ss_alb_dst3(4), ss_alb_dst3(5)
!      write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 4: ', &
!         ss_alb_dst4(1), ss_alb_dst4(2), ss_alb_dst4(3), ss_alb_dst4(4), ss_alb_dst4(5)
!      write(iulog,*)
!   ENDIF

end subroutine LakeOptics_init
!-----------------------------------------------------------------------

END MODULE MOD_Lake
