module MossLichenMod

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Functionalities for moss and lichen that can be used for CLM and FATES
  ! SoilFluxes then determines soil/snow and ground temperatures and updates the surface
  ! fluxes for the new ground temperature.
  !
  ! !USES:
  use shr_sys_mod           , only : shr_sys_flush
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use abortutils            , only : endrun
  use clm_varctl            , only : iulog, use_cn, use_lch4, use_c13, use_c14, use_cndv, use_fates, &
                                     use_luna, use_hydrstress
  use clm_varpar            , only : nlevgrnd, nlevsno
  use clm_varcon            , only : namep
  use pftconMod             , only : pftcon
  use decompMod             , only : bounds_type
  use ActiveLayerMod        , only : active_layer_type
  use PhotosynthesisMod     , only : Photosynthesis, PhotoSynthesisHydraulicStress, PhotosynthesisTotal, Fractionation
  use EDAccumulateFluxesMod , only : AccumulateFluxes_ED
  use SoilMoistStressMod    , only : calc_effective_soilporosity, calc_volumetric_h2oliq
  use SoilMoistStressMod    , only : calc_root_moist_stress, set_perchroot_opt
  use SimpleMathMod         , only : array_div_vector
  use SurfaceResistanceMod  , only : do_soilevap_beta,do_soil_resistance_sl14
  use atm2lndType           , only : atm2lnd_type
  use CanopyStateType       , only : canopystate_type
  use EnergyFluxType        , only : energyflux_type
  use FrictionvelocityMod   , only : frictionvel_type
  use OzoneBaseMod          , only : ozone_base_type
  use SoilStateType         , only : soilstate_type
  use SolarAbsorbedType     , only : solarabs_type
  use SurfaceAlbedoType     , only : surfalb_type
  use TemperatureType       , only : temperature_type
  use WaterFluxBulkType         , only : waterfluxbulk_type
  use WaterStateBulkType        , only : waterstatebulk_type
  use WaterDiagnosticBulkType        , only : waterdiagnosticbulk_type
  use Wateratm2lndBulkType        , only : wateratm2lndbulk_type
  use HumanIndexMod         , only : humanindex_type
  use ch4Mod                , only : ch4_type
  use PhotosynthesisMod     , only : photosyns_type
  use GridcellType          , only : grc
  use ColumnType            , only : col
  use PatchType             , only : patch
  use EDTypesMod            , only : ed_site_type
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type

  use EDPftvarcon         , only : EDPftvarcon_inst

  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: MossLichenReadNML     ! Read in namelist settings for moss lichen
  public :: MossLichenPhotosynthesis      ! Calculate photosynthestic rate for moss lichen
  public :: MossLichenAlbedo    ! Calculate albedo and radiation transfer for moss lichen
  public :: MossLichenHydro     ! Calculate water interception, storage and evaporation
  public :: MossLichenHeat      ! Calculate heat flux, storage for moss and lichen
  public :: wrap_MossLichen    ! wrap up input and output variable of moss&lichen between clm and fates

  !------------------------------------------------------------------------------

contains

  subroutine MossLichenReadNML(NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the namelist for Canopy Fluxes
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'MossLichenReadNML'
    character(len=*), parameter :: nmlname = 'canopyfluxes_inparm' ! ???
    !-----------------------------------------------------------------------

    namelist /canopyfluxes_inparm/ use_undercanopy_stability
    namelist /canopyfluxes_inparm/ itmax_canopy_fluxes

    ! Initialize options to default values, in case they are not specified in
    ! the namelist

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=canopyfluxes_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if

       if (itmax_canopy_fluxes < 1) then
          call endrun(msg=' ERROR: expecting itmax_canopy_fluxes > 0 ' // &
            errMsg(sourcefile, __LINE__))
       end if

       call relavu( unitn )
    end if

    call shr_mpi_bcast (use_undercanopy_stability, mpicom)
    call shr_mpi_bcast (itmax_canopy_fluxes, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=canopyfluxes_inparm)
       write(iulog,*) ' '
    end if

  end subroutine MossLichenReadNML

  subroutine MossLichenPhotosynthesis(bounds,  num_exposedvegp, filter_exposedvegp,                  &
       clm_fates, nc, active_layer_inst, atm2lnd_inst, canopystate_inst,                 &
       energyflux_inst, frictionvel_inst, soilstate_inst, solarabs_inst, surfalb_inst,   &
       temperature_inst, waterfluxbulk_inst, waterstatebulk_inst,                        &
       waterdiagnosticbulk_inst, wateratm2lndbulk_inst, ch4_inst, ozone_inst,            &
       photosyns_inst, &
       humanindex_inst, soil_water_retention_curve, &
       downreg_patch, leafn_patch, froot_carbon, croot_carbon)
  end subroutine MossLichenPhotosynthesis

  subroutine MossLichenAlbedo(bounds,  num_exposedvegp, filter_exposedvegp,                  &
       clm_fates, nc, active_layer_inst, atm2lnd_inst, canopystate_inst,                 &
       energyflux_inst, frictionvel_inst, soilstate_inst, solarabs_inst, surfalb_inst,   &
       temperature_inst, waterfluxbulk_inst, waterstatebulk_inst,                        &
       waterdiagnosticbulk_inst, wateratm2lndbulk_inst, ch4_inst, ozone_inst,            &
       photosyns_inst, &
       humanindex_inst, soil_water_retention_curve, &
       downreg_patch, leafn_patch, froot_carbon, croot_carbon)
  end subroutine MossLichenAlbedo

  subroutine MossLichenHydro(bounds,  num_exposedvegp, filter_exposedvegp,                  &
       clm_fates, nc, active_layer_inst, atm2lnd_inst, canopystate_inst,                 &
       energyflux_inst, frictionvel_inst, soilstate_inst, solarabs_inst, surfalb_inst,   &
       temperature_inst, waterfluxbulk_inst, waterstatebulk_inst,                        &
       waterdiagnosticbulk_inst, wateratm2lndbulk_inst, ch4_inst, ozone_inst,            &
       photosyns_inst, &
       humanindex_inst, soil_water_retention_curve, &
       downreg_patch, leafn_patch, froot_carbon, croot_carbon)
  end subroutine MossLichenHydro

  subroutine MossLichenHeat(bounds,  num_exposedvegp, filter_exposedvegp,                  &
       clm_fates, nc, active_layer_inst, atm2lnd_inst, canopystate_inst,                 &
       energyflux_inst, frictionvel_inst, soilstate_inst, solarabs_inst, surfalb_inst,   &
       temperature_inst, waterfluxbulk_inst, waterstatebulk_inst,                        &
       waterdiagnosticbulk_inst, wateratm2lndbulk_inst, ch4_inst, ozone_inst,            &
       photosyns_inst, &
       humanindex_inst, soil_water_retention_curve, &
       downreg_patch, leafn_patch, froot_carbon, croot_carbon)
  end subroutine MossLichenHeat

  subroutine wrap_MossLichen(bounds,  num_exposedvegp, filter_exposedvegp,                  &
       clm_fates, nc, active_layer_inst, atm2lnd_inst, canopystate_inst,                 &
       energyflux_inst, frictionvel_inst, soilstate_inst, solarabs_inst, surfalb_inst,   &
       temperature_inst, waterfluxbulk_inst, waterstatebulk_inst,                        &
       waterdiagnosticbulk_inst, wateratm2lndbulk_inst, ch4_inst, ozone_inst,            &
       photosyns_inst, &
       humanindex_inst, soil_water_retention_curve, &
       downreg_patch, leafn_patch, froot_carbon, croot_carbon)
  end subroutine wrap_MossLichen
