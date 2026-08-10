!
! Copyright (C) 2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution.
!
!----------------------------------------------------------------------------
MODULE caimd_module
  !----------------------------------------------------------------------------
  !! PW-local constrained AIMD support for small strain-and-hold cell driving.
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : e2, au_ps, k_boltzmann_ry, ry_kbar
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  LOGICAL  :: caimd_enabled_ = .FALSE.
  CHARACTER(LEN=32) :: caimd_loading_mode_ = 'legacy_strain'
  INTEGER  :: caimd_stride_ = 1
  REAL(DP) :: caimd_step_strain_(6) = 0.0_DP
  REAL(DP) :: caimd_step_defgrad_(3,3) = 0.0_DP
  LOGICAL  :: caimd_mask_(6) = .FALSE.
  REAL(DP) :: caimd_cumulative_strain_(6) = 0.0_DP
  REAL(DP) :: caimd_reference_cell_(3,3) = 0.0_DP
  LOGICAL  :: caimd_reference_set_ = .FALSE.
  INTEGER  :: caimd_event_count_ = 0
  LOGICAL  :: caimd_last_event_ = .FALSE.
  REAL(DP) :: caimd_last_increment_(6) = 0.0_DP
  LOGICAL  :: caimd_output_initialized_ = .FALSE.
  LOGICAL  :: caimd_defgrad_output_initialized_ = .FALSE.
  LOGICAL  :: caimd_have_temperature_ = .FALSE.
  LOGICAL  :: caimd_cumulative_warning_emitted_ = .FALSE.
  REAL(DP) :: caimd_temperature_ = 0.0_DP
  REAL(DP) :: caimd_last_velocity_jump_max_ = 0.0_DP
  REAL(DP) :: caimd_last_velocity_jump_rms_ = 0.0_DP
  REAL(DP) :: caimd_last_f_before_(3,3) = 0.0_DP
  REAL(DP) :: caimd_last_f_after_(3,3) = 0.0_DP
  REAL(DP) :: caimd_last_a_event_(3,3) = 0.0_DP
  REAL(DP) :: caimd_last_green_lagrange_(3,3) = 0.0_DP
  REAL(DP) :: caimd_last_det_f_before_ = 1.0_DP
  REAL(DP) :: caimd_last_det_f_after_ = 1.0_DP
  REAL(DP) :: caimd_last_det_a_event_ = 1.0_DP
  INTEGER  :: caimd_last_event_index_ = 0
  !
  REAL(DP), PARAMETER :: caimd_max_event_strain = 1.0D-2
  REAL(DP), PARAMETER :: caimd_max_condition = 1.0D12
  REAL(DP), PARAMETER :: caimd_det_floor = 1.0D-18
  REAL(DP), PARAMETER :: caimd_change_tol = 1.0D-14
  REAL(DP), PARAMETER :: caimd_zero_tol = 0.0_DP
  !
  PUBLIC :: caimd_configure, caimd_enabled, caimd_controls_cell
  PUBLIC :: caimd_project_tangent, caimd_apply_final_projection
  PUBLIC :: caimd_write_premove, caimd_note_ionic_temperature
  PUBLIC :: caimd_step_strain, caimd_step_defgrad, caimd_mask, caimd_stride
  PUBLIC :: caimd_loading_mode
  PUBLIC :: caimd_event_count, caimd_cumulative_strain
  !
CONTAINS
  !
  SUBROUTINE caimd_configure( enabled, loading_mode, stride, step_strain,     &
                              step_defgrad, mask, calculation, calc,         &
                              cell_dofree, restart, fix_volume, fix_area )
    !
    LOGICAL, INTENT(IN) :: enabled, restart, fix_volume, fix_area
    INTEGER, INTENT(IN) :: stride
    REAL(DP), INTENT(IN) :: step_strain(6)
    REAL(DP), INTENT(IN) :: step_defgrad(3,3)
    LOGICAL, INTENT(IN) :: mask(6)
    CHARACTER(LEN=*), INTENT(IN) :: loading_mode, calculation, calc, cell_dofree
    INTEGER :: i, j
    REAL(DP) :: event_inc(6)
    LOGICAL :: strain_nonzero, defgrad_nonzero
    !
    caimd_enabled_ = enabled
    SELECT CASE ( TRIM( loading_mode ) )
    CASE ( 'legacy_strain' )
       caimd_loading_mode_ = 'legacy_strain'
    CASE ( 'linear_defgrad' )
       caimd_loading_mode_ = 'linear_defgrad'
    CASE DEFAULT
       CALL errore( 'caimd_configure', &
            'unknown caimd_loading_mode: ' // TRIM( loading_mode ), 1 )
    END SELECT
    caimd_stride_ = stride
    caimd_step_strain_ = step_strain
    caimd_step_defgrad_ = step_defgrad
    caimd_mask_ = mask
    caimd_cumulative_strain_ = 0.0_DP
    caimd_reference_cell_ = 0.0_DP
    caimd_reference_set_ = .FALSE.
    caimd_event_count_ = 0
    caimd_last_event_ = .FALSE.
    caimd_last_increment_ = 0.0_DP
    caimd_output_initialized_ = .FALSE.
    caimd_defgrad_output_initialized_ = .FALSE.
    caimd_have_temperature_ = .FALSE.
    caimd_cumulative_warning_emitted_ = .FALSE.
    caimd_temperature_ = 0.0_DP
    caimd_last_velocity_jump_max_ = 0.0_DP
    caimd_last_velocity_jump_rms_ = 0.0_DP
    caimd_last_f_before_ = 0.0_DP
    caimd_last_f_after_ = 0.0_DP
    caimd_last_a_event_ = 0.0_DP
    caimd_last_green_lagrange_ = 0.0_DP
    DO i = 1, 3
       caimd_last_f_before_(i,i) = 1.0_DP
       caimd_last_f_after_(i,i) = 1.0_DP
       caimd_last_a_event_(i,i) = 1.0_DP
    END DO
    caimd_last_det_f_before_ = 1.0_DP
    caimd_last_det_f_after_ = 1.0_DP
    caimd_last_det_a_event_ = 1.0_DP
    caimd_last_event_index_ = 0
    !
    IF ( stride <= 0 ) CALL errore( 'caimd_configure', &
         'caimd_stride must be positive', 1 )
    !
    DO i = 1, 6
       IF ( .NOT. caimd_finite( step_strain(i) ) ) &
          CALL errore( 'caimd_configure', &
               'caimd_step_strain contains a non-finite value', i )
    END DO
    event_inc = REAL( stride, DP ) * step_strain
    DO i = 1, 6
       IF ( .NOT. caimd_finite( event_inc(i) ) ) &
          CALL errore( 'caimd_configure', &
               'caimd event strain increment is non-finite', i )
    END DO
    DO j = 1, 3
       DO i = 1, 3
          IF ( .NOT. caimd_finite( step_defgrad(i,j) ) ) &
             CALL errore( 'caimd_configure', &
                  'caimd_step_defgrad contains a non-finite value', 3 * (j - 1) + i )
       END DO
    END DO
    strain_nonzero = ANY( ABS( step_strain ) > caimd_zero_tol )
    defgrad_nonzero = ANY( ABS( step_defgrad ) > caimd_zero_tol )
    IF ( TRIM( caimd_loading_mode_ ) == 'legacy_strain' .AND. defgrad_nonzero ) &
       CALL errore( 'caimd_configure', &
            'caimd_step_defgrad requires caimd_loading_mode=''linear_defgrad''', 1 )
    IF ( TRIM( caimd_loading_mode_ ) == 'linear_defgrad' .AND. strain_nonzero ) &
       CALL errore( 'caimd_configure', &
            'caimd_step_strain requires caimd_loading_mode=''legacy_strain''', 1 )
    IF ( MAXVAL( ABS( event_inc ) ) > caimd_max_event_strain ) &
       CALL errore( 'caimd_configure', &
            'max abs(caimd_stride*caimd_step_strain) exceeds 1.0d-2', 1 )
    !
    IF ( .NOT. enabled ) RETURN
    !
    IF ( TRIM( calculation ) /= 'vc-md' .OR. TRIM( calc ) /= 'cd' ) &
       CALL errore( 'caimd_configure', &
            'caimd=.true. requires calculation=''vc-md'' and cell_dynamics=''pr''', 1 )
    IF ( restart ) CALL errore( 'caimd_configure', &
         'caimd legacy_strain and linear_defgrad modes support restart_mode=''from_scratch'' only', 1 )
    IF ( TRIM( cell_dofree ) /= 'all' ) CALL errore( 'caimd_configure', &
         'caimd version 1 requires cell_dofree=''all''', 1 )
    IF ( fix_volume .OR. fix_area ) CALL errore( 'caimd_configure', &
         'caimd version 1 does not support fixed-volume or fixed-area cell modes', 1 )
    !
    CALL caimd_print_configuration()
    !
  END SUBROUTINE caimd_configure
  !
  LOGICAL FUNCTION caimd_enabled()
    caimd_enabled = caimd_enabled_
  END FUNCTION caimd_enabled
  !
  LOGICAL FUNCTION caimd_controls_cell()
    caimd_controls_cell = caimd_enabled_ .AND. &
       ( ( TRIM( caimd_loading_mode_ ) == 'legacy_strain' .AND. &
           ANY( caimd_mask_ ) ) .OR. &
         TRIM( caimd_loading_mode_ ) == 'linear_defgrad' )
  END FUNCTION caimd_controls_cell
  !
  FUNCTION caimd_loading_mode() RESULT( loading_mode )
    CHARACTER(LEN=32) :: loading_mode
    loading_mode = caimd_loading_mode_
  END FUNCTION caimd_loading_mode
  !
  INTEGER FUNCTION caimd_stride()
    caimd_stride = caimd_stride_
  END FUNCTION caimd_stride
  !
  FUNCTION caimd_step_strain() RESULT( step_strain )
    REAL(DP) :: step_strain(6)
    step_strain = caimd_step_strain_
  END FUNCTION caimd_step_strain
  !
  FUNCTION caimd_step_defgrad() RESULT( step_defgrad )
    REAL(DP) :: step_defgrad(3,3)
    step_defgrad = caimd_step_defgrad_
  END FUNCTION caimd_step_defgrad
  !
  FUNCTION caimd_mask() RESULT( mask )
    LOGICAL :: mask(6)
    mask = caimd_mask_
  END FUNCTION caimd_mask
  !
  INTEGER FUNCTION caimd_event_count()
    caimd_event_count = caimd_event_count_
  END FUNCTION caimd_event_count
  !
  FUNCTION caimd_cumulative_strain() RESULT( strain )
    REAL(DP) :: strain(6)
    strain = caimd_cumulative_strain_
  END FUNCTION caimd_cumulative_strain
  !
  SUBROUTINE caimd_project_tangent( h, x, label )
    !
    REAL(DP), INTENT(IN) :: h(3,3)
    REAL(DP), INTENT(INOUT) :: x(3,3)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(DP) :: hinv(3,3), k(3,3), ksym(3,3), vec(6), det
    INTEGER :: i
    !
    IF ( .NOT. caimd_controls_cell() ) RETURN
    IF ( TRIM( caimd_loading_mode_ ) == 'linear_defgrad' ) RETURN
    !
    CALL caimd_inverse3( h, hinv, det, TRIM( label ) // ' inverse cell' )
    k = MATMUL( x, hinv )
    ksym = 0.5_DP * ( k + TRANSPOSE( k ) )
    CALL caimd_mat_to_vec( ksym, vec )
    DO i = 1, 6
       IF ( caimd_mask_(i) ) vec(i) = 0.0_DP
    END DO
    CALL caimd_vec_to_mat( vec, ksym )
    x = MATMUL( ksym, h )
    !
  END SUBROUTINE caimd_project_tangent
  !
  SUBROUTINE caimd_apply_final_projection( ionic_step, nat, h_old, h_candidate, &
                                           ratd, rat2d, rat2di, avecd, avec2d, &
                                           avec2di )
    !
    INTEGER, INTENT(IN) :: ionic_step, nat
    REAL(DP), INTENT(IN) :: h_old(3,3)
    REAL(DP), INTENT(INOUT) :: h_candidate(3,3)
    REAL(DP), INTENT(INOUT) :: ratd(3,nat), rat2d(3,nat), rat2di(3,nat)
    REAL(DP), INTENT(INOUT) :: avecd(3,3), avec2d(3,3), avec2di(3,3)
    REAL(DP) :: h_pr(3,3), h_new(3,3), inc(6)
    LOGICAL :: event_step, changed
    !
    IF ( .NOT. caimd_controls_cell() ) RETURN
    !
    h_pr = h_candidate
    caimd_last_velocity_jump_max_ = 0.0_DP
    caimd_last_velocity_jump_rms_ = 0.0_DP
    CALL caimd_project_cell_increment( ionic_step, h_old, h_pr, h_new, &
                                       event_step, inc )
    changed = MAXVAL( ABS( h_new - h_pr ) ) > caimd_change_tol * &
              MAX( 1.0_DP, caimd_norm_inf( h_pr ) )
    IF ( changed ) CALL caimd_transform_atomic_history( nat, h_pr, h_new, &
                                                        ratd, rat2d, rat2di )
    IF ( event_step .AND. TRIM( caimd_loading_mode_ ) == 'linear_defgrad' ) &
       CALL caimd_transform_cell_history( caimd_last_a_event_, avecd, avec2d, &
                                          avec2di )
    h_candidate = h_new
    CALL caimd_project_tangent( h_candidate, avecd, 'caimd final avecd' )
    CALL caimd_project_tangent( h_candidate, avec2d, 'caimd final avec2d' )
    CALL caimd_project_tangent( h_candidate, avec2di, 'caimd final avec2di' )
    !
    caimd_last_event_ = event_step
    caimd_last_increment_ = 0.0_DP
    IF ( event_step ) THEN
       caimd_last_increment_ = inc
       caimd_cumulative_strain_ = caimd_cumulative_strain_ + inc
       caimd_event_count_ = caimd_event_count_ + 1
       CALL caimd_warn_large_cumulative_strain()
       IF ( TRIM( caimd_loading_mode_ ) == 'linear_defgrad' ) &
          CALL caimd_write_defgrad_event( ionic_step )
    END IF
    !
  END SUBROUTINE caimd_apply_final_projection
  !
  SUBROUTINE caimd_note_ionic_temperature( nat, ionic_ekin )
    !
    INTEGER, INTENT(IN) :: nat
    REAL(DP), INTENT(IN) :: ionic_ekin
    INTEGER :: ndof
    !
    IF ( .NOT. caimd_enabled_ ) RETURN
    ndof = MAX( 3 * ( nat - 1 ), 1 )
    caimd_temperature_ = 2.0_DP * ionic_ekin / &
                         ( REAL( ndof, DP ) * k_boltzmann_ry )
    caimd_have_temperature_ = .TRUE.
    !
  END SUBROUTINE caimd_note_ionic_temperature
  !
  SUBROUTINE caimd_write_premove( ionic_step, dt, h, volume, stress, input_temp )
    !
    USE io_files,  ONLY : tmp_dir, prefix, delete_if_present
    USE io_global, ONLY : ionode
    !
    INTEGER, INTENT(IN) :: ionic_step
    REAL(DP), INTENT(IN) :: dt, h(3,3), volume, stress(3,3), input_temp
    INTEGER :: iun
    CHARACTER(LEN=512) :: filename
    LOGICAL :: event_next
    REAL(DP) :: actual(6), stress_gpa(6), temp_out, time_ps
    !
    IF ( .NOT. caimd_enabled_ ) RETURN
    !
    CALL caimd_ensure_reference( h )
    CALL caimd_actual_strain( h, actual )
    CALL caimd_mat_to_vec( stress, stress_gpa )
    stress_gpa = stress_gpa * ry_kbar * 0.1_DP
    event_next = caimd_controls_cell() .AND. caimd_is_event_step( ionic_step )
    time_ps = REAL( ionic_step - 1, DP ) * dt * e2 * au_ps
    IF ( caimd_have_temperature_ ) THEN
       temp_out = caimd_temperature_
    ELSE
       temp_out = input_temp
    END IF
    !
    IF ( ionode ) THEN
       filename = TRIM( tmp_dir ) // TRIM( prefix ) // '.caimd'
       IF ( .NOT. caimd_output_initialized_ ) THEN
          CALL delete_if_present( TRIM( filename ) )
          OPEN( NEWUNIT = iun, FILE = TRIM( filename ), STATUS = 'NEW', &
                FORM = 'FORMATTED' )
          WRITE( iun, '(A)' ) '# CAIMD_V1 stress-strain output'
          WRITE( iun, '(A,A)' ) '# loading_mode: ', TRIM( caimd_loading_mode_ )
          WRITE( iun, '(A)' ) '# units: time ps; strain dimensionless; stress GPa; volume bohr^3; temperature K'
          WRITE( iun, '(A)' ) '# shear convention: tensor components yz xz xy; engineering shear is 2*tensor shear'
          IF ( TRIM( caimd_loading_mode_ ) == 'linear_defgrad' ) &
             WRITE( iun, '(A)' ) '# linear_defgrad event matrices are written to prefix.caimd_defgrad'
          WRITE( iun, '(A)' ) '# columns: ionic_step time_ps imposed_events_done following_move_is_event ' // &
                              'nominal_xx nominal_yy nominal_zz nominal_yz nominal_xz nominal_xy ' // &
                              'actual_xx actual_yy actual_zz actual_yz actual_xz actual_xy ' // &
                              'stress_xx_GPa stress_yy_GPa stress_zz_GPa stress_yz_GPa stress_xz_GPa stress_xy_GPa ' // &
                              'temperature_K volume_bohr3 velocity_jump_max velocity_jump_rms'
          caimd_output_initialized_ = .TRUE.
       ELSE
          OPEN( NEWUNIT = iun, FILE = TRIM( filename ), STATUS = 'OLD', &
                POSITION = 'APPEND', FORM = 'FORMATTED' )
       END IF
       WRITE( iun, '(I8,1X,ES24.16,1X,I8,1X,L1,22(1X,ES24.16))' ) &
          ionic_step, time_ps, caimd_event_count_, event_next, &
          caimd_cumulative_strain_, actual, stress_gpa, temp_out, volume, &
          caimd_last_velocity_jump_max_, caimd_last_velocity_jump_rms_
       CLOSE( UNIT = iun, STATUS = 'KEEP' )
    END IF
    !
  END SUBROUTINE caimd_write_premove
  !
  SUBROUTINE caimd_project_cell_increment( ionic_step, h_old, h_pr, h_new, &
                                           event_step, inc )
    !
    INTEGER, INTENT(IN) :: ionic_step
    REAL(DP), INTENT(IN) :: h_old(3,3), h_pr(3,3)
    REAL(DP), INTENT(OUT) :: h_new(3,3), inc(6)
    LOGICAL, INTENT(OUT) :: event_step
    REAL(DP) :: href_inv(3,3), fpr(3,3), lpr(3,3), epr(3,3), enew(3,3)
    REAL(DP) :: vec(6), target(6), det, f_before_inv(3,3)
    INTEGER :: i
    !
    CALL caimd_ensure_reference( h_old )
    event_step = caimd_is_event_step( ionic_step )
    inc = 0.0_DP
    IF ( TRIM( caimd_loading_mode_ ) == 'linear_defgrad' ) THEN
       h_new = h_pr
       IF ( event_step ) THEN
          caimd_last_event_index_ = caimd_event_count_ + 1
          CALL caimd_build_f_load( caimd_event_count_, caimd_last_f_before_ )
          CALL caimd_build_f_load( caimd_last_event_index_, caimd_last_f_after_ )
          CALL caimd_inverse3( caimd_last_f_before_, f_before_inv, &
                               caimd_last_det_f_before_, &
                               'caimd linear_defgrad F_load before event' )
          caimd_last_det_f_after_ = caimd_det3( caimd_last_f_after_ )
          IF ( .NOT. caimd_finite( caimd_last_det_f_before_ ) .OR. &
               caimd_last_det_f_before_ <= caimd_det_floor ) &
             CALL errore( 'caimd_project_cell_increment', &
                  'caimd linear_defgrad F_load before event has non-positive determinant', 1 )
          IF ( .NOT. caimd_finite( caimd_last_det_f_after_ ) .OR. &
               caimd_last_det_f_after_ <= caimd_det_floor ) &
             CALL errore( 'caimd_project_cell_increment', &
                  'caimd linear_defgrad F_load after event has non-positive determinant', 1 )
          caimd_last_a_event_ = MATMUL( caimd_last_f_after_, f_before_inv )
          caimd_last_det_a_event_ = caimd_det3( caimd_last_a_event_ )
          IF ( .NOT. caimd_finite( caimd_last_det_a_event_ ) .OR. &
               caimd_last_det_a_event_ <= caimd_det_floor ) &
             CALL errore( 'caimd_project_cell_increment', &
                  'caimd linear_defgrad A_event has non-positive determinant', 1 )
          CALL caimd_green_lagrange( caimd_last_f_after_, &
                                     caimd_last_green_lagrange_ )
          ! QE stores lattice vectors as columns, so the physical affine
          ! relation r_new=A*r_old is H_new=A*H_old for the cell matrix.
          h_new = MATMUL( caimd_last_a_event_, h_pr )
       END IF
       CALL caimd_check_cell( h_new, 'caimd linear_defgrad projected cell' )
       RETURN
    END IF
    !
    CALL caimd_inverse3( caimd_reference_cell_, href_inv, det, &
                         'caimd reference cell inverse' )
    fpr = MATMUL( h_pr, href_inv )
    lpr = fpr
    DO i = 1, 3
       lpr(i,i) = lpr(i,i) - 1.0_DP
    END DO
    epr = 0.5_DP * ( lpr + TRANSPOSE( lpr ) )
    CALL caimd_mat_to_vec( epr, vec )
    !
    IF ( event_step ) THEN
       DO i = 1, 6
          IF ( caimd_mask_(i) ) &
             inc(i) = REAL( caimd_stride_, DP ) * caimd_step_strain_(i)
       END DO
    END IF
    target = caimd_cumulative_strain_ + inc
    DO i = 1, 6
       IF ( caimd_mask_(i) ) vec(i) = target(i)
    END DO
    CALL caimd_vec_to_mat( vec, enew )
    DO i = 1, 3
       enew(i,i) = enew(i,i) + 1.0_DP
    END DO
    h_new = MATMUL( enew, caimd_reference_cell_ )
    CALL caimd_check_cell( h_new, 'caimd projected cell' )
    !
  END SUBROUTINE caimd_project_cell_increment
  !
  LOGICAL FUNCTION caimd_is_event_step( ionic_step )
    !
    INTEGER, INTENT(IN) :: ionic_step
    !
    caimd_is_event_step = caimd_controls_cell() .AND. &
                          MOD( ionic_step - 1, caimd_stride_ ) == 0
    !
  END FUNCTION caimd_is_event_step
  !
  SUBROUTINE caimd_build_f_load( event_index, f_load )
    !
    INTEGER, INTENT(IN) :: event_index
    REAL(DP), INTENT(OUT) :: f_load(3,3)
    INTEGER :: i
    !
    f_load = REAL( event_index, DP ) * REAL( caimd_stride_, DP ) * &
             caimd_step_defgrad_
    DO i = 1, 3
       f_load(i,i) = f_load(i,i) + 1.0_DP
    END DO
    !
  END SUBROUTINE caimd_build_f_load
  !
  SUBROUTINE caimd_green_lagrange( f_load, e_green )
    !
    REAL(DP), INTENT(IN) :: f_load(3,3)
    REAL(DP), INTENT(OUT) :: e_green(3,3)
    INTEGER :: i
    !
    e_green = 0.5_DP * MATMUL( TRANSPOSE( f_load ), f_load )
    DO i = 1, 3
       e_green(i,i) = e_green(i,i) - 0.5_DP
    END DO
    !
  END SUBROUTINE caimd_green_lagrange
  !
  SUBROUTINE caimd_transform_cell_history( a_event, avecd, avec2d, avec2di )
    !
    REAL(DP), INTENT(IN) :: a_event(3,3)
    REAL(DP), INTENT(INOUT) :: avecd(3,3), avec2d(3,3), avec2di(3,3)
    !
    avecd = MATMUL( a_event, avecd )
    avec2d = MATMUL( a_event, avec2d )
    avec2di = MATMUL( a_event, avec2di )
    !
  END SUBROUTINE caimd_transform_cell_history
  !
  SUBROUTINE caimd_write_defgrad_event( ionic_step )
    !
    USE io_files,  ONLY : tmp_dir, prefix, delete_if_present
    USE io_global, ONLY : ionode
    !
    INTEGER, INTENT(IN) :: ionic_step
    INTEGER :: i, j, iun
    CHARACTER(LEN=512) :: filename
    !
    IF ( TRIM( caimd_loading_mode_ ) /= 'linear_defgrad' ) RETURN
    IF ( .NOT. ionode ) RETURN
    !
    filename = TRIM( tmp_dir ) // TRIM( prefix ) // '.caimd_defgrad'
    IF ( .NOT. caimd_defgrad_output_initialized_ ) THEN
       CALL delete_if_present( TRIM( filename ) )
       OPEN( NEWUNIT = iun, FILE = TRIM( filename ), STATUS = 'NEW', &
             FORM = 'FORMATTED' )
       WRITE( iun, '(A)' ) '# CAIMD_LINEAR_DEFGRAD_V1 event output'
       WRITE( iun, '(A)' ) '# F_load(m)=I+m*caimd_stride*caimd_step_defgrad'
       WRITE( iun, '(A)' ) '# A_event(m)=F_load(m)*inverse(F_load(m-1))'
       WRITE( iun, '(A)' ) '# matrix order: row-major 11 12 13 21 22 23 31 32 33'
       WRITE( iun, '(A)' ) '# columns: ionic_step event_index det_F_before det_F_after det_A_event ' // &
                           'F_before_11 F_before_12 F_before_13 F_before_21 F_before_22 F_before_23 ' // &
                           'F_before_31 F_before_32 F_before_33 F_after_11 F_after_12 F_after_13 ' // &
                           'F_after_21 F_after_22 F_after_23 F_after_31 F_after_32 F_after_33 ' // &
                           'A_event_11 A_event_12 A_event_13 A_event_21 A_event_22 A_event_23 ' // &
                           'A_event_31 A_event_32 A_event_33 E_11 E_12 E_13 E_21 E_22 E_23 ' // &
                           'E_31 E_32 E_33 velocity_jump_max velocity_jump_rms'
       caimd_defgrad_output_initialized_ = .TRUE.
    ELSE
       OPEN( NEWUNIT = iun, FILE = TRIM( filename ), STATUS = 'OLD', &
             POSITION = 'APPEND', FORM = 'FORMATTED' )
    END IF
    WRITE( iun, '(I8,1X,I8,41(1X,ES24.16))' ) &
       ionic_step, caimd_last_event_index_, caimd_last_det_f_before_, &
       caimd_last_det_f_after_, caimd_last_det_a_event_, &
       ((caimd_last_f_before_(i,j), j = 1, 3), i = 1, 3), &
       ((caimd_last_f_after_(i,j), j = 1, 3), i = 1, 3), &
       ((caimd_last_a_event_(i,j), j = 1, 3), i = 1, 3), &
       ((caimd_last_green_lagrange_(i,j), j = 1, 3), i = 1, 3), &
       caimd_last_velocity_jump_max_, caimd_last_velocity_jump_rms_
    CLOSE( UNIT = iun, STATUS = 'KEEP' )
    !
  END SUBROUTINE caimd_write_defgrad_event
  !
  SUBROUTINE caimd_print_configuration()
    !
    USE io_global, ONLY : ionode, stdout
    !
    INTEGER :: i
    !
    IF ( .NOT. caimd_enabled_ ) RETURN
    IF ( .NOT. ionode ) RETURN
    WRITE( stdout, '(/,5X,A)' ) 'CAIMD loading configuration:'
    WRITE( stdout, '(8X,A,A)' ) 'caimd_loading_mode = ', TRIM( caimd_loading_mode_ )
    WRITE( stdout, '(8X,A,I8)' ) 'caimd_stride = ', caimd_stride_
    WRITE( stdout, '(8X,A)' ) 'caimd_step_defgrad matrix, per MD step:'
    DO i = 1, 3
       WRITE( stdout, '(10X,3(1X,ES14.6))' ) caimd_step_defgrad_(i,1:3)
    END DO
    WRITE( stdout, '(8X,A)' ) 'event convention: steps 1, 1+stride, 1+2*stride, ...'
    WRITE( stdout, '(8X,A)' ) 'restart support: not implemented; from_scratch only'
    !
  END SUBROUTINE caimd_print_configuration
  !
  SUBROUTINE caimd_transform_atomic_history( nat, h_from, h_to, ratd, rat2d, &
                                             rat2di )
    !
    INTEGER, INTENT(IN) :: nat
    REAL(DP), INTENT(IN) :: h_from(3,3), h_to(3,3)
    REAL(DP), INTENT(INOUT) :: ratd(3,nat), rat2d(3,nat), rat2di(3,nat)
    REAL(DP) :: h_to_inv(3,3), transform(3,3), det
    REAL(DP) :: u_before(3), u_after(3), diff(3), jump2_sum, jump_norm
    INTEGER :: na
    !
    CALL caimd_inverse3( h_to, h_to_inv, det, 'caimd new cell inverse' )
    transform = MATMUL( h_to_inv, h_from )
    caimd_last_velocity_jump_max_ = 0.0_DP
    jump2_sum = 0.0_DP
    DO na = 1, nat
       u_before = MATMUL( h_from, ratd(:,na) )
       ratd(:,na) = MATMUL( transform, ratd(:,na) )
       rat2d(:,na) = MATMUL( transform, rat2d(:,na) )
       rat2di(:,na) = MATMUL( transform, rat2di(:,na) )
       u_after = MATMUL( h_to, ratd(:,na) )
       diff = u_after - u_before
       jump_norm = SQRT( SUM( diff * diff ) )
       caimd_last_velocity_jump_max_ = MAX( caimd_last_velocity_jump_max_, &
                                            jump_norm )
       jump2_sum = jump2_sum + SUM( diff * diff )
    END DO
    caimd_last_velocity_jump_rms_ = SQRT( jump2_sum / REAL( MAX( nat, 1 ), DP ) )
    !
  END SUBROUTINE caimd_transform_atomic_history
  !
  SUBROUTINE caimd_warn_large_cumulative_strain()
    !
    USE io_global, ONLY : ionode, stdout
    !
    REAL(DP) :: max_strain
    INTEGER :: i
    !
    IF ( caimd_cumulative_warning_emitted_ ) RETURN
    IF ( .NOT. caimd_controls_cell() ) RETURN
    max_strain = 0.0_DP
    DO i = 1, 6
       IF ( caimd_mask_(i) ) max_strain = MAX( max_strain, &
                                               ABS( caimd_cumulative_strain_(i) ) )
    END DO
    IF ( max_strain > caimd_max_event_strain ) THEN
       IF ( ionode ) WRITE( stdout, '(/,5X,A,ES12.4,A)' ) &
          'WARNING: CAIMD cumulative masked nominal strain is ', max_strain, &
          '; V1 uses a small-strain convention.'
       caimd_cumulative_warning_emitted_ = .TRUE.
    END IF
    !
  END SUBROUTINE caimd_warn_large_cumulative_strain
  !
  SUBROUTINE caimd_ensure_reference( h )
    !
    REAL(DP), INTENT(IN) :: h(3,3)
    !
    IF ( caimd_reference_set_ ) RETURN
    caimd_reference_cell_ = h
    CALL caimd_check_cell( caimd_reference_cell_, 'caimd reference cell' )
    caimd_reference_set_ = .TRUE.
    !
  END SUBROUTINE caimd_ensure_reference
  !
  SUBROUTINE caimd_actual_strain( h, strain )
    !
    REAL(DP), INTENT(IN) :: h(3,3)
    REAL(DP), INTENT(OUT) :: strain(6)
    REAL(DP) :: href_inv(3,3), ftotal(3,3), ltotal(3,3), etotal(3,3), det
    INTEGER :: i
    !
    CALL caimd_ensure_reference( h )
    CALL caimd_inverse3( caimd_reference_cell_, href_inv, det, &
                         'caimd reference inverse' )
    ftotal = MATMUL( h, href_inv )
    ltotal = ftotal
    DO i = 1, 3
       ltotal(i,i) = ltotal(i,i) - 1.0_DP
    END DO
    etotal = 0.5_DP * ( ltotal + TRANSPOSE( ltotal ) )
    CALL caimd_mat_to_vec( etotal, strain )
    !
  END SUBROUTINE caimd_actual_strain
  !
  SUBROUTINE caimd_mat_to_vec( a, v )
    !
    REAL(DP), INTENT(IN) :: a(3,3)
    REAL(DP), INTENT(OUT) :: v(6)
    !
    v(1) = a(1,1)
    v(2) = a(2,2)
    v(3) = a(3,3)
    v(4) = 0.5_DP * ( a(2,3) + a(3,2) )
    v(5) = 0.5_DP * ( a(1,3) + a(3,1) )
    v(6) = 0.5_DP * ( a(1,2) + a(2,1) )
    !
  END SUBROUTINE caimd_mat_to_vec
  !
  SUBROUTINE caimd_vec_to_mat( v, a )
    !
    REAL(DP), INTENT(IN) :: v(6)
    REAL(DP), INTENT(OUT) :: a(3,3)
    !
    a = 0.0_DP
    a(1,1) = v(1)
    a(2,2) = v(2)
    a(3,3) = v(3)
    a(2,3) = v(4)
    a(3,2) = v(4)
    a(1,3) = v(5)
    a(3,1) = v(5)
    a(1,2) = v(6)
    a(2,1) = v(6)
    !
  END SUBROUTINE caimd_vec_to_mat
  !
  SUBROUTINE caimd_check_cell( h, label )
    !
    REAL(DP), INTENT(IN) :: h(3,3)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(DP) :: hinv(3,3), det, cond
    !
    CALL caimd_inverse3( h, hinv, det, label )
    IF ( det <= caimd_det_floor ) CALL errore( 'caimd_check_cell', &
         TRIM( label ) // ' has non-positive or too-small determinant', 1 )
    cond = caimd_norm_inf( h ) * caimd_norm_inf( hinv )
    IF ( cond > caimd_max_condition ) CALL errore( 'caimd_check_cell', &
         TRIM( label ) // ' is too ill-conditioned for CAIMD v1', 1 )
    !
  END SUBROUTINE caimd_check_cell
  !
  REAL(DP) FUNCTION caimd_det3( a )
    !
    REAL(DP), INTENT(IN) :: a(3,3)
    !
    caimd_det3 = a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) - &
                 a(1,2) * ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) + &
                 a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )
    !
  END FUNCTION caimd_det3
  !
  SUBROUTINE caimd_inverse3( a, ainv, det, label )
    !
    REAL(DP), INTENT(IN) :: a(3,3)
    REAL(DP), INTENT(OUT) :: ainv(3,3), det
    CHARACTER(LEN=*), INTENT(IN) :: label
    !
    det = caimd_det3( a )
    IF ( .NOT. caimd_finite( det ) .OR. ABS( det ) <= caimd_det_floor ) &
       CALL errore( 'caimd_inverse3', TRIM( label ) // ' is singular', 1 )
    !
    ainv(1,1) =  ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
    ainv(1,2) = -( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) / det
    ainv(1,3) =  ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det
    ainv(2,1) = -( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) / det
    ainv(2,2) =  ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
    ainv(2,3) = -( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) / det
    ainv(3,1) =  ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
    ainv(3,2) = -( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) / det
    ainv(3,3) =  ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det
    IF ( ANY( ABS( ainv ) > HUGE( 1.0_DP ) * 0.1_DP ) ) &
       CALL errore( 'caimd_inverse3', TRIM( label ) // ' inverse overflow', 1 )
    !
  END SUBROUTINE caimd_inverse3
  !
  REAL(DP) FUNCTION caimd_norm_inf( a )
    !
    REAL(DP), INTENT(IN) :: a(3,3)
    INTEGER :: i
    !
    caimd_norm_inf = 0.0_DP
    DO i = 1, 3
       caimd_norm_inf = MAX( caimd_norm_inf, SUM( ABS( a(i,:) ) ) )
    END DO
    !
  END FUNCTION caimd_norm_inf
  !
  LOGICAL FUNCTION caimd_finite( x )
    !
    REAL(DP), INTENT(IN) :: x
    !
    caimd_finite = ( x == x ) .AND. ABS( x ) < HUGE( x ) * 0.1_DP
    !
  END FUNCTION caimd_finite
  !
END MODULE caimd_module
