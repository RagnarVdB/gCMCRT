module CIA_tables_interp
  use optools_data_mod
  use optools_aux, only : locate, linear_interp, bilinear_interp, Bezier_interp
  use ieee_arithmetic
  implicit none

  public :: interp_CIA_tables


contains

  subroutine interp_CIA_tables(z,CIA_work)
    use CIA_tables_Hminus
    use CIA_tables_Heminus
    use CIA_tables_H2minus
    use CIA_tables_fake_H2O_special
    implicit none

    integer, intent(in) :: z
    real(kind=dp), dimension(1), intent(out) :: CIA_work
    integer :: s, sn, j
    integer :: iwn, iwn1, iT, iT1
    real(kind=dp) :: xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval, lTGz
    real(kind=dp) :: CIA_spec
    integer l
    
    do l = 1, nwl
      CIA_work(l) = 0.0_dp
    end do
    do s = 1, nCIA
      
      ! Special species check
      if (CIA_tab(s)%form /= 4) then
        select case (CIA_tab(s)%sp)
          
        case ('H-')
          do l = 1, nwl
            call CIA_Hminus(s,l,z,CIA_spec)
            CIA_work(l) = CIA_work(l) + CIA_spec
          end do
        case ('He-')
          
          if (CIA_tab(s)%form == 2) then
            do l = 1, nwl
              call CIA_Heminus_Bell(s,l,z,CIA_spec)
              CIA_work(l) = CIA_work(l) + CIA_spec
            end do
          else
            do l = 1, nwl
              call CIA_Heminus(s,l,z,CIA_spec)
              CIA_work(l) = CIA_work(l) + CIA_spec
            end do
          end if

        case ('H2-')
          do l = 1, nwl
            call CIA_H2minus_Bell(s,l,z,CIA_spec)
            CIA_work(l) = CIA_work(l) + CIA_spec
          end do
        case ('H2O')
          do l = 1, nwl
            call Fake_H2O_special(s,l,z,CIA_spec)
            CIA_work(l) = CIA_work(l) + CIA_spec
          end do
        case default
          
          print*, 'ERROR - CIA species not found in CIA_special - STOPPING'
          print*, 'Species: ', CIA_tab(s)%sp, CIA_tab(s)%form
          stop
          
        end select
        
        ! Check for NaN's from interpolation
        do l = 1, nwl
          if (ieee_is_nan(CIA_spec) .eqv. .True.) then
            print*, 'CIA: NaN in CIA table special: ', s, l, z, CIA_tab(s)%sp
            print*, '--', CIA_spec
          end if
        end do
        cycle
      end if
        
      ! sn = 0
      ! if (CIA_tab(s)%nset > 1) then
      !   do j = 1, CIA_tab(s)%nset
      !     if (wn(l) > CIA_tab(s)%wn_s(j) .and. wn(l) < CIA_tab(s)%wn_e(j)) then
      !       if (TG_lay(z) > CIA_tab(s)%Tmin(j) .and. TG_lay(z) < CIA_tab(s)%Tmax(j)) then
      !         sn = j
      !         exit
      !       end if
      !     end if
      !   end do
      
      !   if (sn == 0) then
      
      !     do j = 1, CIA_tab(s)%nset
      !       if (wn(l) > CIA_tab(s)%wn_s(j) .and. wn(l) < CIA_tab(s)%wn_e(j)) then
      !         if (TG_lay(z) < CIA_tab(s)%Tmin(j)) then
      !           sn = j
      !           exit
      !         end if
      !         if (TG_lay(z) > CIA_tab(s)%Tmax(j)) then
      !           sn = j
      !           exit
      !         end if
      !       end if
      !     end do
      
      !     if (sn == 0) then
      !       cycle
      !     end if
      
      !   end if
      
      ! else
      ! sn = 1
      ! end if
      sn = 1
      
      ! Locate required T indexes in CIA wn array for layer temperature
      call locate(CIA_tab(s)%T(sn,1:CIA_tab(s)%nT(sn)),TG_lay(z),iT)
      iT1 = iT + 1
      lTGz = log10(TG_lay(z))

      do l = 1, nwl
        ! Prelocated wn indexes in CIA wn array
        iwn = iwns(s, l)
        iwn1 = iwn + 1

        ! iwn = 1
        ! do while (CIA_tab(s)%wn(sn,iwn) < wn(l) .and. iwn < CIA_tab(s)%irec(sn))
        !   iwn = iwn + 1
        ! end do


        ! Check in wavenumber within bounds
        if ((iwn1 > CIA_tab(s)%irec(sn)) .or. (iwn < 1)) then
          cycle
        end if

        !! Perform temperature edge case check
        if (iT < 1) then

          ! Temperature of layer is outside lower bounds of table
          !print*, 'CIA: TG_lay < minval(T) @: ', CIA_tab(s)%sp, z, TG_lay(z), minval(CIA_tab(s)%T(:)), 'Assuming = minval(T)'

          ! Perform wn linear interp to minval(T)
          xval = lwn(l) ; x0 = CIA_tab(s)%lwn(sn,iwn) ; x1 = CIA_tab(s)%lwn(sn,iwn1)
          y0 = CIA_tab(s)%ltab(sn,iwn,1) ; y1 = CIA_tab(s)%ltab(sn,iwn1,1)

          ! Perform log linear interpolation
          call linear_interp(xval, x0, x1, y0, y1, yval)
          yval = 10.0_dp**(yval)

          ! Check for NaN's from interpolation
          if (ieee_is_nan(yval) .eqv. .True.) then
            print*, 'CIA: NaN in CIA table linear_log_interp: ', l, z, CIA_tab(s)%sp
            print*, '--', xval, yval, x0, x1, y0, y1
          end if

          ! Add to result to work variable in units of [cm-1]
          CIA_work(l) = CIA_work(l) + yval &
            & * VMR_lay(CIA_tab(s)%iVMR(1),z) * N_lay(z) &
            & * VMR_lay(CIA_tab(s)%iVMR(2),z) * N_lay(z)

        else if (iT1 > CIA_tab(s)%nT(sn)) then

          ! Temperature of layer is outside upper bounds of table
          !print*, 'CIA: TG_lay > maxval(T) @: ', CIA_tab(s)%sp, z, TG_lay(z), maxval(CIA_tab(s)%T(:)), 'Assuming = maxval(T)'

          ! Perform wn linear interp to maxval(T)
          xval = lwn(l) ; x0 = CIA_tab(s)%lwn(sn,iwn) ; x1 = CIA_tab(s)%lwn(sn,iwn1)
          y0 = CIA_tab(s)%ltab(sn,iwn,CIA_tab(s)%nT(1)) ; y1 = CIA_tab(s)%ltab(sn,iwn1,CIA_tab(s)%nT(1))

          ! Perform log linear interpolation
          call linear_interp(xval, x0, x1, y0, y1, yval)
          yval = 10.0_dp**(yval)

          ! Check for NaN's from interpolation
          if (ieee_is_nan(yval) .eqv. .True.) then
            print*, 'CIA: NaN in CIA table linear_log_interp: ', l, z, CIA_tab(s)%sp
            print*, '--', xval, yval, x0, x1, y0, y1
          end if

          ! Add to result to work variable in units of [cm-1]
          CIA_work(l) = CIA_work(l) + yval &
            & * VMR_lay(CIA_tab(s)%iVMR(1),z) * N_lay(z) &
            & * VMR_lay(CIA_tab(s)%iVMR(2),z) * N_lay(z)

        else

          !! wn and T are within the table bounds
          ! xval = wn(l) ; x0 = CIA_tab(s)%wn(sn,iwn) ; x1 = CIA_tab(s)%wn(sn,iwn1)
          ! yval = TG_lay(z) ; y0 = CIA_tab(s)%T(sn,iT) ; y1 = CIA_tab(s)%T(sn,iT1)
          xval = lwn(l) ; x0 = CIA_tab(s)%lwn(sn,iwn) ; x1 = CIA_tab(s)%lwn(sn,iwn1)
          yval = lTGz ; y0 = CIA_tab(s)%lT(sn,iT) ; y1 = CIA_tab(s)%lT(sn,iT1)
          ! a00 = CIA_tab(s)%tab(sn,iwn,iT) ; a10 = CIA_tab(s)%tab(sn,iwn1,iT)
          ! a01 = CIA_tab(s)%tab(sn,iwn,iT1) ; a11 = CIA_tab(s)%tab(sn,iwn1,iT1)
          a00 = CIA_tab(s)%ltab(sn,iwn,iT) ; a10 = CIA_tab(s)%ltab(sn,iwn1,iT)
          a01 = CIA_tab(s)%ltab(sn,iwn,iT1) ; a11 = CIA_tab(s)%ltab(sn,iwn1,iT1)

          ! Perform bi-linear interpolation
          ! call bilinear_log_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)
          call bilinear_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)
          aval = 10.0_dp**(aval)
          ! Check for NaN's from bi-linear interpolation
          if (ieee_is_nan(aval) .eqv. .True.) then
            print*, 'CIA: NaN in CIA table bilinear_log_interp: ', l, z, CIA_tab(s)%sp
            print*, '--', xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval
          end if

          ! Add to result to work variable in units of [cm-1]
          CIA_work(l) = CIA_work(l) + aval &
            & * VMR_lay(CIA_tab(s)%iVMR(1),z) * N_lay(z) &
            & * VMR_lay(CIA_tab(s)%iVMR(2),z) * N_lay(z)

        end if

      end do
    end do

  end subroutine interp_CIA_tables

end module CIA_tables_interp
