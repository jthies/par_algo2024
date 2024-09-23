! This file is part of Fortuno.
! Licensed under the BSD-2-Clause Plus Patent license.
! SPDX-License-Identifier: BSD-2-Clause-Patent

!> Test app, collecting and executing the tests
program testapp
  use fortuno_coarray, only : execute_coa_cmd_app
  use m_unit_tests, only : get_unit_tests
  implicit none

  call execute_coa_cmd_app(&
    testitems=[&
      get_unit_tests()&
    ]&
  )

end program testapp
