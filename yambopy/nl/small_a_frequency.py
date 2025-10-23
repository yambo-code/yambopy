function small_a_frequency(W_i,E_field,iErr)
  !
  use pars,      ONLY:SP,cZERO,cI,pi,schlen
  use stderr,    ONLY:STRING_split
  use fields,    ONLY:ext_field,n_fields_defs_max
  !
  implicit none
  !
  complex(SP) :: small_a_frequency
  !
  type(ext_field)   :: E_field
  complex(SP)       :: W_i
  integer           :: iErr
  !
  real(SP)          ::W_0
  complex(SP)       ::local_a(2)
  character(schlen) ::field_defs(n_fields_defs_max)
  !
  iErr=-1
  local_a=cZERO
  !
  field_defs=""
  call STRING_split(trim(E_field%ef_name),field_defs)
  !
  W_0=E_field%frequency
  select case( trim(field_defs(1)) )
  case('SIN')
    iErr=0
    local_a(1)=local_a(1)+(1._SP/(W_i-W_0)                 -1._SP/W_0)/2._SP  ! RES
    local_a(2)=local_a(2)+(               -1._SP/(W_i+W_0) -1._SP/W_0)/2._SP  ! ARES
  case('DELTA')
    iErr=0
    local_a=1._SP/2._SP
  end select
  !
  if(trim(field_defs(2))==    'RES') local_a(2)=0._SP
  if(trim(field_defs(2))=='ANTIRES') local_a(1)=0._SP
  !
  small_a_frequency=local_a(1)+local_a(2)
  !
end function small_a_frequency
