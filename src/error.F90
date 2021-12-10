module error
  use header
implicit none
interface face_error
module procedure print_error1
module procedure print_error2
module procedure print_error3
module procedure print_error4
module procedure print_error5
module procedure print_error6
module procedure print_error7
module procedure print_error8
module procedure print_error9
module procedure print_error10
module procedure print_error11
module procedure print_error12
module procedure print_error13
module procedure print_error14
module procedure print_error15
end interface face_error
interface face_warning
module procedure print_warning1
module procedure print_warning2
module procedure print_warning3
module procedure print_warning4
module procedure print_warning5
module procedure print_warning6
module procedure print_warning7
module procedure print_warning8
module procedure print_warning9
module procedure print_warning10
end interface face_warning
character*256:: message_exit_error="Abnormal exit of FACE..."
contains
! errors
subroutine print_error1(str)
character(*)::str
write(message_exit_error,*)    str
error_status=1
if (enforce_error) then
call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif

end subroutine print_error1

subroutine print_error2(str,r)
real :: r
character(*)::str
write(message_exit_error,*)    str,r
error_status=1
if (enforce_error) then
 
 call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif
end subroutine print_error2

subroutine print_error3(str,i)
integer :: i
character(*)::str
write(message_exit_error,*)    str,i
error_status=1
if (enforce_error) then
 
 call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif
end subroutine print_error3

subroutine print_error4(str,r,str2,r2)
real :: r,r2
character(*)::str,str2
write(message_exit_error,*)    str,r,str2,r2
error_status=1
if (enforce_error) then
 
 call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif
end subroutine print_error4

subroutine print_error5(str,i,str2,i2)
integer :: i,i2
character(*)::str,str2
write(message_exit_error,*)    str,i,str2,i2
error_status=1
if (enforce_error) then
 
 call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif
end subroutine print_error5

subroutine print_error6(str,str2)
character(*)::str,str2
write(message_exit_error,*)    str,str2
error_status=1
if (enforce_error) then
 
 call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif

end subroutine print_error6

subroutine print_error7(str,i,str2,r,str3,r2)
integer ::i
real :: r,r2
character(*)::str,str2,str3
write(message_exit_error,*)    str,i,str2,r,str3,r2
error_status=1
if (enforce_error) then
 
 call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif
end subroutine print_error7

subroutine print_error8(str,i,str2,i2,str3,r,str4,r2)
integer ::i,i2
real :: r,r2
character(*)::str,str2,str3,str4
write(message_exit_error,*)    str,i,str2,i2,str3,r,str4,r2
error_status=1
if (enforce_error) then
 
 call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif
end subroutine print_error8

subroutine print_error9(str,str2,str3)
character(*)::str,str2,str3
write(message_exit_error,*)    str,str2,str3
error_status=1
if (enforce_error) then
 
 call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif
end subroutine print_error9

subroutine print_error10(str,str2,str3,i)
character(*)::str,str2,str3
integer :: i
write(message_exit_error,*)    str,str2,str3,i
error_status=1
if (enforce_error) then
 
 call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif
end subroutine print_error10

subroutine print_error11(str,i,str2,i2,str3,str4,str5,str6)
character(*)::str,str2,str3,str4,str5,str6
integer :: i,i2
write(message_exit_error,*)    str,i,str2,i2,str3,str4,str5,str6
error_status=1
if (enforce_error) then
 
 call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif
end subroutine print_error11

subroutine print_error12(str,i,str2,r)
real :: r
integer::i
character(*)::str,str2
write(message_exit_error,*)    str,i,str2,r
error_status=1
if (enforce_error) then
 
 call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif
end subroutine print_error12

subroutine print_error13(str,str2,r,str3,r2)
real :: r,r2
character(*)::str,str2,str3
write(message_exit_error,*)    str,str2,r,str3,r2
error_status=1
if (enforce_error) then
 
 call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif
end subroutine print_error13
subroutine print_error14(str,i,str2,i2,str3,i3)
character(*)::str,str2,str3
integer::i,i2,i3
write(message_exit_error,*)    str,i,str2,i2,str3,i3
error_status=1
if (enforce_error) then
 
 call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif
end subroutine print_error14
subroutine print_error15(str,str2,str3,i,str4)
character(*)::str,str2,str3,str4
integer::i
write(message_exit_error,*)    str,str2,str3,i,str4
error_status=1
if (enforce_error) then
 
 call kaboom(message_exit_error)
return
else
 write(iout,*) "FACE ERROR:", message_exit_error
endif
end subroutine print_error15
! warning
subroutine print_warning1(str)
character(*)::str
write(message_exit_error,*)  "Warning:", str

end subroutine print_warning1

subroutine print_warning2(str,r)
real :: r
character(*)::str
write(message_exit_error,*)  "Warning:", str,r

end subroutine print_warning2

subroutine print_warning3(str,i)
integer :: i
character(*)::str
write(message_exit_error,*)  "Warning:", str,i

end subroutine print_warning3

subroutine print_warning4(str,r,str2,r2)
real :: r,r2
character(*)::str,str2
write(message_exit_error,*)  "Warning:", str,r,str2,r2

end subroutine print_warning4

subroutine print_warning5(str,i,str2,i2)
integer :: i,i2
character(*)::str,str2
write(message_exit_error,*)  "Warning:", str,i,str2,i2

end subroutine print_warning5

subroutine print_warning6(str,str2)
character(*)::str,str2
write(message_exit_error,*)  "Warning:", str,str2

end subroutine print_warning6

subroutine print_warning7(str,i,str2,r,str3,r2)
integer ::i
real :: r,r2
character(*)::str,str2,str3
write(message_exit_error,*)  "Warning:", str,i,str2,r,str3,r2

end subroutine print_warning7

subroutine print_warning8(str,i,str2,i2,str3,r,str4,r2)
integer ::i,i2
real :: r,r2
character(*)::str,str2,str3,str4
write(message_exit_error,*)  "Warning:", str,i,str2,i2,str3,r,str4,r2

end subroutine print_warning8

subroutine print_warning9(str,i,str2,i2,str3,r)
integer ::i,i2
real :: r
character(*)::str,str2,str3
write(message_exit_error,*)  "Warning:", str,i,str2,i2,str3,r

end subroutine print_warning9

subroutine print_warning10(str,i,str2,r)
integer ::i
real :: r
character(*)::str,str2
write(message_exit_error,*)  "Warning:", str,i,str2,r

end subroutine print_warning10

end module error
