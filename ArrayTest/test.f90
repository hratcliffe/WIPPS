
program test

do i = 0, 100
  call temp
end do



end program


subroutine temp
  integer :: i = 0
  i = i + 1
  write(*,*) i
end subroutine temp
