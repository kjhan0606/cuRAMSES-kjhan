!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character(LEN=5)::nchar
  write(nchar,'(i05)') n
end subroutine title

