module IO


use parameters
use constants


implicit none


private

public IO_array


contains


subroutine IO_array(array,rows,cols)
integer(kind=dp), intent(in) :: rows,cols
real(kind=dp), intent(in) :: array(rows,cols)


print *, array(rows,:)

print *, 'subroutine called ok'

end subroutine IO_array



end module IO
