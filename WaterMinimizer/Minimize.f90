subroutine Minimize(coord,force)
    real,intent(in)     :: force(3,3)
    real,intent(inout)  :: coord(3,3)
    do i=1,3
        coord(i,1) = coord(i,1)-force(i,1)*0.00001
        coord(i,2) = coord(i,2)-force(i,2)*0.00001
        coord(i,3) = coord(i,3)-force(i,3)*0.00001
    enddo
end subroutine Minimize

