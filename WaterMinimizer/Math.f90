! Math routines
subroutine cross_product_3d(v1,v2,v3)
!    O produtro cruzado é o determinante da matriz
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
    implicit none
    real, intent(in)    :: v1(3)
    real, intent(in)    :: v2(3)
    real, intent(out)   :: v3(3)

    v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
    v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
    v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

    return
end subroutine cross_product_3d
