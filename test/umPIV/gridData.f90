program gridgen
    integer i, j, k, knode, kc
    real, dimension(:), allocatable :: a, b, c
    parameter (il=376, jl=147, kl=159, irv=0,&
        & nb=36, ile=9, ite=33,&
        & chord=33.893616,mnod=2000)
    knode=il*jl*kl
    allocate(a(knode), b(knode), c(knode))
    open(unit=7, file='data/gridData.txt', status='old',&
    & form='formatted')
    read(7, *) (a(k),b(k),c(k), k=1,knode)
    close(7)

    open(unit=77, file='gridData.o', status='unknown',&
    & form='unformatted')
    write(77)  il, jl, kl
    write(77) (a(k)*1e-3, b(k)*1e-3, c(k)*1e-3, k=1,knode)
    deallocate(a,b,c)
    close(77)
end program gridgen
