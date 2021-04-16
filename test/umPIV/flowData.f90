program flowdata
    parameter(il=376, jl=147, kl=159, knode=il*jl*kl)
    integer k
    common q(5, knode)
    open(unit=7, file='data/flowData.txt', status='old',&
        & form='formatted')
    read(7, *) ((q(n,k), n=1,5), k=1,knode)
    close(7)
    open(unit=77, file='flowData.o', status='unknown',&
        & form='unformatted')
    write(77) 1.4, 286.9
    write(77) ((q(n,k), n=1,5), k=1,knode)
    close(77)
end program flowdata
