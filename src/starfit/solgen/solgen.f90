program test
    implicit none
end program test

subroutine choice(db_size, number, sol_size, indices)
    implicit none
    integer(8),        intent(in)      :: db_size
    integer(8),        intent(in)      :: number
    integer(8),        intent(in)      :: sol_size

    !f2py integer(8),  intent(out)     :: indices(number, sol_size)
    integer(8),        intent(out)     :: indices(number, sol_size)

    logical                            :: ok
    integer(8)                         :: i, j, k, s, smax

    do i = 1, number
        do j = 1, sol_size
            ok = .False.
            indices(i, j) = -1
            do while (.not.ok)
                s = int(rand(0) * db_size)
                ok = .True.
                do k = 1, j
                    if (s == indices(i, k)) then
                        ok = .False.
                    endif
                enddo
            enddo
            indices(i, j) = s
        enddo
    enddo
end subroutine choice

subroutine gen_slice(start, finish, sol_size, sols)
    implicit none
    integer(8),        intent(in)      :: start, finish
    integer(8),        intent(in)      :: sol_size

    !f2py integer(8),        intent(out)     :: sols(finish - start, sol_size)
    integer(8),        intent(out)     :: sols(finish - start, sol_size)

    integer(8)                         :: i, islice

    i = 1
    islice = start
    call gen_init(islice, sol_size, sols(i, :))

    do islice = start + 1, finish - 1
        i = i + 1
        call gen_hop(sols(i - 1, :), sol_size, sols(i, :))
    enddo
end subroutine gen_slice

subroutine gen_init(x, sol_size, sol)
    implicit none
    integer(8),        intent(in)      :: x
    integer(8),        intent(in)      :: sol_size

    !f2py integer(8),        intent(out)     :: sol(sol_size)
    integer(8),        intent(out)     :: sol(sol_size)

    integer(8)                         :: i, n, n_old, sol_val, xr
    integer(8)                         :: rem_size

    rem_size = sol_size
    xr = x
    do i = 1, sol_size
        n = 0
        sol_val = rem_size - 1
        do while (.true.)
            n_old = n
            call comb(sol_val, rem_size, n)
            if (n > xr) then
                sol(i) = sol_val - 1
                xr = xr - n_old
                exit
            endif
            sol_val = sol_val + 1
        enddo
        rem_size = rem_size - 1
    enddo

end subroutine gen_init

subroutine gen_hop(sol, sol_size, new_sol)
    implicit none
    integer(8),        intent(in)      :: sol_size
    integer(8),        intent(in)      :: sol(sol_size)

    !f2py integer(8),        intent(out)     :: new_sol(sol_size)
    integer(8),        intent(out)     :: new_sol(sol_size)

    integer(8)                         :: i, j
    integer(8)                         :: c(sol_size)

    c(1:sol_size) = sol(sol_size:1:-1)
    i = sol_size + 1
    do j = 2, sol_size
        if (c(j) >= c(j - 1) + 2) then
            i = j
            exit
        endif
    enddo
    c(i - 1) = c(i - 1) + 1
    do j = 1, i - 2
        c(j) = j - 1
    enddo
    new_sol(1:sol_size) = c(sol_size:1:-1)

end subroutine gen_hop

subroutine comb(n, r, nc)
    implicit none
    integer(8),        intent(in)      :: n
    integer(8),        intent(in)      :: r

    !f2py integer(8),        intent(out)     :: nc
    integer(8),        intent(out)     :: nc

    integer(8)                         :: i, d

    if (r == 0) then
        nc = 1
    else
        nc = 1
        do i = n - r + 1, n
            nc = nc * i
        enddo

        d = 1
        do i = 1, r
            d = d * i
        enddo

        nc = nc / d

    endif

end subroutine comb
