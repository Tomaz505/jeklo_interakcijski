module rutine
    private dp,pi
    public area_n

    integer, parameter :: dp = selected_real_kind(15,307)
    real(dp), parameter :: pi = 3.141592654

    contains


    subroutine area_n(xy,n,a,deg)
        integer :: n,deg,i,j
        real(dp) :: a, xy(:,:),x

        do i=1,n-1
            x = 0.0_dp
            do j= 0,deg
                x = x + xy(2,i)**(deg-j)*xy(2,i+1)**(j)
            end do
            a = a + x*(xy(1,i)*xy(2,i+1)-xy(1,i+1)*xy(2,i))/(1+deg)/(2+deg)
        end do
        x = 0.0_dp
        do j= 0,deg
            x = x + xy(2,n)**(deg-j)*xy(2,1)**(j)
        end do
        a = a + x*(xy(1,n)*xy(2,1)-xy(1,1)*xy(2,n))/(1+deg)/(2+deg)
    end subroutine

    subroutine area_ny1x(xy,n,a,deg)
        integer :: n,deg,i,j
        real(dp) :: a, xy(:,:),x

        do i=1,n-1
            x = 0.0_dp
            do j= 0,deg
                x = x + xy(2,i)**(deg-j)*xy(2,i+1)**(j)*(xy(1,i)*(deg-j+1)+xy(1,i+1)*(j+1))
            end do
            a = a + x*(xy(1,i)*xy(2,i+1)-xy(1,i+1)*xy(2,i))/((1+deg)*(2+deg)*(3+deg))
        end do
        x = 0.0_dp
        do j= 0,deg
            x = x + xy(2,n)**(deg-j)*xy(2,1)**(j)*(xy(1,1)*(deg-j+1)+xy(1,n)*(j+1))
        end do
        a = a + x*(xy(1,n)*xy(2,1)-xy(1,1)*xy(2,n))/((1+deg)*(2+deg)*(3+deg))
    end subroutine



    !Zapise potrebne kolicine v datoteko out.js
    subroutine write_js(xy,nv,rd,nrd,an)
        !n-stevilo vozlisc
        !deg-stopnja potence y^deg v integralu
        !a - vhodna količina
        !xy - koordinate vozlisc

        integer :: n,nrd,an
        real(dp) :: xy(:,:),rd(:,:)

        open (unit = 1,file = "out.js",status = "old")
        write(1,*) "analiza = ", an
        write(1,*) "const sec_coor = `"
        do i=1,nv
            write(1,*) xy(:,i)
        end do
        write(1,*) xy(:,1)
        write(1,*) "`;"

        write(1,*) "const intr = `"
        do i=1,nrd
            write(1,*) rd(:,i)
        end do
        write(1,*) "`;"

        close(1)
    end subroutine





    !Doloci del prereza, ki je med y koordinatama ymin in ymax
    !xy_sub je mnogokotnik po katerem lahko integriramo napetosti
    subroutine eff_section(xy,n,ymin,ymax,xy_sub)
        logical :: state
        integer :: n,i1
        real(dp) :: xy(:,:),ymin,ymax,xy_sub(:,:),x


        ! ali je je prvo vozlisce v prerezu
        state = ((xy(2,1) >= ymin) .and. (xy(2,1) <= ymax))
        xy_sub(:,:) = 0.0_dp
        !print*,"state = ",state
        i1 = 1
        if (state) then
            xy_sub(:,1) = xy(:,1)
            i1=i1+1
        end if


        do i=1,n-1
            if (state) then
                !sem v prerezu. preverim stanje naslednjega vozlišča
                if ((xy(2,i+1)>=ymin) .and. (xy(2,i+1)<=ymax)) then
                    xy_sub(:,i1) = xy(:,i+1)
                    i1 = i1+1
                else if (xy(2,i+1)>ymax) then
                    x = (ymax*(xy(1,i+1)-xy(1,i))+xy(2,i+1)*xy(1,i)-xy(2,i)*xy(1,i+1))/(xy(2,i+1)-xy(2,i))
                    xy_sub(:,i1) = (/x,ymax/)
                    i1 = i1+1
                    state = (.not. state)
                else if (xy(2,i+1)<ymin) then
                    x = (ymin*(xy(1,i+1)-xy(1,i))+xy(2,i+1)*xy(1,i)-xy(2,i)*xy(1,i+1))/(xy(2,i+1)-xy(2,i))
                    xy_sub(:,i1) = (/x,ymin/)
                    i1 = i1+1
                    state = (.not. state)
                end if

            else
                if ((xy(2,i)<ymin) .and. (xy(2,i+1)>ymin)) then
                    x = (ymin*(xy(1,i+1)-xy(1,i))+xy(2,i+1)*xy(1,i)-xy(2,i)*xy(1,i+1))/(xy(2,i+1)-xy(2,i))

                    xy_sub(:,i1) = (/x, ymin/)
                    i1 = i1+1

                    if (xy(2,i+1)> ymax) then
                        x = (ymax*(xy(1,i+1)-xy(1,i))+xy(2,i+1)*xy(1,i)-xy(2,i)*xy(1,i+1))/(xy(2,i+1)-xy(2,i))
                        xy_sub(:,i1) = (/x, ymax/)
                        i1 = i1+1
                    else
                        xy_sub(:,i1) = xy(:,i+1)
                        state = (.not. state)
                        i1 = i1+1
                    end if

                else if ((xy(2,i)>ymax) .and. (xy(2,i+1)<ymax)) then
                    x = (ymax*(xy(1,i+1)-xy(1,i))+xy(2,i+1)*xy(1,i)-xy(2,i)*xy(1,i+1))/(xy(2,i+1)-xy(2,i))

                    xy_sub(:,i1) = (/x, ymax/)

                    if (xy(2,i+1)< ymin) then
                        x = (ymin*(xy(1,i+1)-xy(1,i))+xy(2,i+1)*xy(1,i)-xy(2,i)*xy(1,i+1))/(xy(2,i+1)-xy(2,i))
                        xy_sub(:,i1+1) = (/x, ymin/)
                    else
                        xy_sub(:,i1+1) = xy(:,i+1)
                        state = (.not. state)
                    end if
                    i1 = i1+2
                end if
            end if
        end do


            if (state) then
                !sem v prerezu. preverim stanje naslednjega vozlišča
                if ((xy(2,1)>=ymin) .and. (xy(2,1)<=ymax)) then
                    xy_sub(:,i1) = xy(:,1)
                    i1 = i1+1
                else if (xy(2,1)>ymax) then
                    x = (ymax*(xy(1,1)-xy(1,n))+xy(2,1)*xy(1,n)-xy(2,n)*xy(1,1))/(xy(2,1)-xy(2,n))
                    xy_sub(:,i1) = (/x,ymax/)
                    i1 = i1+1
                    state = (.not. state)
                else if (xy(2,1)<ymin) then
                    x = (ymin*(xy(1,1)-xy(1,n))+xy(2,1)*xy(1,n)-xy(2,n)*xy(1,1))/(xy(2,1)-xy(2,n))
                    xy_sub(:,i1) = (/x,ymin/)
                    i1 = i1+1
                    state = (.not. state)
                end if

            else
                if ((xy(2,n)<ymin) .and. (xy(2,1)>ymin)) then
                    x = (ymin*(xy(1,1)-xy(1,n))+xy(2,1)*xy(1,n)-xy(2,n)*xy(1,1))/(xy(2,1)-xy(2,n))

                    xy_sub(:,i1) = (/x, ymin/)
                    i1 = i1+1

                    if (xy(2,1)> ymax) then
                        x = (ymax*(xy(1,1)-xy(1,n))+xy(2,1)*xy(1,n)-xy(2,n)*xy(1,1))/(xy(2,1)-xy(2,n))
                        xy_sub(:,i1) = (/x, ymax/)
                        i1 = i1+1
                    else
                        xy_sub(:,i1) = xy(:,1)
                        state = (.not. state)
                        i1 = i1+1
                    end if

                else if ((xy(2,n)>ymax) .and. (xy(2,1)<ymax)) then
                    x = (ymax*(xy(1,1)-xy(1,n))+xy(2,1)*xy(1,n)-xy(2,n)*xy(1,1))/(xy(2,1)-xy(2,n))

                    xy_sub(:,i1) = (/x, ymax/)

                    if (xy(2,1)< ymin) then
                        x = (ymin*(xy(1,1)-xy(1,n))+xy(2,1)*xy(1,n)-xy(2,n)*xy(1,1))/(xy(2,1)-xy(2,n))
                        xy_sub(:,i1+1) = (/x, ymin/)
                    else
                        xy_sub(:,i1+1) = xy(:,1)
                        state = (.not. state)
                    end if
                    i1 = i1+2
                end if
            end if
        xy_sub(:,i1) = xy_sub(:,1)

    end subroutine



end module


