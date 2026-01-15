program main
    use rutine


    integer, parameter :: dp = selected_real_kind(15,307)
    real(dp), parameter :: pi = 3.141592654

    integer :: nv, ned, analiza,u  !st vozlisc, stevilo armaturnih palic, stevilo kablov
    real(dp) :: fy, es
    real(dp), allocatable :: xy(:,:),ed(:,:)

    real(dp) :: a0,a2



    !Branje podatkov
    open(newunit = u, file = "in.txt",status = "old")
    read(u,*) analiza, fy, es, nv,ned
    close(u)

    if (ned == 0) then
        allocate(xy(2,nv))
        open(newunit = u, file = "in.txt",status = "old")
        read(u,*) analiza, fy, es, nv,ned,xy
        close(u)
    else
        allocate(xy(2,nv),ed(3,ned))
        open(newunit = u, file = "in.txt",status = "old")
        read(u,*) analiza, fy, es, nv,ned,xy,ed
        close(u)
    end if



    !RAČUN KARAKTERISTIK IN PREMIK TEZISCA V IZODISCE
    block
        real(dp) :: sx=0.0_dp,sy = 0.0_dp,xyt(2,nv)
        a0 = 0.0_dp
        a2 = 0.0_dp


        xyt(1,:) = xy(2,:)
        xyt(2,:) = xy(1,:)

        call area_n(xy,nv,a0,0)
        call area_n(xy,nv,sx,1)
        call area_n(xyt,nv,sy,1)

        xy(2,:) = xy(2,:)-sx/a0
        xy(1,:) = xy(1,:)+sy/a0

        call area_n(xy,nv,a2,2)

    end block

    goto (10,20,30), analiza




    !RAČUN ELASTICNEGA
 10 block
        real(dp) :: rd(3,20402),phi, r, def_pl(2),y_extr(2),y_y, xyt(2,nv),i0,i1,i2,h
        real(dp) :: eps_y
        eps_y = fy/es


        !print *, eps_y

        do i=0,100
            phi = 2.0_dp*pi*(i/100.0_dp)

            xyt(1,:) = xy(1,:)*cos(phi) + xy(2,:) *sin(phi)
            xyt(2,:) = xy(2,:)*cos(phi) - xy(1,:) *sin(phi)

            y_extr = (/minval(xyt(2,:)),maxval(xyt(2,:))/)
            h = y_extr(2)-y_extr(1)

            !Karakteristike rotiranega prereza
            i0 = 0.0_dp
            i1 = 0.0_dp
            i2 = 0.0_dp

            call area_n(xyt,nv,i0,0)
            call area_n(xyt,nv,i1,1)
            call area_n(xyt,nv,i2,2)

            do j = 0,100

                def_pl(2) = 2.0_dp*eps_y/h*((100-j)/100.0_dp)
                def_pl(1) = eps_y-def_pl(2)*y_extr(2)

                rd(:,(i)*100+(j+1)) = es*(/def_pl(1)*i0+def_pl(2)*i1, sin(phi)*(def_pl(1)*i1+def_pl(2)*i2) , cos(phi)*( def_pl(1)*i1+def_pl(2)*i2)/)

                rd(:,10201+(i)*100+(j+1)) = -rd(:,(i)*100+(j+1))
            end do
        end do


        call write_js(xy,nv,rd,20402,analiza)
    end block
    goto 110


    !RAČUN PLASTICNEGA
20  block
        real(dp) :: rd(3,20401),phi, r, def_pl(2),y_extr(2),y_y, xyt(2,nv),i0,i1,i2,xy_eff(2,2*nv)
        real(dp) :: eps_y

        eps_y = fy/es
        xy_eff(:,:) = 0.0_dp


        do i=0,100
            phi = 2.0_dp*pi*(i/100.0_dp)
            !print *, phi
            xyt(1,:) = xy(1,:)*cos(phi) + xy(2,:) *sin(phi)
            xyt(2,:) = xy(2,:)*cos(phi) - xy(1,:) *sin(phi)

            y_extr = (/minval(xyt(2,:)),maxval(xyt(2,:))/)

            do j = 0,100
                r = y_extr(1)*(j/100.0_dp) + y_extr(2)*((100-j)/100.0_dp)
                i0 = 0.0_dp
                i1 = 0.0_dp


                xy_eff(:,:) = 0.0_dp
                call eff_section(xyt,nv,r,y_extr(2),xy_eff)

                call area_n(xy_eff,nv*2,i0,0)
                call area_n(xy_eff,nv*2,i1,1)


                rd(:,(i)*100+(j+1)) = fy*(/i0, cos(phi)*(i1) , sin(phi)*( i1)/)


                xy_eff(:,:) = 0.0_dp
                call eff_section(xyt,nv,y_extr(1),r,xy_eff)




                i0 = 0.0_dp
                i1 = 0.0_dp

                call area_n(xy_eff,nv*2,i0,0)
                call area_n(xy_eff,nv*2,i1,1)


                rd(:,(i)*100+(j+1)) = rd(:,(i)*100+(j+1)) - fy*(/i0, cos(phi)*(i1) , sin(phi)*( i1)/)
            end do
        end do

        call write_js(xy,nv,rd,20401,analiza)
    end block
goto 110


30  block
        real(dp) :: rd(6,20401),phi, r, def_pl(2), y_extr(2), y_y, xyt(2,nv), i0, i1, i2, xy_eff(2,2*nv)
        real(dp) :: i0el, i1el, i2el, h, eps_y

        eps_y = fy/es
        xy_eff(:,:) = 0.0_dp


        do i=0,100
            phi = 2.0_dp*pi*(i/100.0_dp)
            !print *, phi
            xyt(1,:) = xy(1,:)*cos(phi) + xy(2,:) *sin(phi)
            xyt(2,:) = xy(2,:)*cos(phi) - xy(1,:) *sin(phi)

            y_extr = (/minval(xyt(2,:)),maxval(xyt(2,:))/)
            h = y_extr(2)-y_extr(1)

            i0el = 0.0_dp
            i1el = 0.0_dp
            i2el = 0.0_dp

            call area_n(xyt,nv,i0el,0)
            call area_n(xyt,nv,i1el,1)
            call area_n(xyt,nv,i2el,2)


            do j = 0,100


            end do


            do j = 0,100

                !ELASTICNA
                def_pl(2) = 2.0_dp*eps_y/h*((100-j)/100.0_dp)
                def_pl(1) = eps_y-def_pl(2)*y_extr(2)
                rd(1:3,(i)*100+(j+1)) = es*(/def_pl(1)*i0el+def_pl(2)*i1el, sin(phi)*(def_pl(1)*i1el+def_pl(2)*i2el) , cos(phi)*( def_pl(1)*i1el+def_pl(2)*i2el)/)
                rd(1:3,10201+(i)*100+(j+1)) = -rd(1:3,(i)*100+(j+1))



                !PLASTICNA
                r = y_extr(1)*(j/100.0_dp) + y_extr(2)*((100-j)/100.0_dp)

                i0 = 0.0_dp
                i1 = 0.0_dp
                xy_eff(:,:) = 0.0_dp
                call eff_section(xyt,nv,r,y_extr(2),xy_eff)
                call area_n(xy_eff,nv*2,i0,0)
                call area_n(xy_eff,nv*2,i1,1)
                rd(4:6,(i)*100+(j+1)) = fy*(/i0, cos(phi)*(i1) , sin(phi)*( i1)/)


                xy_eff(:,:) = 0.0_dp
                call eff_section(xyt,nv,y_extr(1),r,xy_eff)
                i0 = 0.0_dp
                i1 = 0.0_dp
                call area_n(xy_eff,nv*2,i0,0)
                call area_n(xy_eff,nv*2,i1,1)

                rd(4:6,(i)*100+(j+1)) = rd(4:6,(i)*100+(j+1)) - fy*(/i0, cos(phi)*(i1) , sin(phi)*( i1)/)
                !d(4:6,10201+(i)*100+(j+1)) = -rd(4:6,(i)*100+(j+1))
            end do
        end do

        call write_js(xy,nv,rd,20401,analiza)
    end block
goto 110


110     print *," "
        print*,"    Račun končan."
        !print*," "








end program
