program main
    use rutine


    integer, parameter :: dp = selected_real_kind(15,307)
    real(dp), parameter :: pi = 3.141592654

    integer :: nv, ned, analiza,u  !st vozlisc, stevilo armaturnih palic, stevilo kablov
    real(dp) :: fy
    real(dp), allocatable :: xy(:,:),ed(:,:)

    real(dp) :: a0,a2



    !Branje podatkov
    open(newunit = u, file = "in.txt",status = "old")
    read(u,*) analiza, fy, nv,ned
    close(u)

    if (ned == 0) then
        allocate(xy(2,nv))
        open(newunit = u, file = "in.txt",status = "old")
        read(u,*) analiza, fy, nv,ned,xy
        close(u)
    else
        allocate(xy(2,nv),ed(3,ned))
        open(newunit = u, file = "in.txt",status = "old")
        read(u,*) analiza, fy, nv,ned,xy,ed
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
        real(dp) :: rd(3,20402),phi, r, def_pl(2),y_extr(2),y_y, xyt(2,nv),i0,i2,h,ix1,rot(2,2)
        real(dp) :: eps_y
        eps_y = fy

        !print *, eps_y
        i0 = 0.0_dp
        call area_n(xy,nv,i0,0)
        do i=0,100
            phi = 2.0_dp*pi*(i/100.0_dp)

            rot(1,:)  = (/cos(phi), sin(phi)/)
            rot(2,:)  = (/-sin(phi), cos(phi)/)

            xyt = matmul(rot,xy)

            y_extr = (/minval(xyt(2,:)),maxval(xyt(2,:))/)
            h = y_extr(2)-y_extr(1)

            !Karakteristike rotiranega prereza
            i2 = 0.0_dp
            ix1 = 0.0_dp

            call area_n(xyt,nv,i2,2)
            call area_ny1x(xyt,nv,ix1,1)

            do j = 0,100
                def_pl(2) = 2.0_dp*eps_y/h*((100-j)/100.0_dp)
                def_pl(1) = eps_y-def_pl(2)*y_extr(2)
                 rd(:,(i)*100+(j+1)) = (/def_pl(1)*i0, def_pl(2)*i2 , def_pl(2)*ix1 /)

            end do
            rd(2:3,i*100+1:i*100+101) = matmul(rot,rd(2:3,i*100+1:i*100+101))
            rd(:,10202+(i)*100:10302+i*100) = -rd(:,i*100+1:i*100+101)
        end do


        call write_js(xy,nv,rd,20402,analiza)
    end block
    goto 110


    !RAČUN PLASTICNEGA
20  block
        real(dp) :: rd(3,20401),phi, r, def_pl(2),y_extr(2),y_y
        real(dp) :: eps_y,i0,i1,ix0,ia
        real(dp) ::  xyt(2,nv),xy_eff(2,2*nv), rot(2,2)

        eps_y = fy/es

        ia = 0.0_dp
        call area_n(xy,nv,ia,0)



        do i=0,100
            phi = 2.0_dp*pi*(i/100.0_dp)
            rot(1,:)  = (/cos(phi), sin(phi)/)
            rot(2,:)  = (/-sin(phi), cos(phi)/)

            xyt = matmul(rot,xy)
            y_extr = (/minval(xyt(2,:)),maxval(xyt(2,:))/)






            do j=0,100
                !PLASTICNA
                r = y_extr(1)*(j/100.0_dp) + y_extr(2)*((100-j)/100.0_dp)

                i0 = 0.0_dp
                ix0 = 0.0_dp
                i1 = 0.0_dp
                xy_eff(:,:) = 0.0_dp
                call eff_section(xyt,nv,r,y_extr(2),xy_eff)
                call area_n(xy_eff,nv*2,i0,0)
                call area_n(xy_eff,nv*2,i1,1)
                call area_ny1x(xy_eff,nv*2,ix0,0)
                print *,i0,ia,i1,ix0

                rd(:,(i)*100+(j+1)) = fy*(/2*i0-ia, i1*2 , ix0*2/)
            end do
            rd(2:3,i*100+1:i*100+101) = matmul(rot,rd(2:3,i*100+1:i*100+101))
        end do

        call write_js(xy,nv,rd,20401,analiza)
    end block
goto 110


30  block
        real(dp) :: rd(6,40602),phi, r, def_pl(2), y_extr(2), y_y, xyt(2,nv), i0, i1, i2, xy_eff(2,2*nv),rot(2,2)
        real(dp) :: i0el, i1el, i2el, h, eps_y,ix0,ix1

        eps_y = fy
        xy_eff(:,:) = 0.0_dp

        i0el = 0.0_dp
        call area_n(xy,nv,i0el,0)
        do i=0,200
            phi = 2.0_dp*pi*(i/200.0_dp)
            rot(1,:)  = (/cos(phi), sin(phi)/)
            rot(2,:)  = (/-sin(phi), cos(phi)/)

            xyt = matmul(rot,xy)

            y_extr = (/minval(xyt(2,:)),maxval(xyt(2,:))/)
            h = y_extr(2)-y_extr(1)



            !ELASTICNA

            i2 = 0.0_dp
            ix1 = 0.0_dp
            call area_n(xyt,nv,i2,2)
            call area_ny1x(xyt,nv,ix1,1)
            do j = 0,100
                def_pl(2) = 2.0_dp*eps_y/h*((100-j)/100.0_dp)
                def_pl(1) = eps_y-def_pl(2)*y_extr(2)
                rd(1:3,(i)*100+(j+1)) = (/def_pl(1)*i0el,def_pl(2)*i2,def_pl(2)*ix1/)
            end do

            do j=0,201
                !PLASTICNA
                r = y_extr(1)*(j/201.0_dp) + y_extr(2)*((201-j)/201.0_dp)

                i0 = 0.0_dp
                i1 = 0.0_dp
                ix0 = 0.0_dp
                xy_eff(:,:) = 0.0_dp
                call eff_section(xyt,nv,r,y_extr(2),xy_eff)
                call area_n(xy_eff,nv*2,i0,0)
                call area_n(xy_eff,nv*2,i1,1)
                call area_ny1x(xy_eff,nv*2,ix0,0)

                rd(4:6,(i)*200+(j+1)) = fy*(/2*i0-i0el, i1*2 , ix0*2/)
            end do
            rd(5:6,i*200+1:i*200+202) = matmul(rot,rd(5:6,i*200+1:i*200+202))
            rd(2:3,i*100+1:i*100+101) = matmul(rot,rd(2:3,i*100+1:i*100+101))
            rd(1:3,20302+(i)*100:20402+i*100) = -rd(1:3,i*100+1:i*100+101)

        end do


        call write_js(xy,nv,rd,40402,analiza)
    end block
goto 110


110     print *," "
        print*,"    Račun končan."
        !print*," "








end program
