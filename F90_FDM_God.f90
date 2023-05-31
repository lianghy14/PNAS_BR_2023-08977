#include "fintrf.h" !必须有的头文件,里面有mxGetPr, mxGetM, mxGetN,mxCreateDoubleMatrix等函数的申明

Subroutine mexFunction(OutSum,OutVar,InSum,InVar)!函数接口名称必须为mexFunction,
    !OutSum:输出参数个数, OutVar:输出参数数组指针, InSum:输入参数个数, InVar:输入参数数组指针

    Integer InSum,OutSum
    mwPointer InVar(*),OutVar(*)                           !mwPointer专门用于表示指针变量,这个不能随意用Integer代替
    mwPointer mxGetPr, mxGetM, mxGetN, mxCreateDoubleMatrix, mxCreateNumericArray,mxGetDimensions,mxGetNumberOfDimensions !这个对返回指针函数的再次申明,
    integer, parameter              ::  acc = 8
    Real(kind = acc),Allocatable    ::  q1(:,:,:,:),q2(:,:,:,:),boundary(:,:,:),panic(:,:,:),slope(:,:),Ve(:,:,:,:),P2(:,:),F_P2(:,:,:),q_in(:,:)
    real(kind = acc)                ::  h,tstep,alpha_LF,nit_r,bj_out(2)
    Integer                         ::  ny,nx,n_OD,nit,dim_Q(4)

    If(InSum/=8)Then
        call mexErrMsgIdAndTxt('MATLAB:InputTooBig','Less than 8 inputs.')
        Return
    EndIf

    CALL mxCopyPtrToInteger8(mxGetDimensions(InVar(1)), dim_Q,mxGetNumberOfDimensions(InVar(1)))
    n_OD    = dim_Q(4)
    ny      = dim_Q(1)
    nx      = dim_Q(2)
    Allocate(q1(ny,nx,3,n_OD),q2(ny,nx,3,n_OD),boundary(ny,nx,n_OD),panic(ny,nx,n_OD),slope(ny,nx),Ve(ny,nx,2,n_OD),P2(ny,nx),F_P2(ny,nx,2),q_in(3,n_OD))

    Call mxCopyPtrToReal8(mxGetPr(InVar(1)),q1,ny*nx*3*n_OD)
    Call mxCopyPtrToReal8(mxGetPr(InVar(2)),bj_out,2)
    Call mxCopyPtrToReal8(mxGetPr(InVar(3)),boundary,ny*nx*n_OD)
    Call mxCopyPtrToReal8(mxGetPr(InVar(4)),panic,ny*nx*n_OD)
    Call mxCopyPtrToReal8(mxGetPr(InVar(5)),slope,ny*nx)
    Call mxCopyPtrToReal8(mxGetPr(InVar(6)),q_in,3*n_OD)
    Call mxCopyPtrToReal8(mxGetPr(InVar(7)),h,1)
    Call mxCopyPtrToReal8(mxGetPr(InVar(8)),tstep,1)

    Call FDM_euler(q1,q2,bj_out,boundary,panic,slope,q_in,h,tstep,ny,nx,n_OD,alpha_LF,nit,P2,F_P2,Ve)
    ! q2 = q1*2d0

    nit_r = real(nit)

    OutVar(1) = mxCreateNumericArray(4,dim_Q, mxClassIDFromClassName('double'),0);
    OutVar(2) = mxCreateNumericArray(4,(/ny,nx,2,n_OD/), mxClassIDFromClassName('double'),0);
    OutVar(3) = mxCreateNumericArray(3,(/ny,nx,2/), mxClassIDFromClassName('double'),0);
    OutVar(4) = mxCreateDoubleMatrix(ny,nx,0);
    OutVar(5) = mxCreateDoubleMatrix(1,1,0)!给返回参数分配内存
    OutVar(6) = mxCreateDoubleMatrix(1,1,0)!给返回参数分配内存

    Call mxCopyReal8ToPtr(q2,mxGetPr(OutVar(1)),ny*nx*3*n_OD)
    Call mxCopyReal8ToPtr(Ve,mxGetPr(OutVar(2)),ny*nx*n_OD*2)
    Call mxCopyReal8ToPtr(F_P2,mxGetPr(OutVar(3)),ny*nx*2)
    Call mxCopyReal8ToPtr(P2,mxGetPr(OutVar(4)),ny*nx)
    Call mxCopyReal8ToPtr(alpha_LF,mxGetPr(OutVar(5)),1)
    Call mxCopyReal8ToPtr(nit_r,mxGetPr(OutVar(6)),1)

    DeAllocate(q1,q2,Ve,P2,F_P2,boundary,panic,slope,q_in)!释放临时分配的内存

    Return

End SubRoutine


SubRoutine FDM_euler(q1,q2,bj_out,bj_real,panic,slope,q_in,h,tstep,ny,nx,n_OD,alpha_LF,nit,P2,F_P2,Ve)

    implicit none
    integer, parameter          ::  acc = 8
    real(kind=acc),parameter    ::  vel_f = 1.034d0, gama_c = -0.075d0, gama_p = -0.045d0, gama_2 = -0.019d0,&
                                    &c0 = 0.6d0, rho_0 = 6d0, rho_1 = 7d0, rho_m = 10d0, tao_c = 5d0, tao_p = 0.5d0,&
                                    &m_ave = 65d0, p2capacity = 200d0

    integer, Intent(In)         ::  n_OD, ny, nx

    real(kind=acc), Intent(In)  ::  q1(ny,nx,3,n_OD), q_in(3,n_OD), h, tstep, bj_out(2), bj_real(ny,nx,n_OD), panic(ny,nx,n_OD),slope(ny,nx)

    real(kind=acc), Intent(Out) ::  q2(ny,nx,3,n_OD), alpha_LF, P2(ny,nx), F_P2(ny,nx,2)

    real(kind=acc)              ::  gama_1(ny,nx),tao(ny,nx,n_OD),cost(ny,nx),potential(ny,nx),p_x(ny,nx),p_y(ny,nx),vel_mag(ny,nx),&
                                    &F_cell(ny,nx,3),G_cell(ny,nx,3),Fcell_mid(ny,nx,3),Gcell_mid(ny,nx,3),alpha_F,alpha_G,alpha_FG(n_OD),q_xex(ny,nx,3),&
                                    &q_yex(ny,nx,3),rho_ex(ny,nx,3,n_OD),ve_dir(ny,nx,2,n_OD),&
                                    &Ve(ny,nx,2,n_OD),rho_com(ny,nx), ve_com(ny,nx,2), cosangle(ny,nx), c_P1x(ny,nx),c_P1y(ny,nx), &
                                    &P1_xex(ny,nx), P1_yex(ny,nx), rho_x, rho_y, angle, alpha(ny,nx), conflict(ny,nx), pforce(ny,nx), debug(10,10),&
                                    &q_out, v_out

    integer                     ::  i_OD, j_OD, i, j, n, bound_all(ny,nx,n_OD), boundary(ny,nx), nit, region(ny,nx), bound_out(2)

    bound_all       =   NINT(bj_real)
    bound_out       =   NINT(bj_out)
    rho_com         =   0d0
    ve_com          =   0d0
    P2              =   0d0
    
    do i_OD = 1,n_OD,1
        boundary              =   bound_all(:,:,i_OD)
        rho_ex(:,:,:,i_OD)    =   q1(:,:,:,i_OD)
        do i = 1,ny,1
            do j = 1,nx,1
                ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
                if ((boundary(i,j)>2).AND.(i>3).AND.(j>3).AND.(i<ny-2).AND.(j<ny-2)) then
                    if (boundary(i,j+1)==2) then
                        rho_ex(i, j+1:j+3:1, 1,i_OD)      = q_in(1,i_OD)
                        rho_ex(i, j+1:j+3:1, 2,i_OD)      = q_in(2,i_OD)
                        rho_ex(i, j+1:j+3:1, 3,i_OD)      = q_in(3,i_OD)
                    end if
                    if (boundary(i,j-1)==2) then
                        rho_ex(i, j-3:j-1:1, 1,i_OD)      = q_in(1,i_OD)
                        rho_ex(i, j-3:j-1:1, 2,i_OD)      = q_in(2,i_OD)
                        rho_ex(i, j-3:j-1:1, 3,i_OD)      = q_in(3,i_OD)
                    end if
                    if (boundary(i+1,j)==2) then
                        rho_ex(i+1:i+3:1, j ,1,i_OD)      = q_in(1,i_OD)
                        rho_ex(i+1:i+3:1, j ,2,i_OD)      = q_in(2,i_OD)
                        rho_ex(i+1:i+3:1, j ,3,i_OD)      = q_in(3,i_OD)
                    end if
                    if (boundary(i-1,j)==2) then
                        rho_ex(i-3:i-1:1, j, 1,i_OD)      = q_in(1,i_OD)
                        rho_ex(i-3:i-1:1, j, 2,i_OD)      = q_in(2,i_OD)
                        rho_ex(i-3:i-1:1, j, 3,i_OD)      = q_in(3,i_OD)
                    end if
                end if
            end do
        end do
        rho_com =   rho_com  + rho_ex(:,:,1,i_OD)
    end do
    ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
    if (bound_out(1)>0) then
        do i = 1,ny,1
            do j = 1,nx,1
                if ((bound_all(i,j,1)==1)) then
                    rho_com(i,j) = rho_0
                end if
            end do
        end do
    end if

    do i_OD = 1,n_OD,1

        boundary        =   bound_all(:,:,i_OD)
        gama_1          =   gama_c*(1d0-panic(:,:,i_OD)) + gama_p*panic(:,:,i_OD)
        tao(:,:,i_OD)   =   tao_c*(1d0-panic(:,:,i_OD)) + tao_p*panic(:,:,i_OD)

        vel_mag     =   vel_f * exp( gama_1 * (max(rho_com,0d0)**2) )

        do j_OD = 1,n_OD,1
            cosangle = (rho_ex(:,:,2,i_OD)*rho_ex(:,:,2,j_OD) + rho_ex(:,:,3,i_OD)*rho_ex(:,:,3,j_OD)) /&
            & max(1d-12,(sqrt(rho_ex(:,:,2,i_OD)**2+rho_ex(:,:,3,i_OD)**2) * sqrt(rho_ex(:,:,2,j_OD)**2+rho_ex(:,:,3,j_OD)**2)))

            vel_mag = vel_mag * exp( gama_2*(1d0-cosangle)*max(rho_ex(:,:,1,j_OD),0d0)**2);
        end do
        
        cost    =   1d0 / vel_mag + 0.02d0 * (max(rho_com,0d0))**2

        call F90_fsmgod(cost,boundary,potential,p_x,p_y,h,ny,nx,nit)

        ve_dir(:,:,1,i_OD) = - p_x / max(1d-12,(p_x**2 + p_y**2)**0.5)
        ve_dir(:,:,2,i_OD) = - p_y / max(1d-12,(p_x**2 + p_y**2)**0.5)
        ve_com(:,:,1) = ve_com(:,:,1) + rho_ex(:,:,1,i_OD) * ve_dir(:,:,1,i_OD) / rho_com
        ve_com(:,:,2) = ve_com(:,:,2) + rho_ex(:,:,1,i_OD) * ve_dir(:,:,2,i_OD) / rho_com

        Ve(:,:,1,i_OD)   =   ve_dir(:,:,1,i_OD) * vel_mag
        Ve(:,:,2,i_OD)   =   ve_dir(:,:,2,i_OD) * vel_mag


        q_xex   =   rho_ex(:,:,:,i_OD)
        q_yex(:,:,1)   =   rho_ex(:,:,1,i_OD)
        q_yex(:,:,2)   =   rho_ex(:,:,3,i_OD)
        q_yex(:,:,3)   =   rho_ex(:,:,2,i_OD)
        ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
        do i = 1,ny,1
            do j = 1,nx,1
                if (boundary(i,j)>2) then
                    if (boundary(i,j+1)==1) then
                        q_xex(i,j+1,:)      = q_xex(i,j,:)
                    end if
                    if (boundary(i,j-1)==1) then
                        q_xex(i,j-1,:)      = q_xex(i,j,:)
                    end if
                    if (boundary(i+1,j)==1) then
                        q_yex(i+1,j,:)      = q_yex(i,j,:)
                    end if
                    if (boundary(i-1,j)==1) then
                        q_yex(i-1,j,:)      = q_yex(i,j,:)
                    end if

                    if (bound_out(1)>0) then
                        if ((i_OD == 1).OR.(i_OD == 2)) then
                            if ((boundary(i-1,j)==1)) then
                                q_out               =   q_yex(i,j,1) / max(1d-12, rho_com(i,j)) * rho_0
                                v_out               =   vel_f * exp( gama_c * (rho_0**2) )
                                q_yex(i-1,j,1)      =   q_out
                                q_yex(i-1,j,2)      =   q_out * v_out
                                q_yex(i-1,j,3)      =   0
                            end if
                        end if
                    end if   

                end if
            end do
        end do

        !   P1
        do i = 1,ny,1
            do j = 1,nx,1
                if (q_xex(i,j,1)<=rho_0) then
                    c_P1x(i,j)   = c0
                    P1_xex(i,j)     = (c0 ** 2) * max(0d0,q_xex(i,j,1))
                else
                    if ((q_xex(i,j,1)>rho_0) .AND. (q_xex(i,j,1)<rho_1)) then
                        c_P1x(i,j)   = c0 / 2d0
                        P1_xex(i,j)     = ((c0**2) / 4d0 * q_xex(i,j,1) + (c0**2) * rho_0 - (c0**2) / 4d0 * rho_0 )
                    else
                        c_P1x(i,j)   = 0d0
                        P1_xex(i,j)     = ((c0**2) / 4d0 * rho_1 + (c0**2) * rho_0 - (c0**2) / 4d0 * rho_0 )
                    end if
                end if
                if (q_yex(i,j,1)<=rho_0) then
                    c_P1y(i,j)   = c0
                    P1_yex(i,j)     = (c0 ** 2) * max(0d0,q_yex(i,j,1))
                else
                    if ((q_yex(i,j,1)>rho_0) .AND. (q_yex(i,j,1)<rho_1)) then
                        c_P1y(i,j)   = c0 / 2d0
                        P1_yex(i,j)     = ((c0**2) / 4d0 * q_yex(i,j,1) + (c0**2) * rho_0 - (c0**2) / 4d0 * rho_0 )
                    else
                        c_P1y(i,j)   = 0d0
                        P1_yex(i,j)     = ((c0**2) / 4d0 * rho_1 + (c0**2) * rho_0 - (c0**2) / 4d0 * rho_0 )
                    end if
                end if
            end do
        end do

        !  F_vector
        F_cell(:,:,1) = q_xex(:,:,2)
        F_cell(:,:,2) = q_xex(:,:,2)**2 / max(1d-12,abs(q_xex(:,:,1)))
        F_cell(:,:,3) = q_xex(:,:,2) * q_xex(:,:,3) / max(1d-12,abs(q_xex(:,:,1)))
        F_cell(:,:,2) = F_cell(:,:,2) + P1_xex
        ! % G_vector
        G_cell(:,:,1) = q_yex(:,:,2)
        G_cell(:,:,2) = q_yex(:,:,2)**2 / max(1d-12,abs(q_yex(:,:,1)))
        G_cell(:,:,3) = q_yex(:,:,2) * q_yex(:,:,3) / max(1d-12,abs(q_yex(:,:,1)))
        G_cell(:,:,2) = G_cell(:,:,2) + P1_yex

        call LF_ghost(F_cell,G_cell,boundary,ny,nx)

        call LFHE_Res(q_xex, F_cell, Fcell_mid, c_P1x, 2, ny, nx, alpha_F,debug)
        call LFHE_Res(q_yex, G_cell, Gcell_mid, c_P1y, 1, ny, nx, alpha_G,debug)
        alpha_FG(i_OD)    =   max(alpha_F,alpha_G)

        ! boundary value
        ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
        do i = 1,ny,1
            do j = 1,nx-1,1
                if ((boundary(i,j)==2).OR.(boundary(i,j+1)==2)) then
                    Fcell_mid(i,j,1) = q_in(2,i_OD)
                    Fcell_mid(i,j,2) = q_in(2,i_OD)**2 / max(1d-12,q_in(1,i_OD)) + (c0 ** 2) * q_in(1,i_OD)
                    Fcell_mid(i,j,3) = q_in(2,i_OD) * q_in(3,i_OD) / max(1d-12,q_in(1,i_OD))
                end if
                if ((boundary(i,j)==0).OR.(boundary(i,j+1)==0)) then
                    Fcell_mid(i,j,1) = 0d0
                    Fcell_mid(i,j,3) = 0d0
                end if
            end do
        end do
        do i = 2,ny,1
            do j = 1,nx,1
                if ((boundary(i-1,j)==2).OR.(boundary(i,j)==2)) then
                    Gcell_mid(i,j,1) = q_in(3,i_OD)
                    Gcell_mid(i,j,2) = q_in(2,i_OD) * q_in(3,i_OD) / max(1d-12,q_in(1,i_OD))
                    Gcell_mid(i,j,3) = q_in(3,i_OD)**2 / max(1d-12,q_in(1,i_OD)) + (c0 ** 2) * q_in(1,i_OD)
                end if
                if ((boundary(i-1,j)==0).OR.(boundary(i,j)==0)) then
                    Gcell_mid(i,j,1) = 0d0
                    Gcell_mid(i,j,2) = 0d0
                end if
            end do
        end do

        ! boudnary out limit

        if (bound_out(2)>0) then
            if ((i_OD == 5).OR.(i_OD == 6)) then
                do i = 1,ny,1
                    do j = 1,nx-1,1
                        if ((boundary(i,j)==1).OR.(boundary(i,j+1)==1)) then
                            Fcell_mid(i,j,1) = min(0.1d0,max(-0.1d0,Fcell_mid(i,j,1)))
                        end if
                    end do
                end do
                do i = 2,ny,1
                    do j = 1,nx,1
                        if ((boundary(i-1,j)==1).OR.(boundary(i,j)==1)) then
                            Gcell_mid(i,j,1) = min(0.1d0,max(-0.1d0,Gcell_mid(i,j,1)))
                        end if
                    end do
                end do
            end if
        end if

        do i = 1,3,1
            q2(:,:,i,i_OD)  =   q1(:,:,i,i_OD) - 1d0/h * tstep * (Fcell_mid(:,:,i) - eoshift(Fcell_mid(:,:,i),-1,0d0,2)+&
            &Gcell_mid(:,:,i) - eoshift(Gcell_mid(:,:,i), 1,0d0,1))
        end do

        q2(:,:,2,i_OD)  =   q2(:,:,2,i_OD) + tstep * (q1(:,:,1,i_OD) * Ve(:,:,1,i_OD) - q1(:,:,2,i_OD)) / tao(:,:,i_OD)
        q2(:,:,3,i_OD)  =   q2(:,:,3,i_OD) + tstep * (q1(:,:,1,i_OD) * Ve(:,:,2,i_OD) - q1(:,:,3,i_OD)) / tao(:,:,i_OD)
    end do

    alpha_LF = maxval(alpha_FG,1)

    ! Derive density gradient
    boundary = bound_all(:,:,1)
    do j = 1,ny,1
        do i = 1,nx,1
            ! central
            ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
            if (boundary(j,i)>2) then
                if ((boundary(j,i-1)>2) .AND. (boundary(j,i+1)>2)) then
                    rho_x = (rho_com(j,i+1) - rho_com(j,i-1)) / h / 2d0
                else
                    if (((boundary(j,i-1)<=2)) .AND. ((boundary(j,i+1)>2))) then
                        rho_x = (rho_com(j,i+1) - rho_com(j,i)) / h
                    end if
                    if (((boundary(j,i+1)<=2)) .AND. ((boundary(j,i-1)>2))) then
                        rho_x = (rho_com(j,i) - rho_com(j,i-1)) / h
                    end if
                end if
                
                if (((boundary(j-1,i)>2)) .AND. ((boundary(j+1,i)>2))) then
                    rho_y = (rho_com(j-1,i) - rho_com(j+1,i)) / h / 2d0
                else
                    if (((boundary(j-1,i)<=2)) .AND. ((boundary(j+1,i)>2))) then
                        rho_y = (rho_com(j,i) - rho_com(j+1,i)) / h
                    end if
                    if (((boundary(j+1,i)<=2)) .AND. ((boundary(j-1,i)>2))) then
                        rho_y = (rho_com(j-1,i) - rho_com(j,i)) / h
                    end if
                end if
            else
                rho_x = 0
                rho_y = 0
            end if

            if (sqrt(rho_x**2+rho_y**2)>1d-12) then
                rho_x = rho_x / sqrt(rho_x**2+rho_y**2)
                rho_y = rho_y / sqrt(rho_x**2+rho_y**2)
            else
                rho_x = 0
                rho_y = 0
            end if

            angle = acos(ve_com(j,i,1)*rho_x+ve_com(j,i,2)*rho_y)
            if (angle<=(3.1415926535897d0/2d0)) then
                region(j,i) = 1
            else
                region(j,i) = 2
            end if
        end do
    end do
    
    alpha       =   min(1d0,max(0d0, (rho_com - rho_0) / (rho_m - rho_0)))
    conflict    =   sqrt(ve_com(:,:,1)**2+ve_com(:,:,2)**2)
    pforce      =   p2capacity * sqrt(max(0d0, rho_com - rho_0)) * conflict * maxval(panic,3);

    call F90_fsmgod_P2(pforce,alpha,region,bound_all(:,:,1),P2,p_x,p_y,h,ny,nx)
    F_P2(:,:,1) = p_x
    F_P2(:,:,2) = p_y
    ! p_x = 0d0
    ! p_y = 0d0

    do i_OD = 1,n_OD,1
        q2(:,:,2,i_OD)  =   q2(:,:,2,i_OD) - tstep * p_x * rho_ex(:,:,1,i_OD)/max(1d-12,abs(rho_com)) / m_ave
        q2(:,:,3,i_OD)  =   q2(:,:,3,i_OD) - tstep * p_y * rho_ex(:,:,1,i_OD)/max(1d-12,abs(rho_com)) / m_ave 

        q2(:,:,3,i_OD)  =   q2(:,:,3,i_OD) - tstep * slope * alpha * rho_ex(:,:,1,i_OD) / tao(:,:,i_OD)
    end do

    q2((/1:3,ny-2:ny/),:,:,:) = 0d0
    q2(:,(/1:3,nx-2:nx/),:,:) = 0d0

    do i_OD = 1,n_OD,1
        do i = 1,ny,1
            do j = 1,nx,1
                if (q2(i,j,1,i_OD)<=1d-6) then
                    q2(i,j,2,i_OD) = q2(i,j,1,i_OD)*Ve(i,j,1,i_OD)
                    q2(i,j,3,i_OD) = q2(i,j,1,i_OD)*Ve(i,j,2,i_OD)
                end if
            end do
        end do
        q2(:,:,1,i_OD) = merge(0d0,q2(:,:,1,i_OD),(bound_all(:,:,i_OD)==0).OR.(bound_all(:,:,i_OD)==1))
        q2(:,:,2,i_OD) = merge(0d0,q2(:,:,2,i_OD),(bound_all(:,:,i_OD)==0).OR.(bound_all(:,:,i_OD)==1))
        q2(:,:,3,i_OD) = merge(0d0,q2(:,:,3,i_OD),(bound_all(:,:,i_OD)==0).OR.(bound_all(:,:,i_OD)==1))
    end do

    return
end subroutine FDM_euler

subroutine F90_fsmgod(cost,boundary,potential,p_x,p_y,h,ny,nx,nit)
    implicit none
    integer, parameter          ::  acc = 8
    real(kind=acc), parameter   ::  sigma = 1d-12,  p_max = 1d6
    integer, Intent(In)         ::  ny, nx,boundary(ny,nx)
    real(kind=acc), Intent(In)  ::  cost(ny,nx),h
    real(kind=acc), Intent(Out) ::  potential(ny,nx),p_x(ny,nx),p_y(ny,nx)
    integer                     ::  i, j, n, nit_c = 1, Gsit(4,2,2), GSi(4,2),nit, bj_fsm(ny,nx)
    real(kind=acc)              ::  p_old(ny,nx), pot_x, pot_y, p_t, p_temp

    Gsit    = reshape(  (/1 , nx, nx, 1 ,&
                            ny, ny, 1 , 1 ,&
                            nx, 1 , 1 , nx,&
                            1 , 1 , ny, ny/),(/4,2,2/))
    GSi     = reshape(  (/ 1, -1, -1,  1,&
                        -1, -1,  1,  1/),(/4,2/))
    p_x = 0d0
    p_y = 0d0
    bj_fsm = boundary
    potential = p_max
    potential = merge(0d0,potential,boundary==1)
    bj_fsm = merge(0,boundary,((boundary==0).OR.(boundary==1)))
    p_old = 0d0
    nit = 0
    ! Obtain potential by Fast Godunov
    do nit = 1,10000,1
        if ( ( sum(abs(potential - p_old)) ) >= sigma ) then
            p_old = potential
            do n = 1,4,1
                do j = GSit(n,2,1), GSit(n,2,2), GSi(n,2)
                    do i = GSit(n,1,1), GSit(n,1,2), GSi(n,1)

                        ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
                        if (bj_fsm(j,i) > 0) then
                                
                            pot_x = min(potential(j,max(1,i-1)), potential(j,min(nx,i+1)))
                            pot_y = min(potential(max(1,j-1),i), potential(min(ny,j+1),i))
                            if (abs(pot_x-pot_y) > (cost(j,i)*h)) then
                                p_temp = (min(pot_x, pot_y) + cost(j,i)*h)
                            else 
                                p_temp = (pot_x+pot_y + sqrt(2d0 * (cost(j,i)*h)**2 - (pot_x-pot_y)**2)) / 2d0
                            end if
                            
                            if (potential(j,max(1,i-1)) >= potential(j,min(nx,i+1))) then
                                p_x(j,i) = (potential(j,min(nx,i+1)) - potential(j,i)) / h
                            else
                                p_x(j,i) = (potential(j,i) - potential(j,max(1,i-1))) / h
                            end if
                            if (potential(min(ny,j+1),i) >= potential(max(1,j-1),i)) then
                                p_y(j,i) = (potential(max(1,j-1),i) - potential(j,i)) / h
                            else
                                p_y(j,i) = (potential(j,i) - potential(min(ny,j+1),i)) / h
                            end if
                            potential(j,i) = min(potential(j,i),p_temp)
                        end if
                    end do                  
                end do
            end do
        else
            exit
        end if
        
    end	do

    return
end subroutine F90_fsmgod

subroutine LF_ghost(flow_x,flow_y,boundary,ny,nx)
    implicit none
    integer, parameter          ::  acc = 8
    integer, Intent(In)         ::  ny, nx, boundary(ny,nx)
    real(kind=acc)              ::  flow_x(ny,nx,3),flow_y(ny,nx,3)
    integer                     ::  i,j
    ! ghost points
    ! 0 is bound_H
        do i = 1,ny,1
            do j = 1,nx,1
                if (boundary(i,j)>2) then
                    if (boundary(i,j+1)==0) then
                        flow_x(i,j+1,1)           =   - flow_x(i,j,1)
                        flow_x(i,j+1,2)           =   flow_x(i,j,2)
                        flow_x(i,j+1,3)           =   - flow_x(i,j,3)
                    end if
                    if (boundary(i,j-1)==0) then
                        flow_x(i,j-1,1)           =   - flow_x(i,j,1)
                        flow_x(i,j-1,2)           =   flow_x(i,j,2)
                        flow_x(i,j-1,3)           =   - flow_x(i,j,3)
                    end if
                    if (boundary(i+1,j)==0) then
                        flow_y(i+1,j,1)           =   - flow_y(i,j,1)
                        flow_y(i+1,j,2)           =   flow_y(i,j,2)
                        flow_y(i+1,j,3)           =   - flow_y(i,j,3)
                    end if
                    if (boundary(i-1,j)==0) then
                        flow_y(i-1,j,1)           =   - flow_y(i,j,1)
                        flow_y(i-1,j,2)           =   flow_y(i,j,2)
                        flow_y(i-1,j,3)           =   - flow_y(i,j,3)
                    end if
                end if
            end do
        end do
    
    return
end subroutine LF_ghost

subroutine LFHE_Res(qcell, fcell, fcell_mid,c_P1, dim, ny, nx, alpha_m, debug)
    implicit none
    integer, parameter          ::  acc     =   8
    real(kind=acc), parameter   ::  rho_0   = 6d0, rho_1 = 7d0, c0 = 0.6d0
    integer, Intent(In)         ::  ny, nx, dim
    real(kind=acc)              ::  fcell(ny,nx,3),qcell(ny,nx,3), rho_mid(ny,nx), fcell_mid(ny,nx,3), c_P1(ny,nx), fcell_p1(ny,nx,3),&
    &qcell_p1(ny,nx,3), vx(ny,nx), vx_mid(ny,nx), vy_mid(ny,nx),&
    &alpha_H(ny,nx), alpha_E(ny,nx), ME_mat(ny,nx), M_E(ny,nx),&
    &K_p(3,1), L_p1(3,3), R_p1(3,3), F_s(3,1),Q_s(3,1), F_p1_s(3,1),&
    &Q_p1_s(3,1), F_p1_s_m(3,1), F_s_p(3,1), F_p1_s_m_E(3), alpha_temp(3), &
    &F_s_p_E(3), temp(ny,nx), alpha_m, c_P1_mid, debug(10,10)
    integer                     ::  i, j
    ! Set Variables

    debug   =   0d0
    alpha_m =   0d0
    if (dim == 2)then
        do i = 1,3,1
            fcell_p1(:,:,i) = eoshift(fcell(:,:,i),1,fcell(:,nx,i),2)
            qcell_p1(:,:,i) = eoshift(qcell(:,:,i),1,qcell(:,nx,i),2)
        end do
        rho_mid = (qcell(:,:,1)+qcell_p1(:,:,1))/2d0
        vx = qcell(:,:,2)/max(1d-12,abs(qcell(:,:,1)))

        vx_mid =( qcell(:,:,2)+qcell_p1(:,:,2) ) / max(1d-12,abs(rho_mid)) /2d0
        vy_mid =( qcell(:,:,3)+qcell_p1(:,:,3) ) / max(1d-12,abs(rho_mid)) /2d0
        ! LLF-H
        alpha_H = c_P1 + abs(vx);
        alpha_H = merge(0d0, alpha_H, qcell(:,:,1) > rho_1)
        alpha_H = max(alpha_H,eoshift(alpha_H,1,alpha_H(:,nx),2))
        ! LF-E
        ME_mat = -2d0*abs(vx) + 2d0*sqrt(max(0d0,vx**2-c_P1**2)) + 0.1d0
        ME_mat = merge(0d0,ME_mat,qcell(:,:,1)<=rho_1)
        M_E = spread(maxval(ME_mat,2),2,nx)
        alpha_E = (sqrt(max(0d0,M_E**2-4d0 * M_E * vx+4d0*c_P1**2))-M_E+abs(2d0*vx)) / 2d0
        alpha_E = merge(0d0,alpha_E,qcell(:,:,1)<=rho_1)
        alpha_E = spread(maxval(alpha_E,2),2,nx)
    else if (dim == 1) then
        do i = 1,3,1
            fcell_p1(:,:,i) = eoshift(fcell(:,:,i),-1,fcell(1,:,i),1)
            qcell_p1(:,:,i) = eoshift(qcell(:,:,i),-1,qcell(1,:,i),1)
        end do
        rho_mid = (qcell(:,:,1)+qcell_p1(:,:,1))/2d0
        vx = qcell(:,:,2)/max(1d-12,abs(qcell(:,:,1)))

        vx_mid =( qcell(:,:,2)+qcell_p1(:,:,2) ) / max(1d-12,abs(rho_mid)) /2d0
        vy_mid =( qcell(:,:,3)+qcell_p1(:,:,3) ) / max(1d-12,abs(rho_mid)) /2d0
        
        ! LLF-H
        alpha_H = c_P1 + abs(vx)
        alpha_H = merge(0d0,alpha_H,qcell(:,:,1) > rho_1)
        alpha_H = max(alpha_H,eoshift(alpha_H,-1,alpha_H(1,:),1))
        ! LF-E
        ME_mat = -2d0 * abs(vx) + 2d0 * sqrt(max(0d0,vx**2-c_P1**2)) + 0.1d0
        ME_mat = merge(0d0,ME_mat,qcell(:,:,1)<=rho_1)
        M_E = spread(maxval(ME_mat,1),1,ny)
        alpha_E = (sqrt(max(0d0,M_E**2-4d0 * M_E * vx + 4d0 * c_P1**2))-M_E+abs(2d0*vx)) / 2d0
        alpha_E = merge(0d0,alpha_E,qcell(:,:,1)<=rho_1)
        alpha_E= spread(maxval(alpha_E,1),1,ny)
    end if

    do i = 1,ny,1
        do j = 1,nx,1
            if (rho_mid(i,j)/=0) then
                if (rho_mid(i,j)<=rho_1) then

                    if (rho_mid(i,j)<=rho_0) then
                        c_P1_mid   = c0
                    else
                        c_P1_mid   = c0 / 2d0
                    end if
                    L_p1 = transpose( reshape( (/ (vx_mid(i,j)+c_P1_mid)/2d0/c_P1_mid, -1d0/2d0/c_P1_mid, 0d0, -vy_mid(i,j),0d0,1d0, -(vx_mid(i,j)-c_P1_mid)/2d0/c_P1_mid,1d0/2d0/c_P1_mid,0d0 /), shape(L_p1)) )

                    F_s    = matmul(L_p1,reshape(fcell(i,j,:),(/3,1/) ))
                    Q_s    = matmul(L_p1,reshape(qcell(i,j,:),(/3,1/) ))
                    F_p1_s = matmul(L_p1,reshape(fcell_p1(i,j,:),(/3,1/) ))
                    Q_p1_s = matmul(L_p1,reshape(qcell_p1(i,j,:),(/3,1/) )) 

                    F_p1_s_m = 1d0 / 2d0 * ( F_p1_s - alpha_H(i,j) * Q_p1_s )
                    F_s_p    = 1d0 / 2d0 * ( F_s    + alpha_H(i,j) * Q_s    )
                    K_p = F_s_p + F_p1_s_m

                    R_p1 = transpose( reshape( (/ 1d0,0d0,1d0,vx_mid(i,j)-c_P1_mid,0d0,vx_mid(i,j)+c_P1_mid, vy_mid(i,j),1d0,vy_mid(i,j)/), shape(R_p1)) )
                    fcell_mid(i,j,:) = reshape(matmul(R_p1,K_p), shape(fcell_mid(i,j,:)))
                    alpha_m = max(alpha_m,alpha_H(i,j))

                else

                    alpha_temp = (/alpha_E(i,j)+M_E(i,j), alpha_E(i,j) , alpha_E(i,j) /)
                    F_p1_s_m_E = 1d0 / 2d0 * ( fcell_p1(i,j,:)    - alpha_temp * qcell_p1(i,j,:)  )
                    F_s_p_E    = 1d0 / 2d0 * ( fcell_p1(i,j,:)    + alpha_temp * qcell_p1(i,j,:)  )
                    fcell_mid(i,j,:) = F_s_p_E + F_p1_s_m_E;
                    alpha_m = max(alpha_m,alpha_E(i,j)+M_E(i,j))

                end if
            else
                fcell_mid(i,j,:) = 0
            end if        
        end do 
    end do

    if (dim == 1) then
        temp = fcell_mid(:,:,2)
        fcell_mid(:,:,2) = fcell_mid(:,:,3)
        fcell_mid(:,:,3) = temp
    end if                
    
    return
end subroutine LFHE_Res

subroutine F90_fsmgod_P2(pforce,alpha,region,boundary,potential,p_x,p_y,h,ny,nx)
    implicit none
    integer, parameter          ::  acc = 8
    real(kind=acc), parameter   ::  sigma = 1d-12,  p_max = 1d6
    integer, Intent(In)         ::  ny, nx, boundary(ny,nx), region(ny,nx)
    real(kind=acc), Intent(In)  ::  pforce(ny,nx),alpha(ny,nx),h
    real(kind=acc), Intent(Out) ::  potential(ny,nx),p_x(ny,nx),p_y(ny,nx)
    integer                     ::  i, j, n, nit_c = 1, Gsit(4,2,2), GSi(4,2),nit, bj_fsm(ny,nx)
    real(kind=acc)              ::  pforce2(ny,nx),p_alpha(ny,nx),p_old(ny,nx), pot_x, pot_y, p_t, p_temp

    Gsit    = reshape(  (/1 , nx, nx, 1 ,&
                            ny, ny, 1 , 1 ,&
                            nx, 1 , 1 , nx,&
                            1 , 1 , ny, ny/),(/4,2,2/))
    GSi     = reshape(  (/ 1, -1, -1,  1,&
                        -1, -1,  1,  1/),(/4,2/))
    ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
    pforce2 = pforce/max(alpha,1d-12)
    p_x = 0d0
    p_y = 0d0

    potential   = p_max
    p_alpha     = p_max

    potential   = merge(0d0,potential,alpha<=1d-12)
    potential   = merge(p_max,potential,(boundary==0).OR.(boundary==1))
    p_alpha     = merge(0d0,p_alpha,alpha<=1d-12)
    p_alpha     = merge(p_max,p_alpha,(boundary==0).OR.(boundary==1))
    bj_fsm      = merge(0,boundary,alpha<=1d-12)

    p_old = 0d0
    nit = 0
    ! Obtain potential by Fast Godunov
    do nit = 1,10000,1
        if ( ( sum(abs(potential - p_old)) ) >= sigma ) then
            p_old = potential
            do n = 1,4,1
                do j = GSit(n,2,1), GSit(n,2,2), GSi(n,2)
                    do i = GSit(n,1,1), GSit(n,1,2), GSi(n,1)

                        ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
                        if (bj_fsm(j,i) > 0) then
                            if (region(j,i) == 1) then
                                pot_x = min(potential(j,max(1,i-1)), potential(j,min(nx,i+1)))
                                pot_y = min(potential(max(1,j-1),i), potential(min(ny,j+1),i))
                                if (abs(pot_x-pot_y) > (pforce(j,i)*h)) then
                                    p_temp = (min(pot_x, pot_y) + pforce(j,i)*h)
                                else 
                                    p_temp = (pot_x+pot_y + sqrt(2d0*(pforce(j,i)*h)**2-(pot_x-pot_y)**2)) / 2d0
                                end if
                                if (alpha(j,i) <= 1d-12) then
                                    p_alpha(j,i)    =   0d0
                                    potential(j,i)  =   0d0
                                else
                                    potential(j,i)  =   min(potential(j,i),p_temp)
                                    p_alpha(j,i)    =   potential(j,i)/alpha(j,i);
                                end if
                            else
                                pot_x = min(p_alpha(j,max(1,i-1)), p_alpha(j,min(nx,i+1)))
                                pot_y = min(p_alpha(max(1,j-1),i), p_alpha(min(ny,j+1),i))
                                if (abs(pot_x-pot_y) > (pforce2(j,i)*h)) then
                                    p_temp = (min(pot_x, pot_y) + pforce2(j,i)*h)
                                else 
                                    p_temp = (pot_x+pot_y + sqrt(2d0*(pforce2(j,i)*h)**2-(pot_x-pot_y)**2)) / 2d0
                                end if
                                p_alpha(j,i)    = min(p_alpha(j,i),p_temp)
                                potential(j,i)  = p_alpha(j,i)*alpha(j,i)
                            end if
                            
                            ! if (potential(j,max(1,i-1)) >= potential(j,min(nx,i+1))) then
                            !     p_x(j,i) = (potential(j,min(nx,i+1)) - potential(j,i)) / h
                            ! else
                            !     p_x(j,i) = (potential(j,i) - potential(j,max(1,i-1))) / h
                            ! end if
                            ! if (potential(min(ny,j+1),i) >= potential(max(1,j-1),i)) then
                            !     p_y(j,i) = (potential(max(1,j-1),i) - potential(j,i)) / h
                            ! else
                            !     p_y(j,i) = (potential(j,i) - potential(min(ny,j+1),i)) / h
                            ! end if

                        end if
                    end do                  
                end do
            end do
        else
        ! central
            do j = 1,ny,1
                do i = 1,nx,1
                    if (bj_fsm(j,i) > 0) then
                        if (((boundary(j,i-1)>2)) .AND. ((boundary(j,i+1)>2))) then
                            p_x(j,i) = (potential(j,i+1) - potential(j,i-1)) / h / 2
                        else
                            if (((boundary(j,i-1)<=2)) .AND. ((boundary(j,i+1)>2))) then
                                p_x(j,i) = (potential(j,i+1) - potential(j,i)) / h
                            end if
                            if (((boundary(j,i-1)>2)) .AND. ((boundary(j,i+1)<=2))) then
                                p_x(j,i) = (potential(j,i) - potential(j,i-1)) / h
                            end if
                        end if
                        
                        if (((boundary(j-1,i)>2)) .AND. ((boundary(j+1,i)>2))) then
                            p_y(j,i) = (potential(j-1,i) - potential(j+1,i)) / h / 2
                        else
                            if (((boundary(j-1,i)<=2)) .AND. ((boundary(j+1,i)>2))) then
                                p_y(j,i) = (potential(j,i) - potential(j+1,i)) / h
                            end if
                            if (((boundary(j-1,i)>2)) .AND. ((boundary(j+1,i)<=2))) then
                                p_y(j,i) = (potential(j-1,i) - potential(j,i)) / h
                            end if
                        end if
                    end if
                end do 
            end do
            exit
        end if
        
    end	do

    return
end subroutine F90_fsmgod_P2