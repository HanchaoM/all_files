function Buu = FBuu(M,k,x,QQ,t,phi)
%syms y e xi tmin kpr gg q qp po pmag cth theta sth sthl cthl cthpr sthpr M M2 TAU X Y F1 F2
[row,col] =size(x);
if(isa(x,'double') && isa(phi,'double'))
    Buu = ones(size(x));
else
    Buu = sym(ones(size(x)));
end

 for i = 1:row
    for j = 1:col
        M2 = M*M;
        y = QQ / ( 2 * M * k * x(i,j) ); %% From eq. (23) where gamma is substituted from eq (12c): gamma = 2*M*x(i,j)/Q
        gg = 4* M2 * x(i,j) * x(i,j) / QQ; %% This is gamma^2 [from eq. (12c)] 
        e = ( 1 - y - ( y * y * (gg / 4) ) ) / ( 1 - y + (y * y / 2) + ( y * y * (gg / 4) ) ); %% epsilon eq. (32)
        xi =  x(i,j) * ( ( 1 + t / ( 2 * QQ ) ) / ( 2 - x(i,j) + x(i,j) * t / QQ ) ); %% skewness parameter eq. (12b) note: there is a minus sign on the write up that shouldn't be there
        tmin = ( QQ * ( 1 - sqrt( 1 + gg ) + gg / 2 ) ) / ( x(i,j) * ( 1 - sqrt( 1 + gg ) + gg / ( 2 * x(i,j) ) ) ); %% minimum t eq. (29)
        kpr = k * ( 1 - y ); %% k' from eq. (23)
        qpr = t / 2 / M + k - kpr; %% q' from eq. bellow to eq. (25) that has no numbering. Here nu = k - k' = k * y
        po = M - t / 2 / M; %% This is p'_0 from eq. (28b)
        pmag = sqrt( ( -t ) * ( 1 - t / 4 / M / M ) ); %% p' magnitude from eq. (28b)
        % cth = ; %% This is cos(theta) eq. (26)
        sthl = sqrt( gg ) / sqrt( 1 + gg ) * ( sqrt ( 1 - y - y * y * gg / 4) ); %% sin(theta_l) from eq. (22a)
        cthl = -1 / sqrt( 1 + gg ) * ( 1 + y * gg / 2 ) ; %% cos(theta_l) from eq. (22a)
        cth = -1 / sqrt( 1 + gg ) * ( 1 + gg / 2 * ( 1+ t / QQ ) / ( 1 + x * t / QQ ) );
        theta = acos(cth); %% theta angle from eq. (26)


        %% momenta vectors defined on eq. (21)
        K = [k, k * sthl, 0,  k * cthl];
        KP = [kpr, K(2), 0, k * ( cthl + y * sqrt( 1 + gg ))];
        Q = K - KP;
        QP = [qpr, qpr * sin(theta) * cos( phi(i,j) ), qpr * sin(theta) * sin( phi(i,j)), qpr * cos(theta)];
        P = [M, 0, 0, 0];
        D = Q - QP; %% delta vector eq. (12a)
        PP = P + D; %% p' from eq. (21)


        PPK = P + K;
        PPP = (P + PP) /2;  %%P vector eq. (12a)

        %% 4-vectors products
        KD = mulp(K, D);
        KPD = mulp(KP, D);
        QD = mulp(Q, D);
        QPD = mulp(QP, D);
        KQP = mulp(K, QP);
        KKP = mulp(K, KP);
        KPQP = mulp(KP, QP);
        KPPP = mulp(K, PPP);
        KPPPP = mulp(KP, PPP);
        PPPQ = mulp(PPP, Q);
        PPPQP = mulp(PPP, QP);

        %% Auu and Buu Twist 2 interference coefficients
        KK_T = 0.5 * ( e / ( 1 - e ) ) * QQ;

        KKP_T = KK_T;

        KQP_T = ( QQ / ( sqrt( gg ) * sqrt( 1 + gg ) ) ) * sqrt ( (0.5 * e) / ( 1 - e ) ) * ( 1 + x(i,j) * t / QQ ) * sin(theta) * cos( phi(i,j) );

        KD_T = - KQP_T;

        DD_T = ( 1 - xi * xi ) * ( tmin - t );

        % Interference coefficients given on eq. (241a,b,c)--------------------

        Buu(i,j) = 2 * xi / ( KQP * KPQP) * ( ( QQ + t ) * ( 2* KK_T * ( KD + KPD ) + KQP_T * ( QD - KQP - KPQP + 2*KKP  ) + 2* KQP * KPD - 2* KPQP * KD ) + ( QQ - t + 4* KD ) * ( ( KK_T - 2* KKP ) * QPD - KKP * DD_T - 2* KD_T * KQP ) );
    end
end

end


function [mulp] = mulp(a,b)

mulp = a(1)*b(1);
for i=2:length(a)
             mulp = mulp - a(i)*b(i);
end

end