function gsm_vs_som(rm, w, pr, N)
%gsm_vs_som Porovnani Gauss–Seidel method oproti Successive over-relaxation

disp('===================================================')
disp('====NOVY START=NOVY START=NOVY START=NOVY START====')
disp('===================================================')

%START Zadani ==================================
if isempty(rm) == 1
    rm = 5;
end
if isempty(w) == 1
    w = 1.5;
end
if isempty(pr) == 1
    pr = 1e-05;
end
if isempty(N) == 1
    N = 50;
end
e = ones(rm,1);
A = full(spdiags([-e 2*e -e],-1:1,rm,rm));
b = [1:rm]';
disp('Zadani:')
if rm < 6
    disp('A=')
    disp(A)
    disp('b=')
    disp(b)
else
    disp('Rozmer matice A a vektoru b je vetsi nez 5 pro prehlednost vypisuji jen rozmer.')
    disp(rm)
end
disp('N=')
disp(N)
disp('pr=')
disp(pr)
disp('w= (omega)')
disp(w)

disp('!Spoustim [u_sor, it_sor] = sor(A,b,w,pr,N)')
tic
[u_sor, it_sor] = sor(A,b,w,pr,N);
cas_sor=toc;
disp('!Koncim [u_sor, it_sor] = sor(A,b,w,pr,N)')

disp('!Spoustim [u_gs, it_gs] = gauss_seidel(A, b, pr, N)')
tic
[u_gs, it_gs] = gauss_seidel(A, b, pr, N); 
cas_gs=toc;
disp('!Koncim [u_gs, it_gs] = gauss_seidel(A, b, pr, N)')
%KONEC Zadani ==================================
format long;
disp('Vysledky:')
if rm < 6
    disp('u_sor:');disp(u_sor(1:rm,1:1));
    disp('u_gs:');disp(u_gs(1:rm,1:1));
    disp('Rozdil vysledku:')
    for i= 1:rm
        disp(abs(u_sor(i:i,1:1)-u_gs(i:i,1:1)));
    end
else
    disp('(omezuji vypis pouze na prvnich 5 clenu vysledneho vektoru)')
    disp('u_sor:');disp(u_sor(1:5,1:1));
    disp('u_gs:');disp(u_gs(1:5,1:1));
    disp('Rozdil vysledku:')
    for i= 1:5
        disp(abs(u_sor(i:i,1:1)-u_gs(i:i,1:1)));
    end
end
disp('Casy behu:')
disp('cas_sor: '); disp(cas_sor);
disp('cas_gs: '); disp(cas_gs);
disp('Pocet iteraci:')
disp('it_sor: '); disp(it_sor);
disp('it_gs: '); disp(it_gs);




%gauss_seidel===========================================================
    function [u, it] = gauss_seidel(A, b, pr, N)

        %rozdeleni matice A na L, U a D
        D = diag(diag(A));
        L = tril(-A,-1);
        U = triu(-A,1);

        %matice a konstantni vektor pouzivane pro iteraci
        Tg = inv(D-L)*U; 
        cg = inv(D-L)*b;

        k = 1;
        x = b;			%startovni vektor

        while k <= N
           x(:,k+1) = Tg*x(:,k) + cg;
           if norm(x(:,k+1)-x(:,k)) < pr
              disp('INFO:Pozadovane presnosti bylo dosazeno drive nez maxima iteraci.')
              it=k;
              u=x(:,k+1);
              break
           end
           k = k+1;
        end
        if k > N
           disp('INFO:Dosazeno maxima iteraci - vysledek bude pravdepodobne nepresny.')
           it=k-1;
           u=x(:,k);
        end
    end %function gauss_seidel

%SOR====================================================================
    function [u,it] = sor(A,F,w,pr,N)
        D = diag(diag(A));
        L = -tril(A,-1);
        U = -triu(A,1);
        DL = D - w*L;
        UE = (1-w)*D + w*U;
        FE = w*F;
        u0=F;

        for it = 1:N
            rhv = UE*u0 + FE;
            u = DL\rhv;
            if abs(u-u0)/abs(u0) < pr
                disp('INFO:Pozadovane presnosti bylo dosazeno drive nez maxima iteraci.')
                return
            else
                u0 = u;
            end
        end
        disp('INFO:Dosazeno maxima iteraci - vysledek bude pravdepodobne nepresny.')
    end %function sor

disp('===========================================================')
disp('===KONEC=KONEC=KONEC=KONEC=KONEC=KONEC=KONEC=KONEC=KONEC===')
disp('===========================================================')
disp('=INFO O VSTUPECH:==========================================')
disp('=gsm_vs_som(rozmer_matice, omega, presnost, pocet iteraci)=')
disp('=napriklad:================================================')
disp('=   gsm_vs_som(15, 1.5, 1e-05, 500)   =====================')
disp('===========================================================')
end

