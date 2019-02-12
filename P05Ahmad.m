function [nLevel,height,nItem,wOut,hOut,x,y] = P05Ahmad(w,h,L)
%Correct code
nLevel = 0;
n = size(w,1);
l = L;
I = 0;
j = 1;
r = 1;
sum = 0;


while(norm(w,2) ~= 0)
    L = l;
    u = w~=0;
    for i = 1:n
        z = u(i);
        if(z==1)
            wOut1(i) = w(i);
            hOut1(i) = h(i);
            w(i) = 0;
            h(i) = 0;
            break;
        end
    end
    if(norm(w,2) ~= 0)
        L = L - wOut1(1);
        f = [-w];
        A = [w]';
        b = L;
        u = bintprog(f,A,b);
        
        for i = 2:n
            z = u(i);
            
            if(z == 1)
                wOut1(i) = w(i);
                hOut1(i) = h(i);
                w(i) = 0;
                h(i) = 0;
            else
                wOut1(i) = 0;
                hOut1(i) = 0;
            end
        end
        
    end
   
   
    
    if(r==1)
        v1 = nonzeros(wOut1)';
        [p,k] = size(v1);
        z1 = nonzeros(hOut1)';
        height(j) = z1(1);
        x(1) = 0;
        [p,k] = size(v1);
        sum = 0;
        %k-1
        for i = 1:k-1
            sum = v1(i)+sum;
            x(i+1) = sum;
        end
        sum = 0;
        x1 = x;
        clear x;
        y1 = zeros(1,k);
        
        
        
        wOut1 = nonzeros(wOut1)';
        hOut1 = nonzeros(hOut1)';
        
        [j1,k1] = size(x1);
        [j1,m1] = size(hOut1);

        L = l;
        u = w~= 0;
        
        
        for i = 1:n
           
            z = u(i);
            
            if(z == 1)
            
                c = L - w(i);
                if(c < 0)
                    break;
                end
                
                sum = 0;
                sum2 = 0;
                k = 0;
                m = 0;
                
                for p = 1:k1
                    sum = x1(p);
                    k = k+1;
                    if(sum >= c)
                        break;
                    end
                end

                
                for p = 1:k1
                    sum2 = x1(p);
                    m = m+1;
                    if(sum2 >= c)
                        break;
                    end
                end
                
                if(sum > c)
                    k = k-1;
                    m = m-1;
                end
                
                
                sum = 0;
                
                if((w(i) <= L - x1(k))&&(h(i)<=height(j)-hOut1(m)))
                    wOut2(i) = w(i);
                    hOut2(i) = h(i);
                    x2(i) = L - w(i);
                    y2(i) = height(j) - h(i);
                    w(i) = 0;
                    h(i) = 0;
                    L = c;  
                end
            else
                wOut2(i) = 0;
                hOut2(i) = 0;
                x2(i) = 0;
                y2(i) = 0;
            end

        end
        

        

        
        v2 = nonzeros(wOut2)';
        z2 = nonzeros(hOut2)';
        x2 = nonzeros(x2)';
        y2 = nonzeros(y2)';
        nItem(j) = nnz(v1)+nnz(v2);
        v1 = [v1,v2];
        z1 = [z1,z2];
        x1 = [x1,x2];
        y1 = [y1,y2];
        clear wOut1;
        clear hOut1;
        clear wOut2;
        clear hOut2;
        
        
        
        
    else
        
        v3 = nonzeros(wOut1)';
        [p,k] = size(v3);
        Q = nnz(v3);
        z3 = nonzeros(hOut1)';
        height(j) = z3(1);
        x(1) = 0;
        [p,k] = size(v3);
        sum = 0;
        for i = 1:k-1
            sum = v3(i)+sum;
            x(i+1) = sum;
        end
        sum = 0;
        x3 = x;
        clear x;
        y3 = zeros(1,k);
        v1 = [v1,v3];
        z1 = [z1,z3];
        x1 = [x1,x3];
        y1 = [y1,y3];
        
        
        
        wOut1 = nonzeros(wOut1)';
        hOut1 = nonzeros(hOut1)';
        
        [j1,k1] = size(x3);
        [j1,m1] = size(hOut1);
        
        
        L = l;
        u = w~= 0;
        
        for i = 1:n
           
            z = u(i);
            if(z == 1)

                c = L - w(i);
                if(c < 0)
                    break;
                end
                
                sum = 0;
                sum2 = 0;
                k = 0;
                m = 0;
                
                for p = 1:k1
                    sum = x3(p);
                    k = k+1;
                    if(sum >= c)
                        break;
                    end
                end
                
                for p = 1:k1
                    sum2 = x3(p);
                    m = m+1;
                    if(sum2 >= c)
                        break;
                    end
                end
                
                if(sum > c)
                    k = k-1;
                    m = m-1;
                end

                sum = 0;
                
                if((w(i) <= L - x3(k))&&(h(i)<=height(j)-hOut1(m)))
                    wOut2(i) = w(i);
                    hOut2(i) = h(i);
                    x4(i) = L - w(i);
                    y4(i) = height(j) - h(i);
                    w(i) = 0;
                    h(i) = 0;
                    L = c;
                    
                end
                
            else
                wOut2(i) = 0;
                hOut2(i) = 0;
                x4(i) = 0;
                y4(i) = 0;
            end
            
        end
        

        
        
        
        v4 = nonzeros(wOut2)';
        U = nnz(v4);
        z4 = nonzeros(hOut2)';
        x4 = nonzeros(x4)';
        y4 = nonzeros(y4)';
        nItem(j) = Q+U;
        v1 = [v1,v4];
        z1 = [z1,z4];
        x1 = [x1,x4];
        y1 = [y1,y4];
        clear wOut1;
        clear hOut1;
        clear wOut2;
        clear hOut2;

        
    end
    r = r+1;
    j = j+1;
    nLevel = nLevel + 1;
    clear wOut;
    clear hOut;
    
   
end

wOut = v1;
hOut = z1;
x = x1;
y = y1;
%end

   close all
   ii = 0;
   for I = 1:nLevel
       figure(I)
       for J = 1:nItem(I)
           ii = ii + 1;
           rectangle('Position',[x(ii),y(ii),wOut(ii),hOut(ii)],'facecolor',[0.5 0.5 0.5])
       end
   end
end