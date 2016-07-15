
% Wigner3j.m by David Terr, Raytheon, 6-17-04

% Compute the Wigner 3j symbol using the Racah formula [1].
% EFK 20130326: modified to use pretabulated factorials for speed

function wigner = Wigner3j(j1,j2,j3,m1,m2,m3,factsave)

% error checking
if ( 2*j1 ~= floor(2*j1) || 2*j2 ~= floor(2*j2) || 2*j3 ~= floor(2*j3) ...
        || 2*m1 ~= floor(2*m1) || 2*m2 ~= floor(2*m2) || 2*m3 ~= floor(2*m3) )
    error('All arguments must be integers or half-integers.');
    return;
end

if ( j1 - m1 ~= floor ( j1 - m1 ) )
    error('2*j1 and 2*m1 must have the same parity');
    return;
end

if ( j2 - m2 ~= floor ( j2 - m2 ) )
    error('2*j2 and 2*m2 must have the same parity');
    return;
end

if ( j3 - m3 ~= floor ( j3 - m3 ) )
    error('2*j3 and 2*m3 must have the same parity');
    return;
end

if j3 > j1 + j2 || j3 < abs(j1 - j2)
    error('j3 is out of bounds.');
    return;
end

if abs(m1) > j1
    error('m1 is out of bounds.');
    return;
end

if abs(m2) > j2
    error('m2 is out of bounds.');
    return;
end

if abs(m3) > j3
    error('m3 is out of bounds.');
    return;
end

if m1+m2+m3 ~= 0
    error('m1+m2+m3 is not 0')
    return;
end

t1 = j2 - m1 - j3;
t2 = j1 + m2 - j3;
t3 = j1 + j2 - j3;
t4 = j1 - m1;
t5 = j2 + m2;

tmin = max( 0, max( t1, t2 ) );
tmax = min( t3, min( t4, t5 ) );

wigner = 0;
if nargin > 6
    % use pretabulated factorials
    for t = tmin:tmax
        wigner = wigner + (-1)^t / ( fact(t) * fact(t-t1) * fact(t-t2) ...
            * fact(t3-t) * fact(t4-t) * fact(t5-t) );
    end
    
    wigner = wigner * (-1)^(j1-j2-m3) ...
        * sqrt( fact(j1+j2-j3) * fact(j1-j2+j3) * fact(-j1+j2+j3) / fact(j1+j2+j3+1)...
        * fact(j1+m1) * fact(j1-m1) * fact(j2+m2) * fact(j2-m2) * fact(j3+m3) * fact(j3-m3) );
else
    % calculate factorials from scratch
    for t = tmin:tmax
        wigner = wigner + (-1)^t / ( factorial(t) * factorial(t-t1) * factorial(t-t2) ...
            * factorial(t3-t) * factorial(t4-t) * factorial(t5-t) );
    end
    
    wigner = wigner * (-1)^(j1-j2-m3) ...
        * sqrt( factorial(j1+j2-j3) * factorial(j1-j2+j3) * factorial(-j1+j2+j3) / factorial(j1+j2+j3+1)...
        * factorial(j1+m1) * factorial(j1-m1) * factorial(j2+m2) * factorial(j2-m2) * factorial(j3+m3) * factorial(j3-m3) );
end

% function for including 0 index in pretabulated factorials
    function res = fact(x)
        if (x==0)
            res = 1;
        else
            res = factsave(x);
        end
    end
end
% Reference: Wigner 3j-Symbol entry of Eric Weinstein's Mathworld: http://mathworld.wolfram.com/Wigner3j-Symbol.html